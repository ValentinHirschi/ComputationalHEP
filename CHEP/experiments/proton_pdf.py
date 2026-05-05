import argparse
import math
import os
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from symbolica import Expression

from CHEP.utils import CHEPException, logger


# LHAPDF serves the official PDF grids from this public directory.  The
# experiment downloads the tarball only when the set is not already present in
# the local PDFsets directory.  Keeping the download logic here makes the
# exercise reproducible on a fresh checkout while still letting students inspect
# the installed grids as plain text files.
LHAPDF_CURRENT_URL = "https://lhapdfsets.web.cern.ch/lhapdfsets/current"

# This is the NNPDF4.0 NLO fixed-four-flavour set with alpha_s(M_Z)=0.118.
# The companion PDF+alpha_s variation set is NNPDF40_nlo_nf_4_pdfas, but for a
# first pedagogical run a single central set gives cleaner plots and sum rules.
DEFAULT_PDF_SET = "NNPDF40_nlo_as_01180_nf_4"

# LHAPDF uses PDG particle IDs for flavours.  We expose names in the CLI because
# they are friendlier in a teaching setting.
FLAVOUR_TO_PDG = {
    "g": 21,
    "gluon": 21,
    "d": 1,
    "u": 2,
    "s": 3,
    "c": 4,
    "b": 5,
    "t": 6,
    "dbar": -1,
    "ubar": -2,
    "sbar": -3,
    "cbar": -4,
    "bbar": -5,
    "tbar": -6,
}
PDG_TO_LABEL = {pdg: name for name, pdg in FLAVOUR_TO_PDG.items() if len(name) <= 4}
PDG_TO_LABEL[21] = "g"

# Colour factors for SU(3).  They enter the leading-order DGLAP splitting
# kernels below.  These are deliberately left visible instead of hidden in a
# black-box package.
CF = 4.0 / 3.0
CA = 3.0
TR = 0.5


@dataclass
class PdfSumRules:
    """Container for the numerical sum rules checked at one factorisation scale."""

    q: float
    momentum: float
    u_valence: float
    d_valence: float
    s_valence: float
    c_valence: float

    @property
    def baryon_number(self) -> float:
        # A proton contains three valence quarks, so (u_v+d_v+s_v+c_v)/3 = 1.
        return (self.u_valence + self.d_valence + self.s_valence + self.c_valence) / 3.0

    @property
    def isospin_3(self) -> float:
        # In the (p,n) isospin doublet, I_3 = (u_v-d_v)/2 for the proton.
        return 0.5 * (self.u_valence - self.d_valence)


@dataclass
class DglapComparison:
    """Small summary of the toy DGLAP evolution compared with LHAPDF."""

    start_q: float
    target_q: float
    flavour: int
    mean_relative_difference: float
    max_relative_difference: float
    plot_file: Path


def _parse_float_list(raw: str) -> list[float]:
    values = []
    for item in raw.split(","):
        item = item.strip()
        if item:
            values.append(float(item))
    if not values:
        raise CHEPException("Expected at least one comma-separated floating point value")
    return values


def _parse_flavour_list(raw: str) -> list[int]:
    flavours = []
    for item in raw.split(","):
        name = item.strip().lower()
        if not name:
            continue
        if name.lstrip("-").isdigit():
            flavours.append(int(name))
            continue
        if name not in FLAVOUR_TO_PDG:
            known = ", ".join(sorted(FLAVOUR_TO_PDG))
            raise CHEPException(f"Unknown PDF flavour '{item}'. Known names are: {known}")
        flavours.append(FLAVOUR_TO_PDG[name])
    if not flavours:
        raise CHEPException("Expected at least one PDF flavour")
    return flavours


def _parse_single_flavour(raw: str) -> int:
    name = raw.strip().lower()
    if name.lstrip("-").isdigit():
        return int(name)
    if name in FLAVOUR_TO_PDG:
        return FLAVOUR_TO_PDG[name]
    known = ", ".join(sorted(FLAVOUR_TO_PDG))
    raise CHEPException(f"Unknown PDF flavour '{raw}'. Known names are: {known}")


def _x_grid(x_min: float, x_max: float, n_points: int) -> np.ndarray:
    """Build an x grid with more points at small x, where PDFs vary rapidly."""

    if not 0.0 < x_min < x_max < 1.0:
        raise CHEPException("Need 0 < x_min < x_max < 1 for the PDF x grid")
    if n_points < 20:
        raise CHEPException("Use at least 20 x points for a meaningful PDF plot")

    # A pure linear grid wastes points near x=1, while a pure logarithmic grid
    # can undersample the valence peak.  This hybrid grid gives enough detail in
    # both regions without making the exercise slow.
    split = min(0.1, 0.5 * (x_min + x_max))
    n_log = max(10, int(0.65 * n_points))
    n_lin = max(10, n_points - n_log)
    log_part = np.geomspace(x_min, split, n_log, endpoint=False)
    lin_part = np.linspace(split, x_max, n_lin)
    return np.unique(np.concatenate((log_part, lin_part)))


def _safe_filename(value: float) -> str:
    return f"{value:g}".replace(".", "p").replace("-", "m")


def _import_lhapdf(lhapdf_python_dir: str | None):
    """Import LHAPDF, optionally from the MadGraph/HEPTools installation."""

    if lhapdf_python_dir:
        lhapdf_python_dir = os.path.abspath(lhapdf_python_dir)
        if os.path.isdir(lhapdf_python_dir) and lhapdf_python_dir not in sys.path:
            sys.path.insert(0, lhapdf_python_dir)

    try:
        import lhapdf  # type: ignore
    except ImportError as exc:
        raise CHEPException(
            "Could not import the Python LHAPDF module. Install LHAPDF with "
            "Python bindings, or pass --lhapdf_python_dir pointing at the "
            "site-packages directory containing lhapdf.py."
        ) from exc

    return lhapdf


def _prepend_lhapdf_path(lhapdf, pdfsets_dir: Path) -> None:
    """Tell LHAPDF to look in this checkout before its global PDF directories."""

    pdfsets_dir.mkdir(parents=True, exist_ok=True)
    if hasattr(lhapdf, "pathsPrepend"):
        lhapdf.pathsPrepend(str(pdfsets_dir))
    elif hasattr(lhapdf, "setPaths"):
        # Older bindings expose only setPaths.  Prepending manually keeps the
        # behaviour equivalent to pathsPrepend for this experiment.
        lhapdf.setPaths([str(pdfsets_dir)] + list(lhapdf.paths()))


def _download_pdf_set(pdf_set: str, pdfsets_dir: Path) -> None:
    """Download and unpack one LHAPDF set tarball into the local PDFsets folder."""

    if (pdfsets_dir / pdf_set / f"{pdf_set}.info").is_file():
        logger.info("PDF set %s is already installed in %s", pdf_set, pdfsets_dir)
        return

    url = f"{LHAPDF_CURRENT_URL}/{pdf_set}.tar.gz"
    logger.info("Downloading %s", url)

    pdfsets_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmpdir:
        archive_path = Path(tmpdir) / f"{pdf_set}.tar.gz"
        try:
            with urllib.request.urlopen(url) as response:
                with archive_path.open("wb") as archive:
                    shutil.copyfileobj(response, archive)
        except OSError as exc:
            raise CHEPException(
                f"Could not download {pdf_set} from {url}. Check the network "
                "connection, or download the tarball manually into PDFsets."
            ) from exc

        # Never trust tarball paths blindly.  LHAPDF tarballs should contain a
        # single top-level directory named after the set, but the check below
        # prevents accidental extraction outside pdfsets_dir.
        with tarfile.open(archive_path, "r:gz") as tar:
            target_root = pdfsets_dir.resolve()
            for member in tar.getmembers():
                member_target = (pdfsets_dir / member.name).resolve()
                if target_root not in member_target.parents and member_target != target_root:
                    raise CHEPException(f"Refusing unsafe path in PDF tarball: {member.name}")
            tar.extractall(pdfsets_dir)

    if not (pdfsets_dir / pdf_set / f"{pdf_set}.info").is_file():
        raise CHEPException(
            f"Downloaded {pdf_set}, but did not find the expected .info file "
            f"under {pdfsets_dir / pdf_set}"
        )


def _load_pdf(args: argparse.Namespace):
    """Import LHAPDF, ensure the requested set is installed, and return member 0."""

    lhapdf = _import_lhapdf(args.lhapdf_python_dir)
    pdfsets_dir = Path(args.lhapdf_pdfsets_dir).expanduser().resolve()
    _prepend_lhapdf_path(lhapdf, pdfsets_dir)

    try:
        pdf = lhapdf.mkPDF(args.pdf_set, args.pdf_member)
    except Exception as first_error:
        if not args.download_pdf:
            raise CHEPException(
                f"LHAPDF could not load {args.pdf_set}/{args.pdf_member}; "
                "rerun with --download_pdf or install the set manually."
            ) from first_error

        _download_pdf_set(args.pdf_set, pdfsets_dir)
        _prepend_lhapdf_path(lhapdf, pdfsets_dir)
        try:
            pdf = lhapdf.mkPDF(args.pdf_set, args.pdf_member)
        except Exception as second_error:
            raise CHEPException(
                f"LHAPDF still could not load {args.pdf_set}/{args.pdf_member} "
                "after downloading the grid."
            ) from second_error

    logger.info("Loaded LHAPDF set %s member %s", args.pdf_set, args.pdf_member)
    return lhapdf, pdf


def _xfx(pdf, pdg: int, x: float, q: float) -> float:
    """LHAPDF's basic query: x times f_i(x,Q) for one flavour."""

    return float(pdf.xfxQ(int(pdg), float(x), float(q)))


def plot_flavours(
    pdf,
    output_dir: Path,
    x_values: np.ndarray,
    q_values: list[float],
    flavours: list[int],
    pdf_set: str,
) -> list[Path]:
    """Plot x f_i(x,Q) for all requested flavours at each scale Q."""

    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    files = []
    for q in q_values:
        fig, ax = plt.subplots(figsize=(8.0, 5.5))

        for pdg in flavours:
            y_values = [_xfx(pdf, pdg, x, q) for x in x_values]
            ax.plot(x_values, y_values, label=PDG_TO_LABEL.get(pdg, str(pdg)))

        ax.set_xscale("log")
        ax.set_xlabel("x")
        ax.set_ylabel(r"$x f_i(x,\mu_F^2)$")
        ax.set_title(f"{pdf_set}, Q = {q:g} GeV")
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(ncol=3, fontsize=9)
        fig.tight_layout()

        output_file = output_dir / f"proton_pdf_xfx_Q_{_safe_filename(q)}.png"
        fig.savefig(output_file, dpi=160)
        plt.close(fig)
        files.append(output_file)
        logger.info("Wrote %s", output_file)

    return files


def _integrate_trapezoid(x_values: np.ndarray, y_values: np.ndarray) -> float:
    # numpy.trapezoid was introduced after numpy.trapz.  Use it when available,
    # but keep compatibility with older environments used in tutorials.
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y_values, x_values))
    return float(np.trapz(y_values, x_values))


def check_sum_rules(
    pdf,
    x_values: np.ndarray,
    q_values: list[float],
    active_quarks: list[int],
) -> list[PdfSumRules]:
    """Numerically check momentum, valence, baryon number, and isospin rules."""

    results = []
    for q in q_values:
        # Momentum sum rule:
        #   integral_0^1 dx x [g + sum_q(q + qbar)] = 1.
        # Since LHAPDF gives x*f directly, the integrand is just xfxQ.
        momentum_density = np.array([_xfx(pdf, 21, x, q) for x in x_values])
        for flavour in active_quarks:
            momentum_density += np.array([
                _xfx(pdf, flavour, x, q) + _xfx(pdf, -flavour, x, q)
                for x in x_values
            ])
        momentum = _integrate_trapezoid(x_values, momentum_density)

        # Valence sum rule:
        #   integral_0^1 dx [q(x)-qbar(x)] = N_q.
        # LHAPDF gives x*q, so we divide by x to recover q itself.
        valence = {}
        for flavour in active_quarks:
            density = np.array([
                (_xfx(pdf, flavour, x, q) - _xfx(pdf, -flavour, x, q)) / x
                for x in x_values
            ])
            valence[flavour] = _integrate_trapezoid(x_values, density)

        result = PdfSumRules(
            q=q,
            momentum=momentum,
            u_valence=valence.get(2, 0.0),
            d_valence=valence.get(1, 0.0),
            s_valence=valence.get(3, 0.0),
            c_valence=valence.get(4, 0.0),
        )
        results.append(result)

        logger.info("Sum rules at Q = %.6g GeV:", q)
        logger.info("  momentum = %.6f (target 1)", result.momentum)
        logger.info("  u_v      = %.6f (target 2)", result.u_valence)
        logger.info("  d_v      = %.6f (target 1)", result.d_valence)
        logger.info("  s_v      = %.6f (target 0)", result.s_valence)
        logger.info("  c_v      = %.6f (target 0)", result.c_valence)
        logger.info("  baryon   = %.6f (target 1)", result.baryon_number)
        logger.info("  I_3      = %.6f (target 0.5)", result.isospin_3)

    return results


def write_sum_rule_table(output_dir: Path, results: list[PdfSumRules]) -> Path:
    """Write a small CSV table so students can compare grids or x ranges."""

    output_file = output_dir / "proton_pdf_sum_rules.csv"
    with output_file.open("w", encoding="utf-8") as table:
        table.write(
            "Q,momentum,u_valence,d_valence,s_valence,c_valence,"
            "baryon_number,isospin_3\n"
        )
        for result in results:
            table.write(
                f"{result.q:.12g},{result.momentum:.12g},"
                f"{result.u_valence:.12g},{result.d_valence:.12g},"
                f"{result.s_valence:.12g},{result.c_valence:.12g},"
                f"{result.baryon_number:.12g},{result.isospin_3:.12g}\n"
            )
    logger.info("Wrote %s", output_file)
    return output_file


@dataclass
class KernelExpression:
    """A readable Symbolica expression plus a fast numerical evaluator."""

    expression: Expression
    values: Callable[[np.ndarray], np.ndarray]


@dataclass
class DglapKernels:
    """Symbolica expressions for the regular parts of the LO kernels."""

    z: Expression
    pqq_regular: KernelExpression
    pqg: KernelExpression
    pgq: KernelExpression
    pgg_plus_regular: KernelExpression
    pgg_regular: KernelExpression


def _build_lo_kernels() -> DglapKernels:
    """Build the splitting kernels with Symbolica expressions.

    DGLAP evolution is an integro-differential equation.  For teaching, it is
    useful to keep the kernels as explicit symbolic formulae before turning them
    into numerical functions.  Symbolica is not doing magic here; it simply gives
    us readable formula objects that we then evaluate in the quadrature.
    """

    z = Expression.symbol("z")

    # P_qq has a plus distribution:
    #   C_F [(1+z^2)/(1-z)]_+ + 3/2 C_F delta(1-z).
    # We store only the singular regular function; the plus prescription and
    # delta term are implemented explicitly in _plus_convolution.
    pqq_regular = CF * (1 + z**2) / (1 - z)

    # A gluon splitting to q qbar.  Each quark or anti-quark receives this
    # contribution from the gluon.
    pqg = TR * (z**2 + (1 - z) ** 2)

    # A quark radiating a gluon.
    pgq = CF * (1 + (1 - z) ** 2) / z

    # We use the common decomposition
    #   P_gg = 2 C_A [ z/(1-z)_+ + (1-z)/z + z(1-z)]
    #          + beta_0/2 delta(1-z).
    pgg_plus_regular = 2 * CA * z / (1 - z)
    pgg_regular = 2 * CA * ((1 - z) / z + z * (1 - z))

    return DglapKernels(
        z=z,
        pqq_regular=KernelExpression(
            pqq_regular,
            lambda zz: CF * (1.0 + zz * zz) / (1.0 - zz),
        ),
        pqg=KernelExpression(
            pqg,
            lambda zz: TR * (zz * zz + (1.0 - zz) * (1.0 - zz)),
        ),
        pgq=KernelExpression(
            pgq,
            lambda zz: CF * (1.0 + (1.0 - zz) * (1.0 - zz)) / zz,
        ),
        pgg_plus_regular=KernelExpression(
            pgg_plus_regular,
            lambda zz: 2.0 * CA * zz / (1.0 - zz),
        ),
        pgg_regular=KernelExpression(
            pgg_regular,
            lambda zz: 2.0 * CA * ((1.0 - zz) / zz + zz * (1.0 - zz)),
        ),
    )


def _kernel_eval(kernel: KernelExpression, z_values: np.ndarray) -> np.ndarray:
    """Evaluate one DGLAP kernel on a numpy array."""

    # The Symbolica expression is kept in KernelExpression.expression for
    # readability and inspection.  The vectorized callable is used for the
    # numerical work so that the default exercise finishes quickly.
    return kernel.values(z_values)


def _interpolate_xf(x_grid: np.ndarray, xf_values: np.ndarray, y_values: np.ndarray) -> np.ndarray:
    """Interpolate x*f(x) on a grid and force the PDF to vanish at x=1."""

    x_with_endpoint = np.concatenate((x_grid, np.array([1.0])))
    xf_with_endpoint = np.concatenate((xf_values, np.array([0.0])))
    return np.interp(y_values, x_with_endpoint, xf_with_endpoint)


def _regular_convolution(
    x: float,
    x_grid: np.ndarray,
    xf_values: np.ndarray,
    kernel: KernelExpression,
    n_z: int,
) -> float:
    """Compute integral_x^1 dz P(z) [x/z f(x/z)] for a regular kernel."""

    z_values = np.linspace(x, 1.0 - 1e-8, n_z)
    y_values = x / z_values
    integrand = _kernel_eval(kernel, z_values)
    integrand *= _interpolate_xf(x_grid, xf_values, y_values)
    return _integrate_trapezoid(z_values, integrand)


def _plus_convolution(
    x: float,
    x_grid: np.ndarray,
    xf_values: np.ndarray,
    kernel: KernelExpression,
    integral_0_to_x,
    n_z: int,
) -> float:
    """Compute a plus-distribution convolution.

    The plus prescription is the part that is most often hidden in production
    codes.  For a kernel [g(z)]_+ and a test function h(z)=x/z f(x/z),
    one may write
      int_x^1 dz [g(z)]_+ h(z)
        = int_x^1 dz g(z) [h(z)-h(1)] - h(1) int_0^x dz g(z).
    This form is finite at z=1 because the bracket cancels the 1/(1-z)
    singularity.
    """

    z_values = np.linspace(x, 1.0 - 1e-8, n_z)
    y_values = x / z_values
    h_values = _interpolate_xf(x_grid, xf_values, y_values)
    h_at_one = float(np.interp(x, x_grid, xf_values))
    integrand = _kernel_eval(kernel, z_values)
    integrand *= h_values - h_at_one
    return _integrate_trapezoid(z_values, integrand) - h_at_one * integral_0_to_x(x)


def _pqq_integral_0_to_x(x: float) -> float:
    # Integral_0^x dz (1+z^2)/(1-z) = -2 log(1-x)-x-x^2/2.
    return CF * (-2.0 * math.log(1.0 - x) - x - 0.5 * x * x)


def _pgg_plus_integral_0_to_x(x: float) -> float:
    # Integral_0^x dz z/(1-z) = -log(1-x)-x.
    return 2.0 * CA * (-math.log(1.0 - x) - x)


def _dglap_rhs_lo(
    pdf,
    q: float,
    state: dict[int, np.ndarray],
    x_grid: np.ndarray,
    active_quarks: list[int],
    kernels: DglapKernels,
    n_z: int,
) -> dict[int, np.ndarray]:
    """Right-hand side d(xf)/d ln(Q^2) at leading order."""

    beta0 = (11.0 * CA - 4.0 * TR * len(active_quarks)) / 3.0
    alpha_factor = float(pdf.alphasQ(q)) / (2.0 * math.pi)
    rhs = {pdg: np.zeros_like(values) for pdg, values in state.items()}

    for i_x, x in enumerate(x_grid):
        gluon_xf = state[21]
        p_qg_g = _regular_convolution(
            x, x_grid, gluon_xf, kernels.pqg, n_z
        )

        for flavour in active_quarks:
            for pdg in (flavour, -flavour):
                q_xf = state[pdg]
                p_qq_q = _plus_convolution(
                    x,
                    x_grid,
                    q_xf,
                    kernels.pqq_regular,
                    _pqq_integral_0_to_x,
                    n_z,
                )
                p_qq_q += 1.5 * CF * float(np.interp(x, x_grid, q_xf))
                rhs[pdg][i_x] = alpha_factor * (p_qq_q + p_qg_g)

        p_gg_g = _plus_convolution(
            x,
            x_grid,
            gluon_xf,
            kernels.pgg_plus_regular,
            _pgg_plus_integral_0_to_x,
            n_z,
        )
        p_gg_g += _regular_convolution(
            x, x_grid, gluon_xf, kernels.pgg_regular, n_z
        )
        p_gg_g += 0.5 * beta0 * float(np.interp(x, x_grid, gluon_xf))

        quark_to_gluon = 0.0
        for flavour in active_quarks:
            for pdg in (flavour, -flavour):
                quark_to_gluon += _regular_convolution(
                    x, x_grid, state[pdg], kernels.pgq, n_z
                )

        rhs[21][i_x] = alpha_factor * (p_gg_g + quark_to_gluon)

    return rhs


def _state_add_scaled(
    state: dict[int, np.ndarray],
    increment: dict[int, np.ndarray],
    scale: float,
) -> dict[int, np.ndarray]:
    return {pdg: state[pdg] + scale * increment[pdg] for pdg in state}


def _rk4_step(
    pdf,
    t: float,
    dt: float,
    state: dict[int, np.ndarray],
    x_grid: np.ndarray,
    active_quarks: list[int],
    kernels: DglapKernels,
    n_z: int,
) -> dict[int, np.ndarray]:
    """One Runge-Kutta step in t = ln(Q^2)."""

    def q_from_t(t_value: float) -> float:
        return math.exp(0.5 * t_value)

    k1 = _dglap_rhs_lo(pdf, q_from_t(t), state, x_grid, active_quarks, kernels, n_z)
    k2 = _dglap_rhs_lo(
        pdf,
        q_from_t(t + 0.5 * dt),
        _state_add_scaled(state, k1, 0.5 * dt),
        x_grid,
        active_quarks,
        kernels,
        n_z,
    )
    k3 = _dglap_rhs_lo(
        pdf,
        q_from_t(t + 0.5 * dt),
        _state_add_scaled(state, k2, 0.5 * dt),
        x_grid,
        active_quarks,
        kernels,
        n_z,
    )
    k4 = _dglap_rhs_lo(
        pdf,
        q_from_t(t + dt),
        _state_add_scaled(state, k3, dt),
        x_grid,
        active_quarks,
        kernels,
        n_z,
    )

    return {
        pdg: state[pdg] + dt * (k1[pdg] + 2.0 * k2[pdg] + 2.0 * k3[pdg] + k4[pdg]) / 6.0
        for pdg in state
    }


def compare_dglap_running(
    pdf,
    output_dir: Path,
    pdf_set: str,
    start_q: float,
    target_q: float,
    x_min: float,
    x_max: float,
    n_x: int,
    n_steps: int,
    n_z: int,
    active_quarks: list[int],
    comparison_flavour: int,
) -> DglapComparison:
    """Evolve PDFs with a tiny transparent DGLAP solver and compare to LHAPDF.

    This is intentionally not a replacement for a production evolution code.
    It solves the LO singlet DGLAP equations for x*f(x,Q), using LHAPDF's own
    alpha_s(Q) so that the comparison isolates the PDF evolution.  The scaffold
    is the right place to add the full NLO P_ij^(1) kernels next: the RK4 ODE
    solver, plus-distribution treatment, and LHAPDF comparison are already here.
    """

    if start_q <= 0.0 or target_q <= 0.0 or start_q == target_q:
        raise CHEPException("DGLAP comparison needs positive, distinct start and target Q values")
    if n_steps <= 0:
        raise CHEPException("DGLAP comparison needs at least one RK4 step")

    x_values = _x_grid(x_min, x_max, n_x)
    kernels = _build_lo_kernels()
    state = {21: np.array([_xfx(pdf, 21, x, start_q) for x in x_values])}
    for flavour in active_quarks:
        state[flavour] = np.array([_xfx(pdf, flavour, x, start_q) for x in x_values])
        state[-flavour] = np.array([_xfx(pdf, -flavour, x, start_q) for x in x_values])
    starting_values = state[comparison_flavour].copy()

    uniform_log_state = {pdg: values.copy() for pdg, values in state.items()}
    t = math.log(start_q * start_q)
    target_t = math.log(target_q * target_q)
    dt = (target_t - t) / float(n_steps)

    for _ in range(n_steps):
        uniform_log_state = _rk4_step(
            pdf, t, dt, uniform_log_state, x_values, active_quarks, kernels, n_z
        )
        t += dt

    n_linear_q_steps = 10
    linear_q_state = {pdg: values.copy() for pdg, values in state.items()}
    q_values_for_steps = np.linspace(start_q, target_q, n_linear_q_steps + 1)
    for q_left, q_right in zip(q_values_for_steps[:-1], q_values_for_steps[1:]):
        t_left = math.log(q_left * q_left)
        dt_linear_q = math.log(q_right * q_right) - t_left
        linear_q_state = _rk4_step(
            pdf,
            t_left,
            dt_linear_q,
            linear_q_state,
            x_values,
            active_quarks,
            kernels,
            n_z,
        )

    evolved = uniform_log_state[comparison_flavour]
    evolved_linear_q = linear_q_state[comparison_flavour]
    lhapdf_values = np.array([_xfx(pdf, comparison_flavour, x, target_q) for x in x_values])
    scale = np.maximum(np.abs(lhapdf_values), 1e-10)
    relative_difference = (evolved - lhapdf_values) / scale
    relative_difference_linear_q = (evolved_linear_q - lhapdf_values) / scale
    path_relative_difference = np.abs(evolved - evolved_linear_q) / scale
    log_path_relative_difference = np.log10(np.maximum(path_relative_difference, 1e-18))
    paths_are_bitwise_identical = bool(np.array_equal(evolved, evolved_linear_q))

    os.environ.setdefault("MPLCONFIGDIR", str(output_dir / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(8.0, 9.0),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0, 1.0]},
    )
    axes[0].plot(x_values, starting_values, linestyle=":", label=f"LHAPDF at Q={start_q:g} GeV")
    axes[0].plot(x_values, lhapdf_values, label=f"LHAPDF at Q={target_q:g} GeV")
    axes[0].plot(
        x_values,
        evolved,
        linestyle="--",
        label=f"LO DGLAP ({n_steps} uniform ln Q^2 step(s))",
    )
    axes[0].plot(
        x_values,
        evolved_linear_q,
        linestyle="-.",
        label=(
            "LO DGLAP "
            f"({n_linear_q_steps} linear Q steps, "
            f"Delta Q={abs(target_q - start_q) / n_linear_q_steps:g} GeV)"
        ),
    )
    axes[0].set_xscale("log")
    axes[0].set_ylabel(r"$x f_i(x,\mu_F^2)$")
    axes[0].set_title(
        f"{pdf_set}: Q={start_q:g} GeV to Q={target_q:g} GeV, "
        f"flavour {PDG_TO_LABEL.get(comparison_flavour, comparison_flavour)}"
    )
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend()

    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].plot(x_values, relative_difference, label="uniform ln Q^2 steps")
    axes[1].plot(x_values, relative_difference_linear_q, linestyle="-.", label="linear Q steps")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("x")
    axes[1].set_ylabel("(ODE - LHAPDF) / LHAPDF")
    axes[1].grid(True, which="both", alpha=0.25)
    axes[1].legend()

    axes[2].plot(x_values, log_path_relative_difference)
    axes[2].set_xscale("log")
    axes[2].set_xlabel("x")
    axes[2].set_ylabel(r"$\log_{10}(|A-B|/|\mathrm{LHAPDF}|)$")
    axes[2].grid(True, which="both", alpha=0.25)
    fig.tight_layout()

    output_file = output_dir / "proton_pdf_dglap_lo_vs_lhapdf.png"
    fig.savefig(output_file, dpi=160)
    plt.close(fig)

    comparison = DglapComparison(
        start_q=start_q,
        target_q=target_q,
        flavour=comparison_flavour,
        mean_relative_difference=float(np.mean(np.abs(relative_difference))),
        max_relative_difference=float(np.max(np.abs(relative_difference))),
        plot_file=output_file,
    )
    logger.info(
        "DGLAP toy comparison for flavour %s: mean |rel diff| = %.4g, max = %.4g",
        PDG_TO_LABEL.get(comparison_flavour, str(comparison_flavour)),
        comparison.mean_relative_difference,
        comparison.max_relative_difference,
    )
    logger.info(
        "Uniform-lnQ2 and linear-Q RK4 paths bitwise identical? %s; "
        "max |A-B|/|LHAPDF| = %.4g",
        paths_are_bitwise_identical,
        float(np.max(path_relative_difference)),
    )
    logger.info("Wrote %s", output_file)
    return comparison


def proton_pdf(args: argparse.Namespace):
    """Run the proton PDF experiment."""

    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    _, pdf = _load_pdf(args)

    q_values = _parse_float_list(args.q_values)
    plot_flavours_list = _parse_flavour_list(args.flavours)
    active_quarks = [1, 2, 3, 4]
    x_values = _x_grid(args.x_min, args.x_max, args.x_points)

    plot_files = plot_flavours(
        pdf,
        output_dir,
        x_values,
        q_values,
        plot_flavours_list,
        args.pdf_set,
    )
    sum_rules = check_sum_rules(pdf, x_values, q_values, active_quarks)
    sum_rule_file = write_sum_rule_table(output_dir, sum_rules)

    dglap_comparison = None
    if not args.skip_dglap:
        dglap_comparison = compare_dglap_running(
            pdf,
            output_dir,
            args.pdf_set,
            args.dglap_start_q,
            args.dglap_target_q,
            args.dglap_x_min,
            args.dglap_x_max,
            args.dglap_x_points,
            args.dglap_steps,
            args.dglap_z_points,
            active_quarks,
            _parse_single_flavour(args.dglap_flavour),
        )

    logger.info("Generated %d PDF plot(s)", len(plot_files))
    logger.info("Sum-rule table: %s", sum_rule_file)
    if dglap_comparison is not None:
        logger.info("DGLAP comparison plot: %s", dglap_comparison.plot_file)
