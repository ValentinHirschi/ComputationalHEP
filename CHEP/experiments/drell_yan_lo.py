import argparse
import math
from pathlib import Path

import numpy as np
from symbolica import NumericalIntegrator, RandomNumberGenerator, Sample

from CHEP.experiments.proton_pdf import DEFAULT_PDF_SET, _load_pdf
from CHEP.matrix_elements.madgraph.model.parameters import ModelParameters
from CHEP.phase_space_generators.phase_space_generators import FlatPhaseSpace
from CHEP.utils import CHEPException, logger


GEV_TO_PB = 0.389379338e9


# Electric charges and weak isospin assignments for the quarks we allow in the
# incoming proton.  The experiment is intentionally neutral-current only:
# q qbar -> Z, followed by an optional Z -> l+ l- branching-ratio factor.
QUARK_DATA = {
    1: ("d", -1.0 / 3.0, -0.5),
    2: ("u", +2.0 / 3.0, +0.5),
    3: ("s", -1.0 / 3.0, -0.5),
    4: ("c", +2.0 / 3.0, +0.5),
    5: ("b", -1.0 / 3.0, -0.5),
}

QUARK_NAME_TO_PDG = {
    "d": 1,
    "u": 2,
    "s": 3,
    "c": 4,
    "b": 5,
}

LEPTON_DATA = {
    "e": (-1.0, -0.5),
    "mu": (-1.0, -0.5),
    "tau": (-1.0, -0.5),
}


def _parse_quark_flavours(raw: str) -> list[int]:
    flavours = []
    for item in raw.split(","):
        name = item.strip().lower()
        if not name:
            continue
        if name.isdigit():
            pdg = int(name)
        else:
            if name not in QUARK_NAME_TO_PDG:
                known = ", ".join(sorted(QUARK_NAME_TO_PDG))
                raise CHEPException(f"Unknown quark flavour '{item}'. Known names: {known}")
            pdg = QUARK_NAME_TO_PDG[name]
        if pdg not in QUARK_DATA:
            raise CHEPException(f"Quark PDG id {pdg} is not supported in this LO exercise")
        flavours.append(pdg)

    if not flavours:
        raise CHEPException("At least one quark flavour is needed for Drell-Yan")
    return flavours


def _z_vector_axial_couplings(charge: float, weak_isospin: float, sw2: float) -> tuple[float, float]:
    """Return the standard dimensionless Z vector and axial couplings.

    We use the convention
      Z f fbar vertex = i e/(2 s_w c_w) gamma^mu (v_f - a_f gamma5),
    with v_f = T3_f - 2 Q_f s_w^2 and a_f = T3_f.
    """

    vector = weak_isospin - 2.0 * charge * sw2
    axial = weak_isospin
    return vector, axial


def _tree_level_z_to_lepton_branching_ratio(model: ModelParameters, lepton: str) -> float:
    """Estimate BR(Z -> l+ l-) from the tree-level partial width.

    This keeps the experiment visibly tied to the same vector/axial couplings as
    the production channels while using the model's total Z width for the
    denominator.  QED/QCD corrections and lepton masses are deliberately omitted.
    """

    if lepton not in LEPTON_DATA:
        known = ", ".join(sorted(LEPTON_DATA))
        raise CHEPException(f"Unknown charged lepton '{lepton}'. Known names: {known}")

    charge, weak_isospin = LEPTON_DATA[lepton]
    vector, axial = _z_vector_axial_couplings(charge, weak_isospin, model.mdl_sw2)
    partial_width = (
        model.mdl_Gf
        * model.mdl_MZ**3
        * (vector * vector + axial * axial)
        / (6.0 * math.sqrt(2.0) * math.pi)
    )
    return partial_width / model.mdl_WZ


def _pdf_density(pdf_hook, pdg: int, x: float, mu_f: float) -> float:
    """Convert LHAPDF's x*f_i(x,Q) return value into f_i(x,Q)."""

    return float(pdf_hook.xfxQ(pdg, x, mu_f)) / x


def _z_production_matrix_element_squared(
    model: ModelParameters,
    quark_pdg: int,
    m_ll: float,
) -> float:
    """Spin- and colour-averaged |M|^2 for q qbar -> Z at shat = m_ll^2.

    For massless quarks and an on-shell vector boson,

      |M|^2_bar = (1/N_c) * (1/4 spin average)
                  * sum_polarizations |M(q qbar -> Z)|^2
                = g_Z^2 shat (v_q^2 + a_q^2) / 12.

    This is the only hard matrix element in the experiment.  Everything else is
    the proton structure encoded in the PDFs and the 2 -> 1 phase-space
    Jacobian.
    """

    _, charge, weak_isospin = QUARK_DATA[abs(quark_pdg)]
    vector, axial = _z_vector_axial_couplings(charge, weak_isospin, model.mdl_sw2)
    g_z = model.mdl_ee / (model.mdl_sw * model.mdl_cw)
    return g_z * g_z * m_ll * m_ll * (vector * vector + axial * axial) / 12.0


def _channel_contributions(
    pdf_hook,
    model: ModelParameters,
    x1: float,
    x2: float,
    mu_f: float,
    m_ll: float,
    quark_flavours: list[int],
    branching_ratio: float,
) -> dict[str, float]:
    """Return d sigma_hat factors for every q qbar beam assignment.

    The returned numbers do not yet include the 2 -> 1 phase-space/Bjorken
    Jacobian.  That Jacobian is supplied by FlatPhaseSpace.get_PS_point().
    """

    contributions = {}
    for quark_pdg in quark_flavours:
        quark_name = QUARK_DATA[quark_pdg][0]
        matrix_element_squared = _z_production_matrix_element_squared(model, quark_pdg, m_ll)

        # The one-particle phase space gives pi/shat * |M|^2 delta(shat-m_ll^2).
        # The generator has already integrated the delta function over tau,
        # leaving the remaining 1/S in its Jacobian.
        partonic_prefactor = math.pi * matrix_element_squared / (m_ll * m_ll)
        partonic_prefactor *= branching_ratio

        q_from_beam_1 = _pdf_density(pdf_hook, quark_pdg, x1, mu_f)
        qb_from_beam_2 = _pdf_density(pdf_hook, -quark_pdg, x2, mu_f)
        qb_from_beam_1 = _pdf_density(pdf_hook, -quark_pdg, x1, mu_f)
        q_from_beam_2 = _pdf_density(pdf_hook, quark_pdg, x2, mu_f)

        contributions[f"{quark_name} qbar"] = (
            q_from_beam_1 * qb_from_beam_2 * partonic_prefactor
        )
        contributions[f"{quark_name}bar q"] = (
            qb_from_beam_1 * q_from_beam_2 * partonic_prefactor
        )

    return contributions


def _drell_yan_weight(
    pdf_hook,
    model: ModelParameters,
    ps_generator: FlatPhaseSpace,
    random_variables,
    mu_f: float,
    m_ll: float,
    quark_flavours: list[int],
    branching_ratio: float,
) -> tuple[float, float, dict[str, float]]:
    """Evaluate one Monte Carlo point.

    The phase-space generator supplies x1, x2 and the full Jacobian for
    d tau d y including the on-shell delta function for a 2 -> 1 final state.
    """

    ps_point, jacobian, (x1, _), (x2, _) = ps_generator.get_PS_point(random_variables)
    if ps_point is None:
        return 0.0, 0.0, {}

    y_ll = 0.5 * math.log(x1 / x2)
    contributions = _channel_contributions(
        pdf_hook,
        model,
        x1,
        x2,
        mu_f,
        m_ll,
        quark_flavours,
        branching_ratio,
    )
    weight = sum(contributions.values()) * jacobian * GEV_TO_PB
    return weight, y_ll, {name: value * jacobian * GEV_TO_PB for name, value in contributions.items()}


def integrand(
    pdf_hook,
    model: ModelParameters,
    ps_generator: FlatPhaseSpace,
    mu_f: float,
    m_ll: float,
    quark_flavours: list[int],
    branching_ratio: float,
    samples_batch: list[Sample],
) -> list[float]:
    """Batch integrand passed to Symbolica's numerical integrator."""

    evaluations = []
    for sample in samples_batch:
        weight, _, _ = _drell_yan_weight(
            pdf_hook,
            model,
            ps_generator,
            sample.c,
            mu_f,
            m_ll,
            quark_flavours,
            branching_ratio,
        )
        evaluations.append(weight)
    return evaluations


def _rapidity_range(m_ll: float, e_cm: float) -> tuple[float, float]:
    y_max = math.log(e_cm / m_ll)
    return -y_max, y_max


def write_rapidity_distribution(
    pdf_hook,
    model: ModelParameters,
    ps_generator: FlatPhaseSpace,
    output_dir: Path,
    mu_f: float,
    m_ll: float,
    e_cm: float,
    quark_flavours: list[int],
    branching_ratio: float,
    n_bins: int,
) -> tuple[Path, Path]:
    """Write and plot a deterministic d sigma / d y table.

    This uses the same 2 -> 1 phase-space generator as the Monte Carlo
    integration, evaluated at bin centres.  Dividing by the y range converts the
    generator's d sigma / d u weight, where u is the unit-hypercube y variable,
    back into d sigma / d y.
    """

    os_environ_output = output_dir / ".matplotlib"
    os_environ_output.mkdir(parents=True, exist_ok=True)

    import os
    os.environ.setdefault("MPLCONFIGDIR", str(os_environ_output))
    import matplotlib
    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    y_min, y_max = _rapidity_range(m_ll, e_cm)
    y_range = y_max - y_min

    centers = []
    dsigma_dy = []
    channel_names = []
    channel_rows = []

    for i_bin in range(n_bins):
        u = (i_bin + 0.5) / n_bins
        weight, y_ll, channels = _drell_yan_weight(
            pdf_hook,
            model,
            ps_generator,
            [u],
            mu_f,
            m_ll,
            quark_flavours,
            branching_ratio,
        )
        centers.append(y_ll)
        dsigma_dy.append(weight / y_range)

        if not channel_names:
            channel_names = sorted(channels)
        channel_rows.append([channels.get(name, 0.0) / y_range for name in channel_names])

    table_file = output_dir / "drell_yan_lo_rapidity.csv"
    with table_file.open("w", encoding="utf-8") as table:
        table.write("y,dsigma_dy_pb")
        for name in channel_names:
            table.write(f",{name}_pb")
        table.write("\n")
        for y_ll, value, channel_values in zip(centers, dsigma_dy, channel_rows):
            table.write(f"{y_ll:.12g},{value:.12g}")
            for channel_value in channel_values:
                table.write(f",{channel_value:.12g}")
            table.write("\n")

    plot_file = output_dir / "drell_yan_lo_rapidity.png"
    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    ax.plot(centers, dsigma_dy)
    ax.set_xlabel(r"$y_{\ell\ell}$")
    ax.set_ylabel(r"$d\sigma/dy_{\ell\ell}$ [pb]")
    ax.set_title(f"LO Drell-Yan, m_ll = {m_ll:g} GeV, Q = {mu_f:g} GeV")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(plot_file, dpi=160)
    plt.close(fig)

    logger.info("Wrote %s", table_file)
    logger.info("Wrote %s", plot_file)
    return table_file, plot_file


def drell_yan_lo(args: argparse.Namespace):
    """Run pp -> Z -> l+ l- in the 2 -> 1 narrow-width approximation."""

    model = ModelParameters(None)
    m_ll = args.m_ll if args.m_ll is not None else model.mdl_MZ
    mu_f = args.mu_f if args.mu_f is not None else m_ll
    quark_flavours = _parse_quark_flavours(args.quark_flavours)
    output_dir = Path(args.output_dir).expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if m_ll <= 0.0 or m_ll >= args.e_cm:
        raise CHEPException("Need 0 < m_ll < e_cm for the 2 -> 1 Drell-Yan phase space")

    _, pdf_hook = _load_pdf(args)

    if args.branching_ratio is None:
        branching_ratio = _tree_level_z_to_lepton_branching_ratio(model, args.lepton)
    else:
        branching_ratio = args.branching_ratio
    if branching_ratio < 0.0:
        raise CHEPException("The branching ratio must be non-negative")

    ps_generator = FlatPhaseSpace(
        [0.0, 0.0],
        [m_ll],
        beam_Es=(args.e_cm / 2.0, args.e_cm / 2.0),
        beam_types=(1, 1),
    )
    n_dimensions = len(ps_generator.dimensions)

    logger.info("Running LO Drell-Yan with %s phase-space dimension(s)", n_dimensions)
    logger.info("m_ll = %.6g GeV, e_cm = %.6g GeV, mu_f = %.6g GeV", m_ll, args.e_cm, mu_f)
    logger.info("Using Z -> %s+ %s- branching ratio %.6g", args.lepton, args.lepton, branching_ratio)

    parallel_rngs = [
        RandomNumberGenerator(seed=args.seed, stream_id=i_core)
        for i_core in range(1)
    ]
    integrator = NumericalIntegrator.continuous(n_dimensions)

    avg = err = chi_sq = 0.0
    for i_iteration in range(args.n_iterations):
        samples = integrator.sample(args.n_points_per_iteration, parallel_rngs[0])
        res = integrand(
            pdf_hook,
            model,
            ps_generator,
            mu_f,
            m_ll,
            quark_flavours,
            branching_ratio,
            samples,
        )
        integrator.add_training_samples(samples, res)
        avg, err, chi_sq = integrator.update(
            discrete_learning_rate=0.15,
            continuous_learning_rate=0.15,
        )  # type: ignore
        logger.info(
            "Iteration {}: {:.6} +- {:.6} pb, chi={:.6}".format(
                i_iteration, avg, err, chi_sq
            )
        )

    write_rapidity_distribution(
        pdf_hook,
        model,
        ps_generator,
        output_dir,
        mu_f,
        m_ll,
        args.e_cm,
        quark_flavours,
        branching_ratio,
        args.n_bins,
    )

    logger.info("Integrated LO Drell-Yan cross section: %.6g +- %.6g pb", avg, err)
