import argparse
import math
from symbolica import NumericalIntegrator, Sample, RandomNumberGenerator
from ..phase_space_generators.phase_space_generators import FlatPhaseSpace
from ..matrix_elements.madgraph.processes.all_processes import Matrix_1_ddx_wpwm_no_za, Matrix_1_uux_wpwm_no_za
from ..matrix_elements.madgraph.model.parameters import ModelParameters
from CHEP.utils import logger

GEV_TO_PB = 0.389379338e9


def integrand(ps_generator, pdf_hook, mu_f, model, E_cm: float, processes, samples_batch: list[Sample]) -> list[float]:

    evaluations: list[float] = []
    for sample in samples_batch:
        ps_point, jacobian, (x1, _), (x2,
                                      _) = ps_generator.get_PS_point(sample.c)
        matrix_element_evaluation = 0
        for (proc, (pdg_left, pdg_right)) in processes:
            pdf_left = pdf_hook.xfxQ(pdg_left, x1, mu_f) / x1
            pdf_right = pdf_hook.xfxQ(pdg_right, x2, mu_f) / x2
            matrix_element_evaluation += pdf_left * \
                pdf_right * proc.smatrix(ps_point, model)

        initial_state_flux = 1.0/(8*math.pi**2*x1*x2*E_cm**2)
        evaluations.append(
            matrix_element_evaluation
            * jacobian
            * initial_state_flux
            * GEV_TO_PB
        )
    return evaluations


def pp_wpwm_fixed_order_LO(args: argparse.Namespace):

    model = ModelParameters(None)

    process_uux = Matrix_1_ddx_wpwm_no_za()
    process_ddx = Matrix_1_uux_wpwm_no_za()

    external_masses = process_uux.get_external_masses(model)

    E_cm = 13000.0

    ps_generator = FlatPhaseSpace(
        external_masses[0], external_masses[1],
        beam_Es=(E_cm/2., E_cm/2.),
        beam_types=(1, 1)
    )

    n_dimensions = ps_generator.nDimPhaseSpace()+2

    N_CORES = 1  # Parallelization not implemented yet
    DISCRETE_LEARNING_RATE = 0.15
    CONTINUOUS_LEARNING_RATE = 0.15

    # Load PDF set
    import sys
    if args.lhpadf_python_dir not in sys.path:
        sys.path.insert(0, args.lhpadf_python_dir)
    import lhapdf
    print(lhapdf.__file__)
    lhapdf.pathsPrepend(args.lhpadf_pdfsets_dir)
    pdf_hook = lhapdf.mkPDF(args.pdf_set, 0)

    mu_f = model.mdl_MZ

    parallel_rngs = [RandomNumberGenerator(
        seed=args.seed, stream_id=i_core) for i_core in range(N_CORES)]
    integrator = NumericalIntegrator.continuous(n_dimensions)
    for i_iteration in range(args.n_iterations):
        samples = integrator.sample(
            args.n_points_per_iteration, parallel_rngs[0])
        res = integrand(ps_generator, pdf_hook, mu_f, model, E_cm,
                        [
                            (process_uux, (2, -2)),
                            (process_uux, (-2, 2)),
                            (process_ddx, (1, -1)),
                            (process_ddx, (-1, 1)),
                        ], samples)
        integrator.add_training_samples(samples, res)
        avg, err, chi_sq = integrator.update(
                discrete_learning_rate=DISCRETE_LEARNING_RATE,
                continuous_learning_rate=CONTINUOUS_LEARNING_RATE)  # type: ignore # nopep8
        logger.info(
            'Iteration {}: {:.6} +- {:.6}, chi={:.6}'.format(i_iteration, avg, err, chi_sq))
