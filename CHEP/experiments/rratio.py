import argparse
import math
from symbolica import NumericalIntegrator, Sample, RandomNumberGenerator
from CHEP.phase_space_generators.phase_space_generators import FlatPhaseSpace
from CHEP.utils.vectors import LorentzVector
from CHEP.matrix_elements.madgraph.processes.all_processes import Matrix_3_epem_ddxg_no_z
from CHEP.matrix_elements.madgraph.model.parameters import ModelParameters
from CHEP.utils import logger

GEV_TO_PB = 0.389379338e9

COSTHETA_CUT = 0.1
GLUON_ENERGY_CUT = 50.0


def pass_cuts(event):
    # debug
    # logger.info(event)
    (pq, pqx, pg) = (event[2], event[3], event[4])
    cosThetaqg = pq.space().dot(pg.space())/(abs(pq.space())*abs(pg.space()))
    cosThetaqxg = pqx.space().dot(pg.space())/(abs(pqx.space())*abs(pg.space()))

    if 1-cosThetaqg < COSTHETA_CUT or 1-cosThetaqxg < COSTHETA_CUT:
        return False

    if pg[0] < GLUON_ENERGY_CUT:
        return False
    return True


def integrand(ps_generator, model, E_cm: float, process, samples_batch: list[Sample]) -> list[float]:

    evaluations: list[float] = []
    for sample in samples_batch:
        ps_point, jacobian = ps_generator.generateKinematics(E_cm, sample.c)

        if pass_cuts(ps_point):
            matrix_element_evaluation = process.smatrix(ps_point, model)
            initial_state_flux = 1.0/(8*math.pi**2*E_cm**2)
            evaluations.append(
                matrix_element_evaluation
                * jacobian
                * initial_state_flux
                * GEV_TO_PB
            )
        else:
            evaluations.append(0.0)
    return evaluations


def rratio(args: argparse.Namespace):

    model = ModelParameters(None)

    process = Matrix_3_epem_ddxg_no_z()
    external_masses = process.get_external_masses(model)

    E_cm = 1000.0  # 1 GeV collision

    ps_generator = FlatPhaseSpace(
        external_masses[0], external_masses[1],
        beam_Es=(E_cm/2., E_cm/2.),
        # We do not consider PDF for e+ e- > l+ l- at fixed-order
        beam_types=(0, 0)
    )

    n_dimensions = ps_generator.nDimPhaseSpace()

    N_CORES = 1  # Parallelization not implemented yet
    DISCRETE_LEARNING_RATE = 0.15
    CONTINUOUS_LEARNING_RATE = 0.15

    parallel_rngs = [RandomNumberGenerator(
        seed=args.seed, stream_id=i_core) for i_core in range(N_CORES)]
    integrator = NumericalIntegrator.continuous(n_dimensions)
    for i_iteration in range(args.n_iterations):
        samples = integrator.sample(
            args.n_points_per_iteration, parallel_rngs[0])
        res = integrand(ps_generator, model, E_cm, process, samples)
        integrator.add_training_samples(samples, res)
        avg, err, chi_sq = integrator.update(
            discrete_learning_rate=DISCRETE_LEARNING_RATE,
            continuous_learning_rate=CONTINUOUS_LEARNING_RATE)  # type: ignore # nopep8
        logger.info(
            'Iteration {}: {:.6} +- {:.6}, chi={:.6}'.format(i_iteration, avg, err, chi_sq))


def rratio_analyze_events(args: argparse.Namespace):
    from CHEP.utils import lhe_parser
    logger.info('TODO ANALYZE MADGRAPH EVENTS')
