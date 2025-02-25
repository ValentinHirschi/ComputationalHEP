
import argparse
import random
import math
import matplotlib.pyplot as plt
from symbolica import NumericalIntegrator, Sample
from CHEP.matrix_elements.epem_lplm.model.parameters import ModelParameters
from CHEP.phase_space_generators.multi_channel_phase_space_generator import SingleChannelPhasespace
from CHEP.phase_space_generators.phase_space_generators import FlatPhaseSpace
from CHEP.utils import logger, ModelInformation


FINAL = True
INITIAL = False

DUMMY = 99


def get_topology(channel_name):

    # A specific sets of s- and t-channels for this example:

    ##############################################################################
    # a) e+(1) e-(2) > mu+(3) mu-(4) a(5) with (p3+p4)^2 sampled as a z-resonance
    ##############################################################################
    topology = None

    if channel_name == '2_to_3_s34':

        topology = (
            # s-channels first:
            [
                {
                    'id': DUMMY,  # Irrelevant
                    'legs': [
                        {
                            'id': -12,
                            'number': 3,
                            'state': FINAL,
                        },
                        {
                            'id': 12,
                            'number': 4,
                            'state': FINAL,
                        },
                        {
                            'id': 23,
                            'number': -1,
                            'state': FINAL,
                        }
                    ]
                },
            ],
            # t-channels then:
            [
                {
                    'id': DUMMY,  # Irrelevant
                    'legs': [
                        {
                            'id': 11,
                            'number': 1,
                            'state': INITIAL,
                        },
                        {
                            'id': 22,
                            'number': 5,
                            'state': FINAL,
                        },
                        {
                            'id': 11,
                            'number': -2,
                            'state': FINAL,
                        }
                    ]
                },
                # The dummy vertex below is just to close the diagram and connect
                # with the number -3 which is to be understood here as the initial state #2.
                {
                    'id': DUMMY,  # Irrelevant
                    'legs': [
                        {
                            'id': 11,
                            'number': -2,
                            'state': FINAL,
                        },
                        {
                            'id': 23,
                            'number': -1,
                            'state': FINAL,
                        },
                        {
                            'id': -11,
                            'number': -3,
                            'state': INITIAL,
                        }
                    ]
                },
            ]
        )

    else:
        print("Channel name '%s' not implemented." % channel_name)

    return topology


class IntegrandForTest():
    """An integrand for this phase-space volume test."""

    def __init__(self, phase_space_generator):
        self.phase_space_generator = phase_space_generator
        self.counter = 0

    def __call__(self, samples: list[Sample]) -> list[float]:
        res = []
        for sample in samples:
            self.counter += 1
            PS_point, wgt, x1, x2 = self.phase_space_generator.get_PS_point(
                sample.c)
            res.append(wgt)
        return res


def sampling_experiment(args: argparse.Namespace):

    model = ModelInformation(ModelParameters(None))
    topology = get_topology('2_to_3_s34')

    random.seed(args.seed)
    # Modify below the mass of the resonance
    # It wil be centered here at 150 GeV
    model.parameters.mdl_MZ = 150.0
    # And will have a width of 60 GeV
    model.parameters.mdl_WZ = 60.0

    E_cm = 1000.0

    my_SCPS_generator = SingleChannelPhasespace([0.]*2, [0.]*3,
                                                beam_Es=(E_cm/2., E_cm/2.), beam_types=(0, 0),
                                                model=model, topology=topology, path=[[0,], []])
    random_variables = [random.random()
                        for _ in range(my_SCPS_generator.nDimPhaseSpace())]
    logger.info("Random variables considered: %s", random_variables)

    logger.info("Considering the following topology:")
    logger.info("-"*10)
    logger.info(my_SCPS_generator.get_topology_string(
        my_SCPS_generator.topology, path_to_print=my_SCPS_generator.path))
    logger.info("-"*10)

    # Now generate a point
    momenta, wgt, _, _ = my_SCPS_generator.get_PS_point(random_variables)

    logger.info("Kinematic configuration generated for random variables %s:" % str(
        random_variables))
    logger.info(momenta)
    logger.info("Jacobian weight: %.16f" % wgt)

    logger.info('')
    logger.info('-'*10)
    logger.info('')

    # Finally integrate
    integrator = NumericalIntegrator.continuous(
        n_dims=my_SCPS_generator.nDimPhaseSpace(),
        n_bins=128,
        train_on_avg=False
    )
    logger.info(
        "Now integrating the phase-space volume with the Single Channel Phase-Space parametrisation")
    (avg, err, chi2) = integrator.integrate(
        integrand=IntegrandForTest(my_SCPS_generator),  # type: ignore # nopep8
        max_n_iter=1,
        min_error=0.01,
        n_samples_per_iter=10000,
        seed=args.seed,
        show_stats=True,
    )
    logger.info(
        'Phase-space volume using single channel phase-space: %.4e +/- %.2e (chi^2 = %.2e)' % (avg, err, chi2))

    my_flat_PS_generator = FlatPhaseSpace(
        [0.]*2, [0.]*3, beam_Es=(E_cm/2., E_cm/2.), beam_types=(0, 0))

    integrator = NumericalIntegrator.continuous(
        n_dims=my_flat_PS_generator.nDimPhaseSpace(),
        n_bins=128,
        train_on_avg=False
    )
    logger.info(
        "Now integrating the phase-space volume with the flat phase-space generator")
    (avg, err, chi2) = integrator.integrate(
        integrand=IntegrandForTest(my_flat_PS_generator),  # type: ignore # nopep8
        max_n_iter=5,
        min_error=0.01,
        n_samples_per_iter=100,
        seed=args.seed,
        show_stats=True,
    )
    logger.info(
        'Phase-space volume using flat phase-space parametrisation: %.4e +/- %.2e (chi^2 = %.2e)' % (avg, err, chi2))
    logger.info('')
    logger.info(
        '>>>>> QUESTION 1: Explain why integrating the PS volume with the SCPS parametrisation is less efficient.')

    logger.info('')
    logger.info('-'*100)
    logger.info('')
    logger.info('Now plot the distributions of the invariant mass (p3+p4)^2')
    n_samples = 10000
    flat_data = []
    scps_data = []
    flat_data_weighted = []
    scps_data_weighted = []
    for i in range(n_samples):
        random_variables = [random.random()
                            for _ in range(my_SCPS_generator.nDimPhaseSpace())]
        flat_momenta, flat_wgt, _, _ = my_flat_PS_generator.get_PS_point(
            random_variables)
        scps_momenta, scps_wgt, _, _ = my_SCPS_generator.get_PS_point(
            random_variables)

        flat_momenta_s34 = math.sqrt((flat_momenta[2]+flat_momenta[3]).square())  # type: ignore # nopep8
        scps_momenta_s34 = math.sqrt(
            (scps_momenta[2]+scps_momenta[3]).square())
        flat_data.append((flat_momenta_s34, 1.0))
        scps_data.append((scps_momenta_s34, 1.0))
        flat_data_weighted.append((flat_momenta_s34, flat_wgt))
        scps_data_weighted.append((scps_momenta_s34, scps_wgt))

    fig, axs = plt.subplots(2, 2)
    values, weights = zip(*flat_data)
    axs[0, 0].hist(values, weights=weights, density=1.0, bins=100)
    axs[0, 0].title.set_text('Flat dw/ds34')
    values, weights = zip(*flat_data_weighted)
    axs[0, 1].hist(values, weights=weights, density=1.0, bins=100)
    axs[0, 1].title.set_text('Flat dw/ds34 weighted')
    values, weights = zip(*scps_data)
    axs[1, 0].hist(values, weights=weights, density=1.0, bins=100)
    axs[1, 0].title.set_text('SCPS dw/ds34')
    values, weights = zip(*scps_data_weighted)
    axs[1, 1].hist(values, weights=weights, density=1.0, bins=100)
    axs[1, 1].title.set_text('SCPS dw/ds34 weighted')

    logger.info('')
    logger.info('>>>>> QUESTION 2: Explain what you see in the four plots, and in particular what is the difference for the weighted histogram.')
    logger.info('')
    plt.show()
