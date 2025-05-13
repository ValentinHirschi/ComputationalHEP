import argparse
import math
from symbolica import NumericalIntegrator, Sample, RandomNumberGenerator
from CHEP.phase_space_generators.phase_space_generators import FlatPhaseSpace
from CHEP.matrix_elements.madgraph_gammastar_processes.processes.all_processes import Matrix_1_z_ddx, Matrix_2_z_ddxg
from CHEP.matrix_elements.madgraph_gammastar_processes.model.parameters import ModelParameters
from CHEP.utils import logger
from CHEP.utils.vectors import LorentzVectorList
from CHEP.utils.lhe_parser import CHEPEventFile, CHEPEvent, CHEPParticle, LegState

GEV_TO_PB = 0.389379338e9

# We are looking at a decay here, so do not consider picobarns but keep GeV.
UNIT_CONVERSION_FACTOR = GEV_TO_PB
UNIT_CONVERSION_FACTOR = 1

# COSTHETA_CUT = 0.3
# GLUON_ENERGY_CUT = 30.0
COSTHETA_CUT = None
GLUON_ENERGY_CUT = None

max_real = 0.


def pass_cuts(event):
    # debug
    # logger.info(event)
    if len(event) == 4:
        (pq, pqx, pg) = (event[1], event[2], event[3])
        cosThetaqg = pq.space().dot(pg.space())/(abs(pq.space())*abs(pg.space()))
        cosThetaqxg = pqx.space().dot(pg.space())/(abs(pqx.space())*abs(pg.space()))
        if (COSTHETA_CUT is not None and abs(cosThetaqxg) > COSTHETA_CUT) and \
           (COSTHETA_CUT is not None and abs(cosThetaqg) > COSTHETA_CUT) and \
           (GLUON_ENERGY_CUT is not None and pg[0] > GLUON_ENERGY_CUT):
            return False
    elif len(event) == 3:
        (pq, pqx) = (event[1], event[2])
        return True
    return True


def get_kinematic_variables(Q, p_d, p_dx, p_g) -> dict[str, float]:

    qsq = Q.dot(Q)
    x_1 = 2*p_d.dot(Q) / qsq
    x_2 = 2*p_dx.dot(Q) / qsq
    x_3 = 2*p_g.dot(Q) / qsq

    var = {
        'qsq': qsq,
        'x_1': x_1,
        'x_2': x_2,
        'x_3': x_3
    }
    for i, p_i in enumerate([p_d, p_dx, p_g]):
        for j, p_j in enumerate([p_d, p_dx, p_g]):
            if i == j:
                continue
            var[f'y_{i+1}{j+1}'] = 2*p_i.dot(p_j)/qsq

    return var


def integrand(event_file: CHEPEvent, costheta_histogram, ps_generator, model, E_cm: float, born_process, real_emission_process, samples_batch: list[Sample]) -> list[float]:

    global max_real

    evaluations: list[float] = []
    for sample in samples_batch:
        ps_point, jacobian = ps_generator.generateKinematics(E_cm, sample.c)

        ps_point = LorentzVectorList([ps_point[0]+ps_point[1]]+ps_point[2:])

        p_gstar, p_d, p_dx, p_g = ps_point[0], ps_point[1], ps_point[2], ps_point[3]
        cosThetaqq = p_d.space().dot(p_dx.space())/(abs(p_d.space())*abs(p_dx.space()))

        Q = p_gstar
        var = get_kinematic_variables(Q, p_d, p_dx, p_g)

        # Correct for the fact that we faked a 1 > 2 generation from a 1 > 3 generation
        jacobian *= 1. / (16.*math.pi**5)

        # Build reduced kinematics
        D13_2_kin = [p_gstar,]+[
            Q - p_dx*(1/var["x_2"]),
            p_dx*(1/var["x_2"])
        ]

        D23_1_kin = [p_gstar,]+[
            Q - p_d*(1/var["x_1"]),
            p_d*(1/var["x_1"])
        ]

        # jacobian is equal to : (Q^2 / 16 pi^2) dx_1 dx_2 Theta (x1+x2-1)
        # but it is already included in our phase-space parameterization jacobian

        real_me_normalization_factor = (
            4.0/3.0)*(8.0 * math.pi * model.aS)/var["qsq"]

        subtraction_phase_space_measure_factor = 1/2.

        # Subtract both ends of the dipole
        S13_2 = - born_process.smatrix(D13_2_kin, model) * real_me_normalization_factor * (
            1/(1-var["x_2"])*(2 / (2 - var["x_1"] - var["x_2"]) -
                              (1 + var["x_1"])) + (1-var["x_1"])/var["x_2"]
        )

        S23_1 = - born_process.smatrix(D23_1_kin, model) * real_me_normalization_factor * (
            1/(1-var["x_1"])*(2 / (2 - var["x_1"] - var["x_2"]) -
                              (1 + var["x_2"])) + (1-var["x_2"])/var["x_1"]
        )

        # For a 2 > N process, the flux factor is as given below
        # initial_state_flux = 1.0/(8*math.pi**2*E_cm**2)
        # For a 1 > N decay process, the effective flux factor is just
        initial_state_flux = 1.0/(2*E_cm)

        does_pass_real_emission_cut = pass_cuts(ps_point)
        does_pass_s13_2_cut = pass_cuts(D13_2_kin)
        does_pass_s23_1_cut = pass_cuts(D23_1_kin)

        first_event = True

        if does_pass_real_emission_cut:
            real_emission_matrix_element_evaluation = real_emission_process.smatrix(
                ps_point, model)

            # You can compare the above evaluation to the analytical expression below
            recomputed_real_emission_wgt = real_me_normalization_factor * born_process.smatrix(D23_1_kin, model) * \
                ((var["x_1"]**2 + var["x_2"]**2) /
                 ((1.0-var["x_1"])*(1.0-var["x_2"])))
            recomputed_real_emission_wgt *= subtraction_phase_space_measure_factor * jacobian * \
                initial_state_flux * UNIT_CONVERSION_FACTOR

            real_emission_wgt = subtraction_phase_space_measure_factor * real_emission_matrix_element_evaluation * \
                jacobian * initial_state_flux * UNIT_CONVERSION_FACTOR

            event = CHEPEvent()
            event.nexternal = 4
            event.wgt = real_emission_wgt
            event.aqed = 1
            event.aqcd = 1
            event.extend([
                CHEPParticle(event, LegState.INITIAL, 23, p_gstar[1], p_gstar[2], p_gstar[3], p_gstar[0], mass=model.mdl_MZ),  # nopep8
                CHEPParticle(event, LegState.FINAL, 1, p_d[1], p_d[2], p_d[3], p_d[0], mass=0.0),  # nopep8
                CHEPParticle(event, LegState.FINAL, -1, p_dx[1], p_dx[2], p_dx[3], p_dx[0], mass=0.0),  # nopep8
                CHEPParticle(event, LegState.FINAL, 21, p_g[1], p_g[2], p_g[3], p_g[0], mass=0.0),  # nopep8
            ])
            if first_event:
                event_file.write("<eventgroup>\n")
                first_event = False
            event_file.write_events(event)

        else:
            real_emission_wgt = 0.
            recomputed_real_emission_wgt = 0.

        if does_pass_s13_2_cut:
            s13_2_wgt = S13_2 * subtraction_phase_space_measure_factor * \
                jacobian * initial_state_flux * UNIT_CONVERSION_FACTOR

            # Add the integrated counterterm and virtual weight
            born_wgt = born_process.smatrix(
                D23_1_kin, model) * (jacobian / (var["qsq"]/(16.0*math.pi**2))) * initial_state_flux * UNIT_CONVERSION_FACTOR

            I_V_wgt = born_wgt * ((4.0/3.0) * model.aS / math.pi)

            # Add the Born weight
            event = CHEPEvent()
            event.nexternal = 3
            event.wgt = born_wgt
            event.aqed = 1
            event.aqcd = 0
            event.extend([
                CHEPParticle(event, LegState.INITIAL, 23, p_gstar[1], p_gstar[2], p_gstar[3], p_gstar[0], mass=model.mdl_MZ),  # nopep8
                CHEPParticle(event, LegState.FINAL, 1, D13_2_kin[1][1], D13_2_kin[1][2], D13_2_kin[1][3], D13_2_kin[1][0], mass=0.0),  # nopep8
                CHEPParticle(event, LegState.FINAL, -1, D13_2_kin[2][1], D13_2_kin[2][2], D13_2_kin[2][3], D13_2_kin[2][0], mass=0.0),  # nopep8
            ])
            if first_event:
                event_file.write("<eventgroup>\n")
                first_event = False
            event_file.write_events(event)

            # Add the counter event s13_2
            event = CHEPEvent()
            event.nexternal = 3
            event.wgt = s13_2_wgt + I_V_wgt
            event.aqed = 1
            event.aqcd = 1
            event.extend([
                CHEPParticle(event, LegState.INITIAL, 23, p_gstar[1], p_gstar[2], p_gstar[3], p_gstar[0], mass=model.mdl_MZ),  # nopep8
                CHEPParticle(event, LegState.FINAL, 1, D13_2_kin[1][1], D13_2_kin[1][2], D13_2_kin[1][3], D13_2_kin[1][0], mass=0.0),  # nopep8
                CHEPParticle(event, LegState.FINAL, -1, D13_2_kin[2][1], D13_2_kin[2][2], D13_2_kin[2][3], D13_2_kin[2][0], mass=0.0),  # nopep8
            ])
            if first_event:
                event_file.write("<eventgroup>\n")
                first_event = False
            event_file.write_events(event)
        else:
            s13_2_wgt = 0.
            I_V_wgt = 0.

        if does_pass_s23_1_cut:
            s23_1_wgt = S23_1 * subtraction_phase_space_measure_factor * \
                jacobian * initial_state_flux * UNIT_CONVERSION_FACTOR

            # Add the counter event s23_1
            event = CHEPEvent()
            event.nexternal = 3
            event.wgt = s13_2_wgt
            event.aqed = 1
            event.aqcd = 1
            event.extend([
                CHEPParticle(event, LegState.INITIAL, 23, p_gstar[1], p_gstar[2], p_gstar[3], p_gstar[0], mass=model.mdl_MZ),  # nopep8
                CHEPParticle(event, LegState.FINAL, 1, D23_1_kin[1][1], D23_1_kin[1][2], D23_1_kin[1][3], D23_1_kin[1][0], mass=0.0),  # nopep8
                CHEPParticle(event, LegState.FINAL, -1, D23_1_kin[2][1], D23_1_kin[2][2], D23_1_kin[2][3], D23_1_kin[2][0], mass=0.0),  # nopep8
            ])
            if first_event:
                event_file.write("<eventgroup>\n")
                first_event = False
            event_file.write_events(event)
        else:
            s23_1_wgt = 0.

        if not first_event:
            event_file.write("</eventgroup>\n")

        # bin_id = int((1+cosThetaqq)/2*len(costheta_histogram))
        # costheta_histogram[bin_id][0] += wgt
        # costheta_histogram[bin_id][1] += 1

        # if max_real < abs(real_emission_wgt):
        #     from pprint import pprint
        #     max_real = abs(real_emission_wgt)
        #     pprint(var)
        #     print(ps_point)
        #     print("real_emission_wgt=", real_emission_wgt)
        #     print("recomputed_real_emission_wgt=",
        #           recomputed_real_emission_wgt)
        #     if recomputed_real_emission_wgt != 0.:
        #         print("real_emission_wgt / recomputed_real_emission_wgt = ",
        #               real_emission_wgt / recomputed_real_emission_wgt)
        #     print("s13_2_wgt=", s13_2_wgt)
        #     print("s23_1_wgt=", s23_1_wgt)
        #     print("s13_2_wgt+s23_1_wgt=", s13_2_wgt+s23_1_wgt)
        #     print("(s13_2_wgt+s23_1_wgt) / real_emission_wgt =",
        #           (s13_2_wgt+s23_1_wgt)/real_emission_wgt)

        final_wgt = born_wgt + real_emission_wgt + s13_2_wgt + s23_1_wgt + I_V_wgt
        # final_wgt = born_wgt
        evaluations.append(final_wgt)
        # evaluations.append(I_V_wgt)
        # evaluations.append(born_wgt)

    return evaluations


def rratio_subtracted(args: argparse.Namespace):

    model = ModelParameters(None)
    # print(model.aS)
    # print(model.aEWM1)
    # stop
    real_emission_process = Matrix_2_z_ddxg()
    born_process = Matrix_1_z_ddx()
    external_masses = real_emission_process.get_external_masses(model)

    E_cm = model.mdl_MZ

    ps_generator = FlatPhaseSpace(
        [external_masses[0]/2., external_masses[0]/2.], external_masses[1],
        beam_Es=(E_cm/2., E_cm/2.),
        # We do not consider PDF for e+ e- > l+ l- at fixed-order
        beam_types=(0, 0)
    )

    n_dimensions = ps_generator.nDimPhaseSpace()

    N_CORES = 1  # Parallelization not implemented yet
    DISCRETE_LEARNING_RATE = 0.15
    CONTINUOUS_LEARNING_RATE = 0.15

    event_file = CHEPEventFile("./epem_ddxg_subtracted.lhe", mode='w')

    costheta_histo = [[0.0, 0] for _ in range(200)]
    parallel_rngs = [RandomNumberGenerator(
        seed=args.seed, stream_id=i_core) for i_core in range(N_CORES)]
    integrator = NumericalIntegrator.continuous(n_dimensions)
    for i_iteration in range(args.n_iterations):
        samples = integrator.sample(
            args.n_points_per_iteration, parallel_rngs[0])
        res = integrand(event_file, costheta_histo, ps_generator,
                        model, E_cm, born_process, real_emission_process, samples)
        integrator.add_training_samples(samples, res)
        avg, err, chi_sq = integrator.update(
        discrete_learning_rate=DISCRETE_LEARNING_RATE,
        continuous_learning_rate=CONTINUOUS_LEARNING_RATE)  # type: ignore # nopep8
        logger.info(
            'Iteration {}: {:.6} +- {:.6}, chi={:.6}'.format(i_iteration, avg, err, chi_sq))

    event_file.write("</LesHouchesEvents>\n")

    banner = event_file.get_banner()
    banner.modify_init_cross({1: avg})
    event_file.seek(0)
    banner.write(event_file, close_tag=False)
    event_file.close()
    # unweighted_event_file = CHEPEventFile("./epem_ddxg.lhe", mode='r')
    # unweighted_event_file.unweight("./unweighted_epem_ddxg.lhe")
    # unweighted_event_file.close()

    normalize_histogram = [((tot_wgt/n_events, n_events) if n_events > 0 else (0., 0)) for (
        tot_wgt, n_events) in costheta_histo]

    # use matplotlib to plot the histogram
    import matplotlib.pyplot as plt
    import numpy as np
    bins = np.linspace(-1, 1, len(normalize_histogram))
    plt.bar(bins, [wgt for (wgt, _)
                   in normalize_histogram], width=0.01)
    plt.xlabel("Cosine Theta")
    plt.ylabel("Weight")
    plt.title("Cosine Theta Distribution")
    plt.savefig("epem_ddxg_costhetaqq.png")
