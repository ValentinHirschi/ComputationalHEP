import argparse
import math

from symbolica import NumericalIntegrator, RandomNumberGenerator, Sample

from CHEP.matrix_elements.madgraph.model.parameters import ModelParameters
from CHEP.matrix_elements.madgraph.processes.all_processes import (
    Matrix_1_epem_ddx_no_z,
    Matrix_3_epem_ddxg_no_z,
)
from CHEP.phase_space_generators.phase_space_generators import FlatPhaseSpace
from CHEP.shower.durham import leading_jet_angle_23
from CHEP.shower.event import make_epem_event, shower_event_from_chep_event
from CHEP.shower.histograms import AngleHistogram
from CHEP.shower.qcd import AlphaS
from CHEP.shower.shower import Shower
from CHEP.utils import CHEPException, logger
from CHEP.utils.lhe_parser import CHEPEventFile

GEV_TO_PB = 0.389379338e9

def process_setup(process_name: str):
    match process_name:
        case "epem_a_ddx":
            return {
                "matrix": Matrix_1_epem_ddx_no_z(),
                "final_pids": [1, -1],
                "final_colors": [[501, 0], [0, 501]],
                "aqcd": 0,
                "fixed_order_file": "epem_a_ddx_fixed_order.lhe",
                "showered_file": "epem_a_ddx_showered.lhe",
                "histogram_prefix": "epem_a_ddx_theta23",
                "title": r"$e^+e^- \to \gamma^* \to d\bar{d}$ + shower",
            }
        case "epem_a_ddxg":
            return {
                "matrix": Matrix_3_epem_ddxg_no_z(),
                "final_pids": [1, -1, 21],
                "final_colors": [[502, 0], [0, 501], [501, 502]],
                "aqcd": 1,
                "fixed_order_file": "epem_a_ddxg_fixed_order.lhe",
                "showered_file": "epem_a_ddxg_showered.lhe",
                "histogram_prefix": "epem_a_ddxg_theta23",
                "title": r"$e^+e^- \to \gamma^* \to d\bar{d}g$ + shower",
            }
        case _:
            raise CHEPException(f"Shower process {process_name} not implemented")


def event_weight(ps_point, jacobian: float, model, process, e_cm: float) -> float:
    matrix_element = process.smatrix(ps_point, model)
    initial_state_flux = 1.0 / (8.0 * math.pi**2 * e_cm**2)
    return matrix_element * jacobian * initial_state_flux * GEV_TO_PB


def finalize_lhe(event_file: CHEPEventFile, cross_section: float) -> None:
    event_file.write("</LesHouchesEvents>\n")
    banner = event_file.get_banner()
    banner.modify_init_cross({1: cross_section})
    event_file.seek(0)
    banner.write(event_file, close_tag=False)
    event_file.close()


def integrand(
    fixed_order_file: CHEPEventFile,
    showered_file: CHEPEventFile,
    histogram: AngleHistogram,
    ps_generator,
    model,
    process,
    process_info,
    shower_generator: Shower,
    e_cm: float,
    step_shower: bool,
    samples_batch: list[Sample],
) -> list[float]:
    evaluations: list[float] = []
    ecm2 = e_cm**2

    for sample in samples_batch:
        ps_point, jacobian = ps_generator.generateKinematics(e_cm, sample.c)
        weight = event_weight(ps_point, jacobian, model, process, e_cm)

        fixed_order_event = make_epem_event(
            ps_point,
            process_info["final_pids"],
            process_info["final_colors"],
            weight,
            aqed=1,
            aqcd=process_info["aqcd"],
            scale=e_cm,
        )
        fixed_order_file.write_events(fixed_order_event)

        fixed_angle = leading_jet_angle_23(list(ps_point[2:]), ecm2)
        histogram.fill_fixed_order(fixed_angle, weight)

        shower_event = shower_event_from_chep_event(fixed_order_event)
        shower_generator.run(shower_event, hard_scale2=ecm2, step=step_shower)
        showered_file.write_events(shower_event.to_chep_event())

        shower_angle = leading_jet_angle_23(list(shower_event.final_momenta()), ecm2)
        histogram.fill_showered(shower_angle, weight)

        evaluations.append(weight)
    return evaluations


def shower(args: argparse.Namespace):
    process_info = process_setup(args.shower_process)

    if args.step_shower and args.n_cores != 1:
        logger.info("--step_shower forces n_cores = 1")
        args.n_cores = 1

    model = ModelParameters(None)
    process = process_info["matrix"]
    external_masses = process.get_external_masses(model)

    e_cm = args.e_cm
    ps_generator = FlatPhaseSpace(
        external_masses[0],
        external_masses[1],
        beam_Es=(e_cm / 2.0, e_cm / 2.0),
        beam_types=(0, 0),
    )
    n_dimensions = ps_generator.nDimPhaseSpace()

    fixed_order_path = args.event_file or process_info["fixed_order_file"]
    showered_path = args.showered_event_file or process_info["showered_file"]
    histogram_prefix = args.histogram_prefix or process_info["histogram_prefix"]

    fixed_order_file = CHEPEventFile(fixed_order_path, mode="w")
    showered_file = CHEPEventFile(showered_path, mode="w")
    histogram = AngleHistogram(args.n_bins)

    alpha_s = AlphaS(91.1876, model.aS)
    shower_generator = Shower(alpha_s, t0=args.shower_scale**2, seed=args.seed)

    discrete_learning_rate = 0.15
    continuous_learning_rate = 0.15
    parallel_rngs = [
        RandomNumberGenerator(seed=args.seed, stream_id=i_core)
        for i_core in range(args.n_cores)
    ]
    integrator = NumericalIntegrator.continuous(n_dimensions)

    avg = 0.0
    for i_iteration in range(args.n_iterations):
        samples = integrator.sample(args.n_points_per_iteration, parallel_rngs[0])
        res = integrand(
            fixed_order_file,
            showered_file,
            histogram,
            ps_generator,
            model,
            process,
            process_info,
            shower_generator,
            e_cm,
            args.step_shower,
            samples,
        )
        integrator.add_training_samples(samples, res)
        avg, err, chi_sq = integrator.update(
            discrete_learning_rate=discrete_learning_rate,
            continuous_learning_rate=continuous_learning_rate,
        )
        logger.info(
            "Iteration {}: {:.6} +- {:.6}, chi={:.6}".format(
                i_iteration, avg, err, chi_sq
            )
        )

    finalize_lhe(fixed_order_file, avg)
    finalize_lhe(showered_file, avg)
    histogram.save(histogram_prefix, process_info["title"])
