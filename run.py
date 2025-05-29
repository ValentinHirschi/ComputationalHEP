#!/usr/bin/env python3
import argparse
from CHEP.experiments.epem_lplm_fixed_order_LO import epem_lplm_fixed_order_LO
from CHEP.experiments.sampling_experiment import sampling_experiment
from CHEP.experiments.rratio import rratio
from CHEP.experiments.rratio_differential import rratio_differential
from CHEP.experiments.rratio_subtracted import rratio_subtracted
from CHEP.experiments.pp_wpwm_fixed_order_LO import pp_wpwm_fixed_order_LO, pdf_constraints_test

from CHEP.utils import logger, CHEPException, setup_logging
import os

root_path = os.path.dirname(__file__)
default_lhapdf_python_dir = os.path.normpath(os.path.abspath(os.path.join(
    root_path, os.path.pardir, 'HEPTools', 'lhapdf6_py3', 'lib', 'python3.12', 'site-packages')))
#default_lhapdf_python_dir = "/Users/vjhirsch/HEP_programs/HEPTools/lhapdf6_py3/lib/python3.12/site-packages"

default_lhapdf_pdfsets_dir = os.path.normpath(os.path.abspath(
    os.path.join(root_path, 'PDFsets')))

parser = argparse.ArgumentParser(prog='CHEP')
subparsers = parser.add_subparsers(
    title="experiment to run", dest="experiment", required=True,
    help='Various experiments available to run')

# Create the parser for the "epem_lplm_fixed_order_LO" experiment
parser_epem_lplm_fixed_order_LO = subparsers.add_parser(
    'epem_lplm_fixed_order_LO', help='Running e+ e- > l+ l- at fixed-order.')
parser_epem_lplm_fixed_order_LO.add_argument('--n_iterations', '-ni', type=int, default=10,
                                             help='Number of iterations to run')
parser_epem_lplm_fixed_order_LO.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                             help='Number of points per iteration to consider ')
parser_epem_lplm_fixed_order_LO.add_argument('--seed', '-s', type=int, default=0,
                                             help='Random number generator seed')

# Create the parser for the "pp_zz_fixed_order_LO" experiment
parser_pp_wpwm_fixed_order_LO = subparsers.add_parser(
    'pp_wpwm_fixed_order_LO', help='Running p p > z z at fixed-order.')
parser_pp_wpwm_fixed_order_LO.add_argument('--n_iterations', '-ni', type=int, default=10,
                                           help='Number of iterations to run')
parser_pp_wpwm_fixed_order_LO.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                           help='Number of points per iteration to consider ')
parser_pp_wpwm_fixed_order_LO.add_argument('--seed', '-s', type=int, default=0,
                                           help='Random number generator seed')
parser_pp_wpwm_fixed_order_LO.add_argument('--lhpadf_python_dir', type=str, default=default_lhapdf_python_dir,
                                           help='Installation directory for the python3 LHAPDF module')
parser_pp_wpwm_fixed_order_LO.add_argument('--lhpadf_pdfsets_dir', type=str, default=default_lhapdf_pdfsets_dir,
                                           help='Directory containing PDF sets data')
parser_pp_wpwm_fixed_order_LO.add_argument('--pdf_set', type=str, default='NNPDF23_nlo_as_0119',
                                           help='Selected PDF set for the run')

parser_pdf_constraints_test = subparsers.add_parser(
    'pdf_constraints_test', help='Testing baryon number and momentum conservation in PDFS.')
parser_pdf_constraints_test.add_argument('--n_iterations', '-ni', type=int, default=10,
                                         help='Number of iterations to run')
parser_pdf_constraints_test.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                         help='Number of points per iteration to consider ')
parser_pdf_constraints_test.add_argument('--seed', '-s', type=int, default=0,
                                         help='Random number generator seed')
parser_pdf_constraints_test.add_argument('--lhpadf_python_dir', type=str, default=default_lhapdf_python_dir,
                                         help='Installation directory for the python3 LHAPDF module')
parser_pdf_constraints_test.add_argument('--lhpadf_pdfsets_dir', type=str, default=default_lhapdf_pdfsets_dir,
                                         help='Directory containing PDF sets data')
parser_pdf_constraints_test.add_argument('--pdf_set', type=str, default='NNPDF23_nlo_as_0119',
                                         help='Selected PDF set for the run')

# Create the parser for the "sampling_experiment" experiment
parser_sampling_experiment = subparsers.add_parser(
    'sampling_experiment', help='Experiment momenta distribution from various samples.')
parser_sampling_experiment.add_argument('--seed', '-s', type=int, default=0,
                                        help='Random number generator seed')

parser_rratio_experiment = subparsers.add_parser(
    'rratio', help='Compute the R-ratio.')
parser_rratio_experiment.add_argument('--n_iterations', '-ni', type=int, default=10,
                                      help='Number of iterations to run')
parser_rratio_experiment.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                      help='Number of points per iteration to consider ')
parser_rratio_experiment.add_argument('--seed', '-s', type=int, default=0,
                                      help='Random number generator seed')

parser_rratio_differential_experiment = subparsers.add_parser(
    'rratio_differential', help='Generate events and distributions for the differential R-ratio.')
parser_rratio_differential_experiment.add_argument('--n_iterations', '-ni', type=int, default=10,
                                                   help='Number of iterations to run')
parser_rratio_differential_experiment.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                                   help='Number of points per iteration to consider ')
parser_rratio_differential_experiment.add_argument('--seed', '-s', type=int, default=0,
                                                   help='Random number generator seed')
parser_rratio_differential_experiment.add_argument('--event_file', '-ef', type=str, default='rratio_differential_events.lhe.gz',
                                                   help='Specify the path to the event file to write into.')

parser_rratio_subtracted_experiment = subparsers.add_parser(
    'rratio_subtracted', help='Generate events and distributions for the subtracted R-ratio.')
parser_rratio_subtracted_experiment.add_argument('--n_iterations', '-ni', type=int, default=10,
                                                 help='Number of iterations to run')
parser_rratio_subtracted_experiment.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                                 help='Number of points per iteration to consider ')
parser_rratio_subtracted_experiment.add_argument('--seed', '-s', type=int, default=0,
                                                 help='Random number generator seed')
parser_rratio_subtracted_experiment.add_argument('--event_file', '-ef', type=str, default='rratio_subtracted_events.lhe.gz',
                                                 help='Specify the path to the event file to write into.')


parser_rratio_analyze_events_experiment = subparsers.add_parser(
    'rratio_analyze_events', help='Analyze madgraph events to compute differential R-ratio quantities.')
parser_rratio_analyze_events_experiment.add_argument('--event_file', '-ef', type=str, default=None,
                                                     help='Specify the path to the madgraph event file to analyze')

if __name__ == "__main__":

    setup_logging()

    args = parser.parse_args()

    logger.info(f"Running experiment {args.experiment}")

    match args.experiment:

        case 'epem_lplm_fixed_order_LO':
            epem_lplm_fixed_order_LO(args)

        case 'sampling_experiment':
            sampling_experiment(args)

        case 'rratio':
            rratio(args)

        case 'rratio_differential':
            rratio_differential(args)

        case 'rratio_subtracted':
            rratio_subtracted(args)

        case 'pp_wpwm_fixed_order_LO':
            pp_wpwm_fixed_order_LO(args)

        case 'pdf_constraints_test':
            pdf_constraints_test(args)

        case _:
            raise CHEPException(
                f"Experiment {args.experiment} not implemented")
