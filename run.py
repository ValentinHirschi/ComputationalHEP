#!/usr/bin/env python3
import argparse
import os
import sys
import tempfile

# Some tutorial environments have a non-writable default Matplotlib cache
# directory.  Set a writable cache before importing experiment modules, because
# a few of them import pyplot at module-import time.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "chep-matplotlib"))
os.environ.setdefault("MPLBACKEND", "Agg")

from CHEP.experiments.epem_lplm_fixed_order_LO import epem_lplm_fixed_order_LO
from CHEP.experiments.sampling_experiment import sampling_experiment
from CHEP.experiments.rratio import rratio
from CHEP.experiments.rratio_differential import rratio_differential
from CHEP.experiments.rratio_subtracted import rratio_subtracted
from CHEP.experiments.pp_wpwm_fixed_order_LO import pp_wpwm_fixed_order_LO, pdf_constraints_test
from CHEP.experiments.proton_pdf import DEFAULT_PDF_SET, proton_pdf
from CHEP.experiments.drell_yan_lo import drell_yan_lo
from CHEP.experiments.shower import shower

from CHEP.utils import logger, CHEPException, setup_logging

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

parser_proton_pdf = subparsers.add_parser(
    'proton_pdf',
    help='Study proton PDFs, sum rules, and a small pedagogical DGLAP evolution.')
parser_proton_pdf.add_argument('--lhapdf_python_dir', '--lhpadf_python_dir',
                               dest='lhapdf_python_dir', type=str,
                               default=default_lhapdf_python_dir,
                               help='Installation directory for the python3 LHAPDF module')
parser_proton_pdf.add_argument('--lhapdf_pdfsets_dir', '--lhpadf_pdfsets_dir',
                               dest='lhapdf_pdfsets_dir', type=str,
                               default=default_lhapdf_pdfsets_dir,
                               help='Directory containing PDF set data')
parser_proton_pdf.add_argument('--pdf_set', type=str, default=DEFAULT_PDF_SET,
                               help='Selected PDF set for the run')
parser_proton_pdf.add_argument('--pdf_member', type=int, default=0,
                               help='LHAPDF member to use. Member 0 is the central set.')
parser_proton_pdf.add_argument('--download_pdf', dest='download_pdf', action='store_true',
                               default=True,
                               help='Download the requested LHAPDF grid if it is missing')
parser_proton_pdf.add_argument('--no_download_pdf', dest='download_pdf', action='store_false',
                               help='Do not download missing LHAPDF grids')
parser_proton_pdf.add_argument('--output_dir', type=str, default='proton_pdf_output',
                               help='Directory where plots and tables are written')
parser_proton_pdf.add_argument('--x_min', type=float, default=1.0e-5,
                               help='Smallest x value used in plots and sum-rule quadrature')
parser_proton_pdf.add_argument('--x_max', type=float, default=0.999,
                               help='Largest x value used in plots and sum-rule quadrature')
parser_proton_pdf.add_argument('--x_points', type=int, default=240,
                               help='Number of x points for plots and sum-rule quadrature')
parser_proton_pdf.add_argument('--q_values', type=str, default='2.0,10.0,91.1876,1000.0',
                               help='Comma-separated Q=mu_F values in GeV')
parser_proton_pdf.add_argument('--flavours', type=str,
                               default='g,u,d,s,c,ubar,dbar,sbar,cbar',
                               help='Comma-separated flavours to plot')
parser_proton_pdf.add_argument('--skip_dglap', action='store_true',
                               help='Skip the pedagogical DGLAP ODE comparison')
parser_proton_pdf.add_argument('--dglap_start_q', type=float, default=2.0,
                               help='Starting Q value in GeV for the DGLAP comparison')
parser_proton_pdf.add_argument('--dglap_target_q', type=float, default=10.0,
                               help='Target Q value in GeV for the DGLAP comparison')
parser_proton_pdf.add_argument('--dglap_x_min', type=float, default=1.0e-3,
                               help='Smallest x value for the DGLAP ODE grid')
parser_proton_pdf.add_argument('--dglap_x_max', type=float, default=0.8,
                               help='Largest x value for the DGLAP ODE grid')
parser_proton_pdf.add_argument('--dglap_x_points', type=int, default=45,
                               help='Number of x-grid points for the DGLAP ODE')
parser_proton_pdf.add_argument('--dglap_z_points', type=int, default=80,
                               help='Number of z-quadrature points per convolution')
parser_proton_pdf.add_argument('--dglap_steps', type=int, default=4,
                               help='Number of Runge-Kutta steps in ln(Q^2)')
parser_proton_pdf.add_argument('--dglap_flavour', type=str, default='u',
                               help='Flavour shown in the DGLAP/LHAPDF comparison')

parser_drell_yan_lo = subparsers.add_parser(
    'drell_yan_lo',
    help='Compute pp -> Z -> l+ l- at LO with a 2->1 Bjorken-x phase space.')
parser_drell_yan_lo.add_argument('--n_iterations', '-ni', type=int, default=10,
                                 help='Number of integration iterations to run')
parser_drell_yan_lo.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                 help='Number of points per iteration')
parser_drell_yan_lo.add_argument('--seed', '-s', type=int, default=0,
                                 help='Random number generator seed')
parser_drell_yan_lo.add_argument('--lhapdf_python_dir', '--lhpadf_python_dir',
                                 dest='lhapdf_python_dir', type=str,
                                 default=default_lhapdf_python_dir,
                                 help='Installation directory for the python3 LHAPDF module')
parser_drell_yan_lo.add_argument('--lhapdf_pdfsets_dir', '--lhpadf_pdfsets_dir',
                                 dest='lhapdf_pdfsets_dir', type=str,
                                 default=default_lhapdf_pdfsets_dir,
                                 help='Directory containing PDF set data')
parser_drell_yan_lo.add_argument('--pdf_set', type=str, default=DEFAULT_PDF_SET,
                                 help='Selected PDF set for the run')
parser_drell_yan_lo.add_argument('--pdf_member', type=int, default=0,
                                 help='LHAPDF member to use. Member 0 is the central set.')
parser_drell_yan_lo.add_argument('--download_pdf', dest='download_pdf', action='store_true',
                                 default=True,
                                 help='Download the requested LHAPDF grid if it is missing')
parser_drell_yan_lo.add_argument('--no_download_pdf', dest='download_pdf', action='store_false',
                                 help='Do not download missing LHAPDF grids')
parser_drell_yan_lo.add_argument('--e_cm', type=float, default=13000.0,
                                 help='Proton-proton collider energy in GeV')
parser_drell_yan_lo.add_argument('--m_ll', type=float, default=None,
                                 help='Dilepton invariant mass in GeV. Defaults to m_Z.')
parser_drell_yan_lo.add_argument('--mu_f', type=float, default=None,
                                 help='Factorisation scale in GeV. Defaults to m_ll.')
parser_drell_yan_lo.add_argument('--quark_flavours', type=str, default='d,u,s,c',
                                 help='Comma-separated incoming quark flavours to include')
parser_drell_yan_lo.add_argument('--lepton', choices=('e', 'mu', 'tau'), default='mu',
                                 help='Charged-lepton flavour used for the Z branching ratio')
parser_drell_yan_lo.add_argument('--branching_ratio', type=float, default=None,
                                 help='Override BR(Z -> l+ l-). Use 1 for inclusive Z production.')
parser_drell_yan_lo.add_argument('--output_dir', type=str, default='drell_yan_lo_output',
                                 help='Directory where the rapidity plot and table are written')
parser_drell_yan_lo.add_argument('--n_bins', type=int, default=80,
                                 help='Number of rapidity bins for the output distribution')

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

parser_shower_experiment = subparsers.add_parser(
    'shower', help='Generate fixed-order e+e- events and run the pedagogical final-state shower.')
parser_shower_experiment.add_argument(
    'shower_process', choices=('epem_a_ddx', 'epem_a_ddxg'),
    help='Hard process to generate before showering.')
parser_shower_experiment.add_argument('--n_iterations', '-ni', type=int, default=10,
                                      help='Number of iterations to run')
parser_shower_experiment.add_argument('--n_points_per_iteration', '-npi', type=int, default=1000,
                                      help='Number of points per iteration to consider')
parser_shower_experiment.add_argument('--seed', '-s', type=int, default=0,
                                      help='Random number generator seed')
parser_shower_experiment.add_argument('--e_cm', type=float, default=1000.0,
                                      help='Center-of-mass energy in GeV')
parser_shower_experiment.add_argument('--shower_scale', type=float, default=1.0,
                                      help='Shower cutoff scale in GeV')
parser_shower_experiment.add_argument('--event_file', '-ef', type=str, default=None,
                                      help='Path to the fixed-order LHE file to write')
parser_shower_experiment.add_argument('--showered_event_file', '-sef', type=str, default=None,
                                      help='Path to the showered LHE file to write')
parser_shower_experiment.add_argument('--histogram_prefix', type=str, default=None,
                                      help='Prefix for the theta(j2,j3) histogram files')
parser_shower_experiment.add_argument('--n_bins', type=int, default=80,
                                      help='Number of histogram bins')
parser_shower_experiment.add_argument('--n_cores', '--n_core', dest='n_cores', type=int, default=1,
                                      help='Number of cores to use. Step shower mode forces this to 1.')
parser_shower_experiment.add_argument('--step_shower', action='store_true',
                                      help='Pause before each accepted shower emission and print the generated kinematics')


parser_rratio_analyze_events_experiment = subparsers.add_parser(
    'rratio_analyze_events', help='Analyze madgraph events to compute differential R-ratio quantities.')
parser_rratio_analyze_events_experiment.add_argument('--event_file', '-ef', type=str, default=None,
                                                     help='Specify the path to the madgraph event file to analyze')

if __name__ == "__main__":

    setup_logging()

    # The rest of this file uses argparse subcommands, whose names cannot begin
    # with "--".  The alias below honours the requested spelling
    # "./run.py --proton_pdf" while keeping the implementation consistent with
    # the other experiments.
    if len(sys.argv) > 1 and sys.argv[1] == '--proton_pdf':
        sys.argv[1] = 'proton_pdf'

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

        case 'shower':
            shower(args)

        case 'pp_wpwm_fixed_order_LO':
            pp_wpwm_fixed_order_LO(args)

        case 'pdf_constraints_test':
            pdf_constraints_test(args)

        case 'proton_pdf':
            proton_pdf(args)

        case 'drell_yan_lo':
            drell_yan_lo(args)

        case _:
            raise CHEPException(
                f"Experiment {args.experiment} not implemented")
