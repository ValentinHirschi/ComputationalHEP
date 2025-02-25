#!/usr/bin/env python3
import argparse
from CHEP.experiments.epem_lplm_fixed_order_LO import epem_lplm_fixed_order_LO
from CHEP.experiments.sampling_experiment import sampling_experiment
from CHEP.utils import logger, CHEPException, setup_logging


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

# Create the parser for the "sampling_experiment" experiment
parser_sampling_experiment = subparsers.add_parser(
    'sampling_experiment', help='Experiment momenta distribution from various samples.')
parser_sampling_experiment.add_argument('--seed', '-s', type=int, default=0,
                                        help='Random number generator seed')


if __name__ == "__main__":

    setup_logging()

    args = parser.parse_args()

    logger.info(f"Running experiment {args.experiment}")

    match args.experiment:

        case 'epem_lplm_fixed_order_LO':
            epem_lplm_fixed_order_LO(args)

        case 'sampling_experiment':
            sampling_experiment(args)

        case _:
            raise CHEPException(
                f"Experiment {args.experiment} not implemented")
