import argparse
import logging

import os

from impl_gurobi.common import get_immediate_subdirectories
from utils.common import epilog, enable_logging, version
from utils.medians import medians_without_singletons


def do_test_job(input_directory, gurobi_log_file):
    logging.info('Let\'s do comparison test for the GMP-ILP. Input directory is {0}'.format(input_directory))

    for path, name in get_immediate_subdirectories(input_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            logging.info("Working with directory {0}".format(subpath))
            ilppath = os.path.join(subpath, "force_ilp")

            genome1 = os.path.join(ilppath, "S1.gen")
            genome2 = os.path.join(ilppath, "S2.gen")
            genome3 = os.path.join(ilppath, "S4.gen")

            result_out_file = os.path.join(ilppath, "result.txt")
            genome_out_file = os.path.join(ilppath, "median.gen")

            medians_without_singletons(genome_files=[genome1, genome2, genome3],
                                       out_result_file=result_out_file,
                                       out_median_file=genome_out_file,
                                       problem="GMP",
                                       gurobi_log_file=gurobi_log_file,
                                       time_limit=7200
                                       )


def test_run():
    parser = argparse.ArgumentParser(description="ILP solver for the Genome Median Problem (GMP)",
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=epilog())

    parser.add_argument("-pt", "--paper_test", dest="paper_test",
                        default=None, metavar="PATH",
                        help="The directory with data from the paper. Please use only "
                             "for reproducibility of the paper results")

    parser.add_argument("-tl", "--time_limit", dest="time_limit",
                        default=7200, required=False,
                        type=int, metavar="NUMBER",
                        help="Time limit for Gurobi solver")

    parser.add_argument("-v", "--version", action="version", version=version())

    args = parser.parse_args()

    if not os.path.isdir(args.paper_test):
        parser.error("The directory must exist and have correct structure")
    args.paper_test = os.path.abspath(args.paper_test)

    args.log_file = os.path.join(args.paper_test, "gmp.log")
    args.gurobi_log_file = os.path.join(args.paper_test, "gurobi_gmp.log")
    enable_logging(args.log_file, overwrite=False)

    do_test_job(args.paper_test, args.gurobi_log_file)


def main():
    parser = argparse.ArgumentParser(description="ILP solver for the Genome Median Problem (GMP)",
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=epilog())

    parser.add_argument("-f", "--files", dest="files",
                        default=None, metavar="PATH", nargs=3,
                        help="Three files in GRIMM format")

    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="PATH", help="Output directory")

    parser.add_argument("-tl", "--time_limit", dest="time_limit",
                        default=7200, required=False,
                        type=int, metavar="NUMBER",
                        help="Time limit for Gurobi solver")

    parser.add_argument("-v", "--version", action="version", version=version())

    args = parser.parse_args()

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)
    args.out_dir = os.path.abspath(args.out_dir)

    args.log_file = os.path.join(args.out_dir, "gmp.log")
    args.gurobi_log_file = os.path.join(args.out_dir, "gurobi_gmp.log")
    enable_logging(args.log_file, overwrite=False)

    path_to_genomes = [os.path.abspath(f) for f in args.files]
    result_out_file = os.path.join(args.out_dir, "result.txt")
    genome_out_file = os.path.join(args.out_dir, "median.gen")

    medians_without_singletons(genome_files=path_to_genomes,
                               out_result_file=result_out_file,
                               out_median_file=genome_out_file,
                               problem="GMP",
                               gurobi_log_file=args.gurobi_log_file,
                               time_limit=args.time_limit)


if __name__ == "__main__":
    main()