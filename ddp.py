import argparse
import logging

import operator
from functools import reduce

import os

from utils.common import epilog, enable_logging, version, get_immediate_subdirectories
from utils.genome import parse_genome_in_grimm_file
from impl_gurobi.double_distance import create_ilp_formulation_for_ddp
from utils.set_definer import complete_genes_multiset, \
    define_equiv_function, vertex_set_from_gene_multiset

logger = logging.getLogger()


class DoubleDistConf(object):
    def __init__(self, ordinary_genome, duplicated_genome, name, log_file, tl, mult=2):
        if mult < 2 or mult > 3:
            raise Exception("Unsupported multiplication of a genome")

        self.multiplicity = mult
        self.genes_of_dupl_genome = duplicated_genome.get_gene_multiset()
        self.genes_of_ord_genome = ordinary_genome.get_gene_multiset()

        self.s_all_genes = complete_genes_multiset(
            self.genes_of_ord_genome.keys() | self.genes_of_dupl_genome.keys(), 1)
        self.ms_all_genes = complete_genes_multiset(
            self.genes_of_ord_genome.keys() | self.genes_of_dupl_genome.keys(), mult)

        # Coding breakpoint graph
        self.bg_vertex_set = vertex_set_from_gene_multiset(self.ms_all_genes)
        self.bg_ind2vertex = [''] + [u for u in self.bg_vertex_set]
        self.bg_vertex2ind = {self.bg_ind2vertex[i]: i for i in range(1, len(self.bg_ind2vertex))}

        bg_A_matching, bg_A_telomers = duplicated_genome.convert_to_genome_graph()
        self.ind_bg_A_vertices = {self.bg_vertex2ind[u] for u in
                                  vertex_set_from_gene_multiset(self.genes_of_dupl_genome)}
        self.ind_bg_A_edges = {tuple(sorted((self.bg_vertex2ind[u], self.bg_vertex2ind[v]))) for u, v in bg_A_matching}
        self.ind_bg_A_telomers = {self.bg_vertex2ind[u] for u in bg_A_telomers}

        # Coding contracted breakpoint graph
        self.cbg_vertex_set = vertex_set_from_gene_multiset(self.s_all_genes)
        self.cbg_ind2vertex = [''] + [u for u in self.cbg_vertex_set]
        self.cbg_vertex2ind = {self.cbg_ind2vertex[i]: i for i in range(1, len(self.cbg_ind2vertex))}

        cbg_R_matching, cbg_R_telomers = ordinary_genome.convert_to_genome_graph()
        self.ind_cbg_R_vertices = {self.cbg_vertex2ind[u] for u in
                                   vertex_set_from_gene_multiset(self.genes_of_ord_genome)}
        self.ind_cbg_R_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in
                                cbg_R_matching}
        self.ind_cbg_R_telomers = {self.cbg_vertex2ind[u] for u in cbg_R_telomers}

        self.equiv_map = define_equiv_function(self.ms_all_genes, self.cbg_vertex2ind, self.bg_vertex2ind)

        # Coding completion of genomes
        self.ind_compl_for_A = set(self.bg_vertex2ind.values()) - self.ind_bg_A_vertices
        self.ind_compl_for_R = set(self.bg_vertex2ind.values()) - reduce(operator.or_, [set(self.equiv_map[v])
                                                                                        for v in
                                                                                        self.ind_cbg_R_vertices], set())

        self.biggest_const = len(self.bg_ind2vertex)
        self.name_model = name
        self.log_file = log_file
        self.time_limit = tl

    def get_biggest_const(self):
        return self.biggest_const


def solve_double_distance_problem(ord_genome_file, all_dupl_genome_file, out_result_file, gurobi_log_file, mult=2,
                                  time_limit=7200):
    logging.info('Start to solve DDP problem with time limit equals {0}'.format(time_limit))

    ord_genome = parse_genome_in_grimm_file(ord_genome_file)
    all_dupl_genome = parse_genome_in_grimm_file(all_dupl_genome_file)

    # ord_genome = remove_singletons_in_ord_wrt_two_dupl(ord_genome, all_dupl_genome)
    # all_dupl_genome = remove_singletons_dupl_wrt_gene_set(all_dupl_genome, set(ord_genome.get_gene_multiset().keys()))

    logging.info('Create ILP config')
    config = DoubleDistConf(ordinary_genome=ord_genome, duplicated_genome=all_dupl_genome, name="DDP",
                            log_file=gurobi_log_file, tl=time_limit, mult=mult)

    answer = create_ilp_formulation_for_ddp(config=config)
    if answer is not None:
        logging.info('Save results.')
        answer.write_stats_file(out_result_file)
    else:
        logging.info('There are no answers. Please, check log file.')


def do_test_job(input_directory, gurobi_log_file):
    logging.info('Let\'s do comparison test for the DDP-ILP. Input directory is {0}'.format(input_directory))

    for path, name in get_immediate_subdirectories(input_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            logging.info("Working with directory {0}".format(subpath))
            ilppath = os.path.join(subpath, "ilp")

            ordinary_genome = os.path.join(ilppath, "target_genome.gen")
            all_dupl_genome = os.path.join(ilppath, "a.gen")
            out_file = os.path.join(ilppath, "result.txt")
            solve_double_distance_problem(ord_genome_file=ordinary_genome, all_dupl_genome_file=all_dupl_genome,
                                          out_result_file=out_file, gurobi_log_file=gurobi_log_file)


def test_run():
    parser = argparse.ArgumentParser(description="ILP solver for the Double Distance Problem (DDP)",
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
        parser.error("The directory must exist and have the correct structure")
    args.paper_test = os.path.abspath(args.paper_test)

    args.log_file = os.path.join(args.paper_test, "double_distance.log")
    args.gurobi_log_file = os.path.join(args.paper_test, "gurobi_double_distance.log")
    enable_logging(args.log_file, overwrite=False)

    do_test_job(args.paper_test, args.gurobi_log_file)


def main():
    parser = argparse.ArgumentParser(description="ILP solver for the Double Distance Problem (DDP)",
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=epilog())

    parser.add_argument("-r", "--ordinary", dest="ordinary",
                        default=None, metavar="PATH", help="Ordinary genome in GRIMM format")

    parser.add_argument("-t", "--two_dupl", dest="duplicated",
                        default=None, metavar="PATH", help="2-duplicated genome in GRIMM format")

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

    args.log_file = os.path.join(args.out_dir, "double_distance.log")
    args.gurobi_log_file = os.path.join(args.out_dir, "gurobi_double_distance.log")
    enable_logging(args.log_file, overwrite=False)

    path_to_ord_genome = os.path.abspath(args.ordinary)
    path_to_dupl_genome = os.path.abspath(args.duplicated)

    result_out_file = os.path.join(args.out_dir, "result.txt")

    solve_double_distance_problem(ord_genome_file=path_to_ord_genome,
                                  all_dupl_genome_file=path_to_dupl_genome,
                                  out_result_file=result_out_file,
                                  gurobi_log_file=args.gurobi_log_file,
                                  time_limit=args.time_limit)


if __name__ == "__main__":
    main()
