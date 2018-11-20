import logging

import operator
from functools import reduce

import os

from genome import parse_genome_in_grimm_file
from impl_gurobi.double_distanse import create_ilp_formulation_for_ddp
from impl_gurobi.utils import remove_singletons_in_ord_wrt_two_dupl, remove_singletons_dupl_wrt_gene_set, \
    create_complete_genes_multiset, define_equiv_function, create_vertex_set_from_gene_multiset, \
    get_immediate_subdirectories

logger = logging.getLogger()


class DoubleDistConf(object):
    def __init__(self, ordinary_genome, duplicated_genome, name, log_file, tl, mult=2):
        if mult < 2 or mult > 3:
            raise Exception("Unsupported multiplication of a genome")

        self.multiplicity = mult
        self.genes_of_dupl_genome = duplicated_genome.get_gene_multiset()
        self.genes_of_ord_genome = ordinary_genome.get_gene_multiset()

        self.s_all_genes = create_complete_genes_multiset(self.genes_of_ord_genome.keys() | self.genes_of_dupl_genome.keys(), 1)
        self.ms_all_genes = create_complete_genes_multiset(self.genes_of_ord_genome.keys() | self.genes_of_dupl_genome.keys(), mult)

        # Coding breakpoint graph
        self.bg_vertex_set = create_vertex_set_from_gene_multiset(self.ms_all_genes)
        self.bg_ind2vertex = [''] + [u for u in self.bg_vertex_set]
        self.bg_vertex2ind = {self.bg_ind2vertex[i]: i for i in range(1, len(self.bg_ind2vertex))}

        bg_A_matching, bg_A_telomers = duplicated_genome.convert_to_genome_graph()
        self.ind_bg_A_vertices = {self.bg_vertex2ind[u] for u in
                                  create_vertex_set_from_gene_multiset(self.genes_of_dupl_genome)}
        self.ind_bg_A_edges = {tuple(sorted((self.bg_vertex2ind[u], self.bg_vertex2ind[v]))) for u, v in bg_A_matching}
        self.ind_bg_A_telomers = {self.bg_vertex2ind[u] for u in bg_A_telomers}

        # Coding contracted breakpoint graph
        self.cbg_vertex_set = create_vertex_set_from_gene_multiset(self.s_all_genes)
        self.cbg_ind2vertex = [''] + [u for u in self.cbg_vertex_set]
        self.cbg_vertex2ind = {self.cbg_ind2vertex[i]: i for i in range(1, len(self.cbg_ind2vertex))}

        cbg_R_matching, cbg_R_telomers = ordinary_genome.convert_to_genome_graph()
        self.ind_cbg_R_vertices = {self.cbg_vertex2ind[u] for u in
                                   create_vertex_set_from_gene_multiset(self.genes_of_ord_genome)}
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


def ddp(ord_genome_file, all_dupl_genome_file, out_result_file, mult=2):
    ord_genome = parse_genome_in_grimm_file(ord_genome_file)
    all_dupl_genome = parse_genome_in_grimm_file(all_dupl_genome_file)

    ord_genome = remove_singletons_in_ord_wrt_two_dupl(ord_genome, all_dupl_genome)
    all_dupl_genome = remove_singletons_dupl_wrt_gene_set(all_dupl_genome, set(ord_genome.get_gene_multiset().keys()))

    config = DoubleDistConf(ordinary_genome=ord_genome, duplicated_genome=all_dupl_genome, name="DDP",
                            log_file="gurobi_double_dist.log", tl=7200, mult=mult)
    answer = create_ilp_formulation_for_ddp(config=config)
    answer.write_stats_file(out_result_file)


def dojob(name_directory):
    logging.info('Let us do the ILP')

    for path, name in get_immediate_subdirectories(name_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            logging.info("Working with directory {0}".format(subpath))
            ilppath = os.path.join(subpath, "ilp")

            ordinary_genome = os.path.join(ilppath, "target_genome.gen")
            all_dupl_genome = os.path.join(ilppath, "a.gen")
            out_file = os.path.join(ilppath, "result.txt")
            ddp(ord_genome_file=ordinary_genome, all_dupl_genome_file=all_dupl_genome, out_result_file=out_file)


if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     sys.stderr.write("USAGE: double_distance.py <destination>\n")
    #     sys.exit(1)
    #
    # if os.path.isfile(sys.argv[1]):
    #     sys.stderr.write("<destination> - is directory with specific hierarchy\n")
    #     sys.exit(1)
    #
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(logging.INFO)

    # log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
    #                                   "%(message)s", "%H:%M:%S")
    # file_log = logging.FileHandler("txt_double_distance.log", mode="w")
    # file_log.setFormatter(log_formatter)

    logger.setLevel(logging.INFO)
    logger.addHandler(console_log)
    # logger.addHandler(file_log)
    #
    # dir_path = os.path.abspath(sys.argv[1])
    # dojob(dir_path)

    ddp("r1.gen", "a1.gen", "result.txt")
