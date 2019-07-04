import logging
import argparse

import itertools
from functools import reduce
import networkx as nx
import operator

import os

import sys

from utils.genome import parse_genome_in_grimm_file
from impl_gurobi.halvings import create_ilp_formulation_for_halvings_without_singletons
from impl_gurobi.common import create_complete_genes_multiset, remove_singletons_dupl_wrt_gene_set, \
    create_observed_edges_from_gene_multiset, create_vertex_set_from_gene_multiset, get_immediate_subdirectories, \
    define_equiv_function

logger = logging.getLogger()


class HalvingConf(object):
    def __init__(self, duplicated_genome, ordinary_genomes, name, log_file, tl, mult=2):
        if len(ordinary_genomes) != 1:
            raise Exception("Incorrect number of genomes for guided halving problems")

        if mult < 2 or mult > 3:
            raise Exception("Unsupported multiplication of a genome")

        self.gene_sets = [set(genome.get_gene_multiset().keys()) for genome in ordinary_genomes]
        self.genes_of_dupl_genome = duplicated_genome.get_gene_multiset()

        self.s_all_genes = create_complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()) |
                                                          self.genes_of_dupl_genome.keys(), 1)
        # Coding contracted breakpoint graph
        self.cbg_vertex_set = create_vertex_set_from_gene_multiset(self.s_all_genes)
        self.cbg_ind2vertex = [''] + [u for u in self.cbg_vertex_set]
        self.cbg_vertex2ind = {self.cbg_ind2vertex[i]: i for i in range(1, len(self.cbg_ind2vertex))}

        obverse_edges = create_observed_edges_from_gene_multiset(self.s_all_genes)
        self.ind_cbg_obverse_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v])))
                                      for u, v in obverse_edges}

        genome_graphs = [genome.convert_to_genome_graph() for genome in ordinary_genomes]
        self.ind_cbg_p_i_vertex_sets = [{self.cbg_vertex2ind[u] for u in
                                        create_vertex_set_from_gene_multiset(create_complete_genes_multiset(gene_set, 1))}
                                        for gene_set in self.gene_sets]
        self.ind_cbg_p_i_edges = [{tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in matching}
                                  for matching, _ in genome_graphs]
        self.ind_cbg_p_i_telomers = [{self.cbg_vertex2ind[u] for u in telomers} for _, telomers in genome_graphs]

        # This contracted genome graph does not contain parallel edges. Maybe a problem.
        cbg_A_matching, cbg_A_telomers = duplicated_genome.convert_to_contracted_genome_graph()
        self.ind_cbg_A_vertices = {self.cbg_vertex2ind[u] for u in
                                  create_vertex_set_from_gene_multiset(create_complete_genes_multiset(self.genes_of_dupl_genome.keys(), 1))}
        self.ind_cbg_A_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in cbg_A_matching}
        self.ind_cbg_A_telomers = {self.cbg_vertex2ind[u] for u in cbg_A_telomers}

        self.ms_all_genes = create_complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()) |
                                                           self.genes_of_dupl_genome.keys(), mult)
        # Coding breakpoint graph
        self.bg_vertex_set = create_vertex_set_from_gene_multiset(self.ms_all_genes)
        self.bg_ind2vertex = [''] + [u for u in self.bg_vertex_set]
        self.bg_vertex2ind = {self.bg_ind2vertex[i]: i for i in range(1, len(self.bg_ind2vertex))}

        bg_A_matching, bg_A_telomers = duplicated_genome.convert_to_genome_graph()
        self.ind_bg_A_vertices = {self.bg_vertex2ind[u] for u in create_vertex_set_from_gene_multiset(self.genes_of_dupl_genome)}
        self.ind_bg_A_edges = {tuple(sorted((self.bg_vertex2ind[u], self.bg_vertex2ind[v]))) for u, v in bg_A_matching}
        self.ind_bg_A_telomers = {self.bg_vertex2ind[u] for u in bg_A_telomers}

        # Connection between graphs
        self.equiv_map = define_equiv_function(self.ms_all_genes, self.cbg_vertex2ind, self.bg_vertex2ind)

        # TODO: HERE need to extend for indels
        self.ind_ancestral_set = self.ind_cbg_p_i_vertex_sets[0]
        self.ancestral_gene_set = self.gene_sets[0]

        self.multiplicity = mult
        self.number_of_genomes = len(ordinary_genomes)
        self.biggest_const = len(self.bg_ind2vertex)
        self.name_model = name
        self.log_file = log_file
        self.time_limit = tl

    def get_biggest_const(self):
        return self.biggest_const


class ClassicHalving(HalvingConf):
    def __init__(self, duplicated_genome, ordinary_genomes, name, log_file, tl, mult=2):
        super().__init__(duplicated_genome, ordinary_genomes, name, log_file, tl, mult)

        self.allowable_ancestral_edges = {tuple(sorted([u, v]))
                                          for u, v in itertools.combinations(self.ind_ancestral_set, 2)}

        flag = False

        if len(self.ind_bg_A_telomers) != 0:
            flag = True

        for i in range(self.number_of_genomes):
            if len(self.ind_cbg_p_i_telomers[i]) != 0:
                flag = True

        self.allowable_ancestral_telomers = {x: flag for x in self.ind_ancestral_set}
        self.connection_ancestral_constrs = {u: {tuple(sorted((u, v))) for v in self.ind_ancestral_set if u != v}
                                   for u in self.ind_ancestral_set}


class ConservedHalving(HalvingConf):
    def __init__(self, duplicated_genome, ordinary_genomes, name, log_file, tl, mult=2):
        super().__init__(duplicated_genome, ordinary_genomes, name, log_file, tl, mult)

        set_all_telomers = reduce(operator.or_, self.ind_cbg_p_i_telomers, set()) | self.ind_cbg_A_telomers
        self.allowable_ancestral_telomers = {x: True if x in set_all_telomers else False for x in self.ind_ancestral_set}

        graph = nx.MultiGraph()
        for edges in self.ind_cbg_p_i_edges:
            for u, v in edges:
                graph.add_edge(u, v)

        for u, v in self.ind_cbg_A_edges:
            graph.add_edge(u, v)

        self.allowable_ancestral_edges = set()
        self.connection_ancestral_constrs = dict()

        for v in self.ind_ancestral_set:
            connect_indices = set()

            if v in set(graph.nodes()):
                for u in graph.neighbors(v):
                    if u != v and u in self.ind_ancestral_set:
                        self.allowable_ancestral_edges.add(tuple(sorted((u, v))))
                        connect_indices.add(tuple(sorted((u, v))))

            self.connection_ancestral_constrs[v] = connect_indices


def halvings_without_singletons(ordinary_genome_file, all_dupl_genome_file, out_result_file, out_predup_file, problem):
    ord_genomes = [parse_genome_in_grimm_file(ordinary_genome_file)]
    all_dupl_genome = parse_genome_in_grimm_file(all_dupl_genome_file)
    # genomes = remove_known_singletons(genomes) # TODO Extend on indels

    if problem == "CGGHP":
        cfg = ConservedHalving(duplicated_genome=all_dupl_genome,
                               ordinary_genomes=ord_genomes,
                               name="CGGHP", log_file="gurobi_cgghp.log", tl=7200, mult=2)
    else:
        cfg = ClassicHalving(duplicated_genome=all_dupl_genome,
                             ordinary_genomes=ord_genomes,
                             name="GGHP", log_file="gurobi_gghp.log", tl=7200, mult=2)

    answer = create_ilp_formulation_for_halvings_without_singletons(cfg)

    answer.write_stats_file(out_result_file)
    answer.write_genome_file(out_predup_file)


def dojob_real(name_directory):
    logging.info('Let us do the ILP')

    for path, name in get_immediate_subdirectories(name_directory):
        if os.path.isfile(path):
            continue

        logging.info("Working with directory {0}".format(path))
        ilppath = os.path.join(path, "heur_ilp")

        ordinary_genome_file = os.path.join(ilppath, "B.gen")
        all_dupl_genome_file = os.path.join(ilppath, "A.gen")

        result_out_file = os.path.join(ilppath, "result.txt")
        genome_out_file = os.path.join(ilppath, "pre_dup.gen")

        halvings_without_singletons(ordinary_genome_file=ordinary_genome_file,
                                        all_dupl_genome_file=all_dupl_genome_file,
                                        out_result_file=result_out_file,
                                        out_predup_file=genome_out_file,
                                        problem="CGGHP")


def dojob(name_directory):
    logging.info('Let us do the ILP')

    for path, name in get_immediate_subdirectories(name_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            logging.info("Working with directory {0}".format(subpath))
            ilppath = os.path.join(subpath, "heur_ilp")

            ordinary_genome_file = os.path.join(ilppath, "B.gen")
            all_dupl_genome_file = os.path.join(ilppath, "A.gen")

            result_out_file = os.path.join(ilppath, "result.txt")
            genome_out_file = os.path.join(ilppath, "pre_dup.gen")

            halvings_without_singletons(ordinary_genome_file=ordinary_genome_file,
                                        all_dupl_genome_file=all_dupl_genome_file,
                                        out_result_file=result_out_file,
                                        out_predup_file=genome_out_file,
                                        problem="CGGHP")


if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="ILP for halving problem")
    # parser.add_argument("-g1", "--genome1", type=str, help="First genome in median problem.")
    # parser.add_argument("-g2", "--genome2", type=str, help="Second genome in median problem.")
    # parser.add_argument("-g3", "--genome3", type=str, help="Third genome in median problem.")
    #
    # console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
    #                                       "%(message)s", "%H:%M:%S")
    # console_log = logging.StreamHandler()
    # console_log.setFormatter(console_formatter)
    # console_log.setLevel(logging.INFO)
    #
    # logger.setLevel(logging.INFO)
    # logger.addHandler(console_log)
    #
    # halvings_without_singletons("r1.gen", "a1.gen", "result.txt", "ordinary.txt", "GGHP")

    # if len(sys.argv) < 2:
    #     sys.stderr.write("USAGE: double_distance.py <destination>\n")
    #     sys.exit(1)
    #
    # if os.path.isfile(sys.argv[1]):
    #     sys.stderr.write("<destination> - is directory with specific hierarchy\n")
    #     sys.exit(1)

    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(logging.INFO)

    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%H:%M:%S")

    file_log = logging.FileHandler("txt_heur_halving_real.log", mode="w")
    file_log.setFormatter(log_formatter)

    logger.setLevel(logging.INFO)
    logger.addHandler(console_log)
    logger.addHandler(file_log)

    dir_path = os.path.abspath(sys.argv[1])
    dojob_real(dir_path)
    # halvings_without_singletons("B.gen", "A.gen", "result.txt", "ordinary.txt", "CGGHP")


# def get_genome_set(ord_genome, all_dupl_genome):
#     s_1_a, s_2_a = all_dupl_genome.get_partitioned_gene_set()
#     return s_2_a | (s_1_a & ord_genome.get_gene_set())


# def remove_known_singletons_for_halving(ordinary_genome_file, all_dupl_genome_file):
#     ord_genome = parse_genome_in_grimm_file(ordinary_genome_file)
#     all_dupl_genome = parse_genome_in_grimm_file(all_dupl_genome_file)
#     S_R = get_genome_set(ord_genome, all_dupl_genome)
#     ord_genome = remove_singletons_in_A_wrt_set_B(ord_genome, S_R)
#     all_dupl_genome = remove_singletons_in_A_wrt_set_B(all_dupl_genome, S_R)
#     return ord_genome, all_dupl_genome


# def get_indexed_sets_for_halving_problem(ord_genome, all_dupl_genome):
#     s_all_genes = create_complete_genes_multiset(ord_genome.get_gene_set() | all_dupl_genome.get_gene_set(), 1)
#     vertex2ind, ind2vertex = enumerate_vertex_multiset(s_all_genes)
#
#     matching, tels = ord_genome.convert_to_matching()
#     B_edges = indexing_graph_edges(matching, vertex2ind)
#     obverse_edges = indexing_obverse_edges_in_n_dupl_genome(s_all_genes, vertex2ind)
#     J_T_B = indexing_vertex_set(tels, vertex2ind)
#
#     J_A = indexing_gene_multiset(create_complete_genes_multiset(all_dupl_genome.get_gene_set(), 1), vertex2ind)
#     J_R = indexing_gene_multiset(create_complete_genes_multiset(ord_genome.get_gene_set(), 1), vertex2ind)
#     J_B = indexing_gene_multiset(ord_genome.counting_dict_of_genes(), vertex2ind)
#
#     ms_all_genes = create_complete_genes_multiset(all_dupl_genome.get_gene_set(), 2)
#     hat_vertex2ind, hat_ind2vertex = enumerate_vertex_multiset(ms_all_genes)
#     bar_J_hat_A = set(hat_vertex2ind.values())
#     a_matching, a_tels = all_dupl_genome.convert_to_matching()
#     A_edges_in_hat_A = indexing_graph_edges(a_matching, hat_vertex2ind)
#     J_T_A_in_hat_A = indexing_vertex_set(a_tels, hat_vertex2ind)
#
#     J_T_A = indexing_vertex_set(a_tels, vertex2ind)
#
#     vert_with_A_edge = set()
#     for u, v in A_edges_in_hat_A:
#         vert_with_A_edge.add(u)
#         vert_with_A_edge.add(v)
#     bar_J_hat_A_compl = bar_J_hat_A - vert_with_A_edge
#
#     equiv_map = define_equiv_function(s_all_genes, 2, vertex2ind, hat_vertex2ind)
#
#     return J_A, J_R, J_B, B_edges, obverse_edges, J_T_A, J_T_B, bar_J_hat_A, A_edges_in_hat_A, J_T_A_in_hat_A, \
#            bar_J_hat_A_compl, equiv_map, vertex2ind, ind2vertex, get_biggest_constant(len(hat_ind2vertex))


