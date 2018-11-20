import logging
import itertools
from functools import reduce
import networkx as nx
import operator

import os

import sys

from genome import parse_genome_in_grimm_file
from impl_gurobi.restricted_versions import create_ilp_formulation_for_restricted_median
from impl_gurobi.utils import create_complete_genes_multiset, remove_singletons_dupl_wrt_gene_set, \
    create_observed_edges_from_gene_multiset, create_vertex_set_from_gene_multiset, get_immediate_subdirectories

logger = logging.getLogger()


class RestrictedMedianConf(object):
    def __init__(self, genomes, name, log_file, tl):
        if len(genomes) != 3:
            raise Exception("Incorrect number of genomes for median problems")

        self.gene_sets = [set(genome.get_gene_multiset().keys()) for genome in genomes]
        ms_all_genes = create_complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()), 1)

        self.bg_vertex_set = create_vertex_set_from_gene_multiset(ms_all_genes)
        self.ind2vertex = [''] + [u for u in self.bg_vertex_set]
        self.vertex2ind = {self.ind2vertex[i]: i for i in range(1, len(self.ind2vertex))}

        obverse_edges = create_observed_edges_from_gene_multiset(ms_all_genes)
        self.ind_obverse_edges = {tuple(sorted((self.vertex2ind[u], self.vertex2ind[v]))) for u, v in obverse_edges}

        genome_graphs = [genome.convert_to_genome_graph() for genome in genomes]
        self.ind_p_i_edges = [{tuple(sorted((self.vertex2ind[u], self.vertex2ind[v])))
                               for u, v in matching} for matching, _ in genome_graphs]
        self.ind_telomers = [{self.vertex2ind[u] for u in telomers} for _, telomers in genome_graphs]
        self.ind_vertex_sets = [{self.vertex2ind[u] for u in
                                 create_vertex_set_from_gene_multiset(create_complete_genes_multiset(gene_set, 1))}
                                for gene_set in self.gene_sets]
        self.ind_median_set = set(reduce(operator.or_,
                                     [self.ind_vertex_sets[i] & self.ind_vertex_sets[j] for i, j in
                                      itertools.combinations(range(len(genomes)), 2)],
                                     set()))

        self.number_of_genomes = len(genomes)
        self.biggest_const = len(self.ind2vertex)
        self.name_model = name
        self.log_file = log_file
        self.time_limit = tl

        flag = False
        for i in range(3):
            if len(self.ind_telomers[i]) != 0:
                flag = True

        self.allowable_telomers = {x: flag for x in self.ind_median_set}

        graph = nx.MultiGraph()
        for edges in self.ind_p_i_edges:
            for u, v in edges:
                graph.add_edge(u, v)

        self.allowable_median_edges = set()
        self.connection_constrs = dict()
        for component in nx.connected_component_subgraphs(graph):
            for v in component.nodes():
                connect_indices = set()
                for u in component.nodes():
                    if u != v:
                        self.allowable_median_edges.add(tuple(sorted((u, v))))
                        connect_indices.add(tuple(sorted((u, v))))
                self.connection_constrs[v] = connect_indices

        self.number_of_cycles = 0
        self.number_of_even_paths = 0

        graph = nx.MultiGraph()
        for u in self.ind_median_set:
            graph.add_node(u)

        for edges in self.ind_p_i_edges[0:2]:
            for u, v in edges:
                graph.add_edge(u, v)

        for component in nx.connected_component_subgraphs(graph):
            isCycle = True

            for v in component.nodes():
                if graph.degree(v) != 2:
                    isCycle = False

            if isCycle:
                self.number_of_cycles += 1
            elif (len(component.nodes()) - 1) % 2 == 0:
                self.number_of_even_paths += 1


def remove_known_singletons(genomes):
    gene_sets = [set(genome.get_gene_multiset().keys()) for genome in genomes]
    median_set = reduce(operator.or_,
                        [gene_sets[i] & gene_sets[j] for i, j in itertools.combinations(range(len(genomes)), 2)],
                        set())
    for i in range(len(genomes)):
        genomes[i] = remove_singletons_dupl_wrt_gene_set(genomes[i], median_set)

    return genomes


def solve_restricted_median(genome_files, out_result_file, out_median_file):
    genomes = [parse_genome_in_grimm_file(genome_file) for genome_file in genome_files]
    genomes = remove_known_singletons(genomes)

    cfg = RestrictedMedianConf(genomes=genomes, name="IGMP", log_file="gurobi_igmp.log", tl=7200)
    answer = create_ilp_formulation_for_restricted_median(cfg)

    answer.write_stats_file(out_result_file)
    answer.write_genome_file(out_median_file)


def dojob(name_directory):
    logging.info('Let us do the ILP')

    for path, name in get_immediate_subdirectories(name_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            logging.info("Working with directory {0}".format(subpath))
            ilppath = os.path.join(subpath, "igmp_ilp")

            genome1 = os.path.join(ilppath, "S1.gen")
            genome2 = os.path.join(ilppath, "S2.gen")
            genome3 = os.path.join(ilppath, "S4.gen")

            result_out_file = os.path.join(ilppath, "result.txt")
            genome_out_file = os.path.join(ilppath, "median.gen")

            solve_restricted_median(genome_files=[genome1, genome2, genome3],
                                    out_result_file=result_out_file,
                                    out_median_file=genome_out_file)


if __name__ == "__main__":
    # if len(sys.argv) < 2:
    #     sys.stderr.write("USAGE: double_distance.py <destination>\n")
    #     sys.exit(1)
    #
    # if os.path.isfile(sys.argv[1]):
    #     sys.stderr.write("<destination> - is directory with specific hierarchy\n")
    #     sys.exit(1)
    #
    # console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
    #                                       "%(message)s", "%H:%M:%S")
    # console_log = logging.StreamHandler()
    # console_log.setFormatter(console_formatter)
    # console_log.setLevel(logging.INFO)
    #
    # log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
    #                                   "%(message)s", "%H:%M:%S")
    #
    # file_log = logging.FileHandler("txt_igmp_median.log", mode="w")
    # file_log.setFormatter(log_formatter)
    #
    # logger.setLevel(logging.INFO)
    # logger.addHandler(console_log)
    # logger.addHandler(file_log)
    #
    # dir_path = os.path.abspath(sys.argv[1])
    # dojob(dir_path)

    # genome1 = """>g1
    #  1 2 -5 -4 -3 @"""
    # genome2 = """>g2
    #  1 -4 -3 -2 5 @"""
    # genome3 = """>g3
    #  -3 -2 -1 4 5 @"""
    #
    solve_restricted_median(["S1.gen", "S2.gen", "S4.gen"], "result.txt", "median.txt")
