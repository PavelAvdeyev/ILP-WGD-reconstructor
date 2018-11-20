import logging
import itertools
from functools import reduce
import networkx as nx
import operator

import os

import sys
sys.path.append('../')

from genome import parse_genome_in_grimm_file
from impl_gurobi.medians import create_ilp_formulation_for_medians_without_singletons, \
    create_ilp_formulation_for_medians_with_singletons
from impl_gurobi.utils import create_complete_genes_multiset, remove_singletons_dupl_wrt_gene_set, \
    create_observed_edges_from_gene_multiset, create_vertex_set_from_gene_multiset, get_immediate_subdirectories

logger = logging.getLogger()


class MedianConf(object):
    def __init__(self, genomes, name, log_file, tl):
        if len(genomes) != 3:
            raise Exception("Incorrect number of genomes for median problems")

        self.gene_sets = [set(genome.get_gene_multiset().keys()) for genome in genomes]
        self.s_all_genes = create_complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()), 1)

        self.cbg_vertex_set = create_vertex_set_from_gene_multiset(self.s_all_genes)
        self.cbg_ind2vertex = [''] + [u for u in self.cbg_vertex_set]
        self.cbg_vertex2ind = {self.cbg_ind2vertex[i]: i for i in range(1, len(self.cbg_ind2vertex))}

        obverse_edges = create_observed_edges_from_gene_multiset(self.s_all_genes)
        self.ind_cbg_obverse_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in obverse_edges}

        genome_graphs = [genome.convert_to_genome_graph() for genome in genomes]
        self.ind_cbg_p_i_vertex_sets = [{self.cbg_vertex2ind[u] for u in
                                         create_vertex_set_from_gene_multiset(create_complete_genes_multiset(gene_set, 1))}
                                        for gene_set in self.gene_sets]
        self.ind_cbg_p_i_edges = [{tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v])))
                                   for u, v in matching} for matching, _ in genome_graphs]
        self.ind_cbg_p_i_telomers = [{self.cbg_vertex2ind[u] for u in telomers} for _, telomers in genome_graphs]


        self.ind_ancestral_set = set(reduce(operator.or_,
                                            [self.ind_cbg_p_i_vertex_sets[i] & self.ind_cbg_p_i_vertex_sets[j] for i, j in
                                             itertools.combinations(range(len(genomes)), 2)],
                                            set()))
        self.ancestral_gene_set = reduce(operator.or_,
                                         [self.gene_sets[i] & self.gene_sets[j] for i, j in itertools.combinations(range(len(self.gene_sets)), 2)],
                                         set())

        self.number_of_genomes = len(genomes)
        self.biggest_const = len(self.cbg_ind2vertex)
        self.name_model = name
        self.log_file = log_file
        self.time_limit = tl


class ClassicMedian(MedianConf):
    def __init__(self, genomes, name, log_file, tl):
        super().__init__(genomes, name, log_file, tl)
        self.allowable_ancestral_edges = {tuple(sorted([u, v]))
                                       for u, v in itertools.combinations(self.ind_ancestral_set, 2)}

        flag = False
        for i in range(self.number_of_genomes):
            if len(self.ind_cbg_p_i_telomers[i]) != 0:
                flag = True

        self.allowable_ancestral_telomers = {x: flag for x in self.ind_ancestral_set}
        self.connection_ancestral_constrs = {u: {tuple(sorted((u, v))) for v in self.ind_ancestral_set if u != v}
                                   for u in self.ind_ancestral_set}


class ConservedMedian(MedianConf):
    def __init__(self, genomes, name, log_file, tl):
        super().__init__(genomes, name, log_file, tl)

        set_all_telomers = reduce(operator.or_, self.ind_cbg_p_i_telomers, set())
        self.allowable_ancestral_telomers = {x: True if x in set_all_telomers else False for x in self.ind_ancestral_set}

        graph = nx.MultiGraph()
        for edges in self.ind_cbg_p_i_edges:
            for u, v in edges:
                graph.add_edge(u, v)

        self.allowable_ancestral_edges = set()
        self.connection_ancestral_constrs = dict()

        for v in self.ind_ancestral_set:
            connect_indices = set()

            if v in graph:
                for u in graph.neighbors(v):
                    if u != v and u in self.ind_ancestral_set:
                        self.allowable_ancestral_edges.add(tuple(sorted((u, v))))
                        connect_indices.add(tuple(sorted((u, v))))

            self.connection_ancestral_constrs[v] = connect_indices


def remove_known_singletons(genomes):
    gene_sets = [set(genome.get_gene_multiset().keys()) for genome in genomes]
    median_set = reduce(operator.or_,
                        [gene_sets[i] & gene_sets[j] for i, j in itertools.combinations(range(len(genomes)), 2)],
                        set())
    for i in range(len(genomes)):
        genomes[i] = remove_singletons_dupl_wrt_gene_set(genomes[i], median_set)

    return genomes


def medians_with_singletons(genome_files, out_result_file, out_median_file, problem):
    genomes = [parse_genome_in_grimm_file(genome_file) for genome_file in genome_files]
    genomes = remove_known_singletons(genomes)

    if problem == "CGMP":
        cfg = ConservedMedian(genomes=genomes, name="CGMP", log_file="gurobi_cgmp.log", tl=7200)
    else:
        cfg = ClassicMedian(genomes=genomes, name="GMP", log_file="gurobi_gmp.log", tl=7200)

    answer = create_ilp_formulation_for_medians_with_singletons(cfg)

    answer.write_stats_file(out_result_file)
    answer.write_genome_file(out_median_file)


def medians_without_singletons(genome_files, out_result_file, out_median_file, problem):
    genomes = [parse_genome_in_grimm_file(genome_file) for genome_file in genome_files]
    genomes = remove_known_singletons(genomes)

    if problem == "CGMP":
        cfg = ConservedMedian(genomes=genomes, name="CGMP", log_file="gurobi_cgmp.log", tl=7200)
    else:
        cfg = ClassicMedian(genomes=genomes, name="GMP", log_file="gurobi_gmp.log", tl=7200)

    answer = create_ilp_formulation_for_medians_without_singletons(cfg)

    answer.write_stats_file(out_result_file)
    answer.write_genome_file(out_median_file)


def dojob(name_directory):
    logging.info('Let us do the ILP')

    for path, name in get_immediate_subdirectories(name_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            logging.info("Working with directory {0}".format(subpath))
            ilppath = os.path.join(subpath, "heur_ilp")

            genome1 = os.path.join(ilppath, "S1.gen")
            genome2 = os.path.join(ilppath, "S2.gen")
            genome3 = os.path.join(ilppath, "S4.gen")

            result_out_file = os.path.join(ilppath, "result.txt")
            genome_out_file = os.path.join(ilppath, "median.gen")

            medians_without_singletons(genome_files=[genome1, genome2, genome3],
                                       out_result_file=result_out_file,
                                       out_median_file=genome_out_file,
                                       problem="CGMP")


if __name__ == "__main__":
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

    file_log = logging.FileHandler("txt_heur_median.log", mode="w")
    file_log.setFormatter(log_formatter)

    logger.setLevel(logging.INFO)
    logger.addHandler(console_log)
    # logger.addHandler(file_log)

    # dir_path = os.path.abspath(sys.argv[1])
    # dojob(dir_path)

    # genome1 = """>g1
    #      1 2 -5 -4 -3 @"""
    # genome2 = """>g2
    #      1 -4 -3 -2 5 @"""
    # genome3 = """>g3
    #      -3 -2 -1 4 5 @"""

    # '''
    # genome1 = """>g1
    #  1 2 -3 $"""
    # genome2 = """>g2
    #  1 -2 3 $"""
    # genome3 = """>g3
    #  -1 2 3 $"""
    # '''
    #
    # '''
    # genome1 = """>g1
    # 2 -3 @"""
    # genome2 = """>g2
    # 1 3 @"""
    # genome3 = """>g3
    # -1 2 @"""
    # '''

    medians_without_singletons(["S1.gen", "S2.gen", "S4.gen"], "result.txt", "median.txt", "CGMP")