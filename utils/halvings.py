import logging

import itertools
from functools import reduce
import networkx as nx
import operator

from utils.genome import parse_genome_in_grimm_file
from impl_gurobi.halvings import create_ilp_formulation_for_halvings_without_singletons
from utils.set_definer import complete_genes_multiset, vertex_set_from_gene_multiset, \
    observed_edges_from_gene_multiset, define_equiv_function

logger = logging.getLogger()


class HalvingConf(object):
    def __init__(self, duplicated_genome, ordinary_genomes, name, log_file, tl, mult=2):
        if len(ordinary_genomes) != 1:
            raise Exception("Incorrect number of genomes for guided halving problems")

        if mult < 2 or mult > 3:
            raise Exception("Unsupported multiplication of a genome")

        self.gene_sets = [set(genome.get_gene_multiset().keys()) for genome in ordinary_genomes]
        self.genes_of_dupl_genome = duplicated_genome.get_gene_multiset()

        self.s_all_genes = complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()) |
                                                   self.genes_of_dupl_genome.keys(), 1)
        # Coding contracted breakpoint graph
        self.cbg_vertex_set = vertex_set_from_gene_multiset(self.s_all_genes)
        self.cbg_ind2vertex = [''] + [u for u in self.cbg_vertex_set]
        self.cbg_vertex2ind = {self.cbg_ind2vertex[i]: i for i in range(1, len(self.cbg_ind2vertex))}

        obverse_edges = observed_edges_from_gene_multiset(self.s_all_genes)
        self.ind_cbg_obverse_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v])))
                                      for u, v in obverse_edges}

        genome_graphs = [genome.convert_to_genome_graph() for genome in ordinary_genomes]
        self.ind_cbg_p_i_vertex_sets = [{self.cbg_vertex2ind[u] for u in
                                         vertex_set_from_gene_multiset(complete_genes_multiset(gene_set, 1))}
                                        for gene_set in self.gene_sets]
        self.ind_cbg_p_i_edges = [{tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in matching}
                                  for matching, _ in genome_graphs]
        self.ind_cbg_p_i_telomers = [{self.cbg_vertex2ind[u] for u in telomers} for _, telomers in genome_graphs]

        # This contracted genome graph does not contain parallel edges. Maybe a problem.
        cbg_A_matching, cbg_A_telomers = duplicated_genome.convert_to_contracted_genome_graph()
        self.ind_cbg_A_vertices = {self.cbg_vertex2ind[u] for u in
                                   vertex_set_from_gene_multiset(
                                       complete_genes_multiset(self.genes_of_dupl_genome.keys(), 1))}
        self.ind_cbg_A_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in
                                cbg_A_matching}
        self.ind_cbg_A_telomers = {self.cbg_vertex2ind[u] for u in cbg_A_telomers}

        self.ms_all_genes = complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()) |
                                                    self.genes_of_dupl_genome.keys(), mult)
        # Coding breakpoint graph
        self.bg_vertex_set = vertex_set_from_gene_multiset(self.ms_all_genes)
        self.bg_ind2vertex = [''] + [u for u in self.bg_vertex_set]
        self.bg_vertex2ind = {self.bg_ind2vertex[i]: i for i in range(1, len(self.bg_ind2vertex))}

        bg_A_matching, bg_A_telomers = duplicated_genome.convert_to_genome_graph()
        self.ind_bg_A_vertices = {self.bg_vertex2ind[u] for u in
                                  vertex_set_from_gene_multiset(self.genes_of_dupl_genome)}
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
        self.allowable_ancestral_telomers = {x: True if x in set_all_telomers else False for x in
                                             self.ind_ancestral_set}

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


def halvings_without_singletons(ordinary_genome_file, all_dupl_genome_file, out_result_file, out_predup_file,
                                problem, gurobi_log_file, time_limit):
    logging.info('Start to solve {0} (excluding singletons) with time limit equals {1}'.format(problem, time_limit))

    ord_genomes = [parse_genome_in_grimm_file(ordinary_genome_file)]
    all_dupl_genome = parse_genome_in_grimm_file(all_dupl_genome_file)

    logging.info('Create ILP config')
    if problem == "CGGHP":
        cfg = ConservedHalving(duplicated_genome=all_dupl_genome,
                               ordinary_genomes=ord_genomes,
                               name="CGGHP",
                               log_file=gurobi_log_file,
                               tl=time_limit,
                               mult=2)
    else:
        cfg = ClassicHalving(duplicated_genome=all_dupl_genome,
                             ordinary_genomes=ord_genomes,
                             name="GGHP",
                             log_file=gurobi_log_file,
                             tl=time_limit,
                             mult=2)

    answer = create_ilp_formulation_for_halvings_without_singletons(cfg)

    if answer is not None:
        logging.info('Save results.')
        answer.write_stats_file(out_result_file)
        answer.write_genome_file(out_predup_file)
    else:
        logging.info('There are no answers. Please, check log file.')
