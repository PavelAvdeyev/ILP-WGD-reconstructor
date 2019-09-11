import logging
import itertools
from functools import reduce
import networkx as nx
import operator

from utils.genome import parse_genome_in_grimm_file
from impl_gurobi.medians import create_ilp_formulation_for_medians_without_singletons, \
    create_ilp_formulation_for_medians_with_singletons
from impl_gurobi.common import complete_genes_multiset, remove_singletons_dupl_wrt_gene_set, \
    observed_edges_from_gene_multiset, vertex_set_from_gene_multiset

logger = logging.getLogger()


class MedianConf(object):
    def __init__(self, genomes, name, log_file, tl):
        if len(genomes) != 3:
            raise Exception("Incorrect number of genomes for median problems")

        self.gene_sets = [set(genome.get_gene_multiset().keys()) for genome in genomes]
        self.s_all_genes = complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()), 1)

        self.cbg_vertex_set = vertex_set_from_gene_multiset(self.s_all_genes)
        self.cbg_ind2vertex = [''] + [u for u in self.cbg_vertex_set]
        self.cbg_vertex2ind = {self.cbg_ind2vertex[i]: i for i in range(1, len(self.cbg_ind2vertex))}

        obverse_edges = observed_edges_from_gene_multiset(self.s_all_genes)
        self.ind_cbg_obverse_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in obverse_edges}

        genome_graphs = [genome.convert_to_genome_graph() for genome in genomes]
        self.ind_cbg_p_i_vertex_sets = [{self.cbg_vertex2ind[u] for u in
                                         vertex_set_from_gene_multiset(complete_genes_multiset(gene_set, 1))}
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


def medians_with_singletons(genome_files, out_result_file, out_median_file, problem, gurobi_log_file, time_limit):
    logging.info('Start to solve {0} (including singletons) with time limit equals {1}'.format(problem, time_limit))

    genomes = [parse_genome_in_grimm_file(genome_file) for genome_file in genome_files]
    genomes = remove_known_singletons(genomes)

    logging.info('Create ILP config')
    if problem == "CGMP":
        cfg = ConservedMedian(genomes=genomes, name="CGMP", log_file=gurobi_log_file, tl=time_limit)
    else:
        cfg = ClassicMedian(genomes=genomes, name="GMP", log_file=gurobi_log_file, tl=time_limit)

    answer = create_ilp_formulation_for_medians_with_singletons(cfg)

    if answer is not None:
        logging.info('Save results.')
        answer.write_stats_file(out_result_file)
        answer.write_genome_file(out_median_file)
    else:
        logging.info('There are no answers. Please, check log file.')


def medians_without_singletons(genome_files, out_result_file, out_median_file, problem, gurobi_log_file, time_limit):
    logging.info('Start to solve {0} (excluding singletons) with time limit equals {1}'.format(problem, time_limit))
    genomes = [parse_genome_in_grimm_file(genome_file) for genome_file in genome_files]

    logging.info('Create ILP config')
    if problem == "CGMP":
        cfg = ConservedMedian(genomes=genomes, name="CGMP", log_file=gurobi_log_file, tl=time_limit)
    else:
        cfg = ClassicMedian(genomes=genomes, name="GMP", log_file=gurobi_log_file, tl=time_limit)

    answer = create_ilp_formulation_for_medians_without_singletons(cfg)

    if answer is not None:
        logging.info('Save results.')
        answer.write_stats_file(out_result_file)
        answer.write_genome_file(out_median_file)
    else:
        logging.info('There are no answers. Please, check log file.')


