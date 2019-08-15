import argparse
import logging

import os

import operator

import itertools
import networkx as nx

from impl_gurobi.common import get_immediate_subdirectories, create_complete_genes_multiset, \
    create_observed_edges_from_gene_multiset, create_vertex_set_from_gene_multiset, define_equiv_function
from impl_gurobi.restricted_halving import create_ilp_formulation_for_restricted_halving
from utils.common import epilog, enable_logging, version
from utils.genome import parse_genome_in_grimm_file
from utils.halvings import halvings_without_singletons


class RestrictedHalvingConf(object):
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

        flag = False
        if len(self.ind_bg_A_telomers) != 0:
            flag = True

        for i in range(self.number_of_genomes):
            if len(self.ind_cbg_p_i_telomers[i]) != 0:
                flag = True
        self.allowable_ancestral_telomers = {x: flag for x in self.ind_ancestral_set}

        graph = nx.MultiGraph()
        graph.add_nodes_from(self.cbg_vertex2ind.values())
        for edges in self.ind_cbg_A_edges:
            for u, v in edges:
                graph.add_edge(u, v)

        self.number_of_even_cycles = 0
        self.number_of_even_paths = 0
        vertex_sets_for_ancestral_edges = []
        meta_vertex_set_for_edges = []
        for component in nx.connected_component_subgraphs(graph):
            isCycle = True

            for v in component.nodes():
                if graph.degree(v) != 2:
                    isCycle = False

            if isCycle and len(component.nodes()) % 2 == 0:
                self.number_of_even_cycles += 1
                vertex_sets_for_ancestral_edges.append(list(component.nodes()))
            elif len(component.nodes()) % 2 != 0 and not isCycle:
                self.number_of_even_paths += 1
                meta_vertex_set_for_edges.extend(component.nodes())
            else:
                meta_vertex_set_for_edges.extend(component.nodes())
        vertex_sets_for_ancestral_edges.append(meta_vertex_set_for_edges)

        self.allowable_ancestral_edges = set()
        self.connection_ancestral_constrs = dict()
        for v_set in vertex_sets_for_ancestral_edges:
            self.allowable_ancestral_edges.update({tuple(sorted([u, v])) for u, v in itertools.combinations(v_set, 2)})
            self.connection_ancestral_constrs.update({u: {tuple(sorted((u, v))) for v in v_set if u != v} for u in v_set})


def solve_restricted_halving_problem(ordinary_genome_file, all_dupl_genome_file, out_result_file, out_predup_file,
                                     gurobi_log_file, time_limit):
    logging.info('Start to solve RGGHP problem with time limit equals {0}'.format(time_limit))

    ord_genomes = [parse_genome_in_grimm_file(ordinary_genome_file)]
    all_dupl_genome = parse_genome_in_grimm_file(all_dupl_genome_file)

    logging.info('Create ILP config')
    cfg = RestrictedHalvingConf(duplicated_genome=all_dupl_genome,
                                ordinary_genomes=ord_genomes,
                                name="RGGHP",
                                log_file=gurobi_log_file,
                                tl=time_limit,
                                mult=2)

    answer = create_ilp_formulation_for_restricted_halving(cfg)

    if answer is not None:
        logging.info('Save results.')
        answer.write_stats_file(out_result_file)
        answer.write_genome_file(out_predup_file)
    else:
        logging.info('There are no answers. Please, check log file.')


def do_test_job(input_directory, gurobi_log_file):
    logging.info('Let\'s do comparison test for the RGGHP-ILP. Input directory is {0}'.format(input_directory))

    for path, name in get_immediate_subdirectories(input_directory):
        for subpath, subname in get_immediate_subdirectories(path):
            if os.path.isfile(subpath):
                continue

            logging.info("Working with directory {0}".format(subpath))
            ilppath = os.path.join(subpath, "rgghp_ilp")

            ordinary_genome_file = os.path.join(ilppath, "B.gen")
            all_dupl_genome_file = os.path.join(ilppath, "A.gen")

            result_out_file = os.path.join(ilppath, "result.txt")
            genome_out_file = os.path.join(ilppath, "pre_dup.gen")

            # TODO
            halvings_without_singletons(ordinary_genome_file=ordinary_genome_file,
                                        all_dupl_genome_file=all_dupl_genome_file,
                                        out_result_file=result_out_file,
                                        out_predup_file=genome_out_file,
                                        problem="RGGHP",
                                        gurobi_log_file=gurobi_log_file,
                                        time_limit=7200
                                        )


def test_run():
    parser = argparse.ArgumentParser(description="ILP solver for the Restricted Guided Genome Halving Problem (RGGHP)",
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

    args.log_file = os.path.join(args.paper_test, "rgghp.log")
    args.gurobi_log_file = os.path.join(args.paper_test, "gurobi_rgghp.log")
    enable_logging(args.log_file, overwrite=False)

    do_test_job(args.paper_test, args.gurobi_log_file)


def main():
    parser = argparse.ArgumentParser(description="ILP solver for the Restricted Guided Genome Halving Problem (RGGHP)",
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

    args.log_file = os.path.join(args.out_dir, "rgghp.log")
    args.gurobi_log_file = os.path.join(args.out_dir, "gurobi_rgghp.log")
    enable_logging(args.log_file, overwrite=False)

    ordinary_genome_file = os.path.abspath(args.ordinary)
    all_dupl_genome_file = os.path.abspath(args.duplicated)

    result_out_file = os.path.join(args.out_dir, "result.txt")
    genome_out_file = os.path.join(args.out_dir, "pre_dup.gen")

    halvings_without_singletons(ordinary_genome_file=ordinary_genome_file,
                                all_dupl_genome_file=all_dupl_genome_file,
                                out_result_file=result_out_file,
                                out_predup_file=genome_out_file,
                                problem="RGGHP",
                                gurobi_log_file=args.gurobi_log_file,
                                time_limit=args.time_limit)


if __name__ == "__main__":
    main()
