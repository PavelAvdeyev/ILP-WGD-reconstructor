import argparse
import logging
import itertools
from functools import reduce
import networkx as nx
import operator
import os

from utils.common import epilog, version, enable_logging, get_immediate_subdirectories
from utils.genome import parse_genome_in_grimm_file
from impl_gurobi.restricted_median import create_ilp_formulation_for_restricted_median
from utils.set_definer import complete_genes_multiset, observed_edges_from_gene_multiset, \
    vertex_set_from_gene_multiset

logger = logging.getLogger()


class RestrictedMedianConf(object):
    def __init__(self, genomes, name, log_file, tl):
        if len(genomes) != 3:
            raise Exception("Incorrect number of genomes for median problem")

        self.gene_sets = [set(genome.get_gene_multiset().keys()) for genome in genomes]
        self.s_all_genes = complete_genes_multiset(reduce(operator.or_, self.gene_sets, set()), 1)

        self.cbg_vertex_set = vertex_set_from_gene_multiset(self.s_all_genes)
        self.cbg_ind2vertex = [''] + [u for u in self.cbg_vertex_set]
        self.cbg_vertex2ind = {self.cbg_ind2vertex[i]: i for i in range(1, len(self.cbg_ind2vertex))}

        obverse_edges = observed_edges_from_gene_multiset(self.s_all_genes)
        self.ind_cbg_obverse_edges = {tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v]))) for u, v in
                                      obverse_edges}

        genome_graphs = [genome.convert_to_genome_graph() for genome in genomes]
        self.ind_cbg_p_i_vertex_sets = [{self.cbg_vertex2ind[u] for u in
                                         vertex_set_from_gene_multiset(
                                             complete_genes_multiset(gene_set, 1))}
                                        for gene_set in self.gene_sets]
        self.ind_cbg_p_i_edges = [{tuple(sorted((self.cbg_vertex2ind[u], self.cbg_vertex2ind[v])))
                                   for u, v in matching} for matching, _ in genome_graphs]
        self.ind_cbg_p_i_telomers = [{self.cbg_vertex2ind[u] for u in telomers} for _, telomers in genome_graphs]

        self.ind_ancestral_set = set(reduce(operator.or_,
                                            [self.ind_cbg_p_i_vertex_sets[i] & self.ind_cbg_p_i_vertex_sets[j]
                                             for i, j in itertools.combinations(range(len(genomes)), 2)],
                                            set()))

        self.ancestral_gene_set = reduce(operator.or_,
                                         [self.gene_sets[i] & self.gene_sets[j] for i, j in
                                          itertools.combinations(range(len(self.gene_sets)), 2)],
                                         set())

        self.number_of_genomes = len(genomes)
        self.biggest_const = len(self.cbg_ind2vertex)
        self.name_model = name
        self.log_file = log_file
        self.time_limit = tl

        flag = False
        for i in range(3):
            if len(self.ind_cbg_p_i_telomers[i]) != 0:
                flag = True

        self.allowable_telomers = {x: flag for x in self.ind_ancestral_set}

        graph = nx.MultiGraph()
        graph.add_nodes_from(self.cbg_vertex2ind.values())
        for edges in self.ind_cbg_p_i_edges[0:2]:
            for u, v in edges:
                graph.add_edge(u, v)

        self.allowable_ancestral_edges = set()
        self.connection_constrs = dict()
        for component in nx.connected_component_subgraphs(graph):
            for v in component.nodes():
                connect_indices = set()
                for u in component.nodes():
                    if u != v:
                        self.allowable_ancestral_edges.add(tuple(sorted((u, v))))
                        connect_indices.add(tuple(sorted((u, v))))
                self.connection_constrs[v] = connect_indices

        self.number_of_cycles = 0
        self.number_of_even_paths = 0

        for component in nx.connected_component_subgraphs(graph):
            isCycle = True

            for v in component.nodes():
                if graph.degree(v) != 2:
                    isCycle = False

            if isCycle:
                self.number_of_cycles += 1
            elif (len(component.nodes()) - 1) % 2 == 0:
                self.number_of_even_paths += 1


def solve_restricted_median_problem(genome_files, out_result_file, out_median_file, gurobi_log_file, time_limit=7200):
    logging.info('Start to solve IGMP problem with time limit equals {0}'.format(time_limit))

    genomes = [parse_genome_in_grimm_file(genome_file) for genome_file in genome_files]

    logging.info('Create ILP config')
    cfg = RestrictedMedianConf(genomes=genomes, name="IGMP", log_file=gurobi_log_file, tl=time_limit)

    answer = create_ilp_formulation_for_restricted_median(cfg)

    if answer is not None:
        logging.info('Save results.')
        answer.write_stats_file(out_result_file)
        answer.write_genome_file(out_median_file)
    else:
        logging.info('There are no answers. Please, check log file.')


def do_test_job(input_directory, gurobi_log_file):
    logging.info('Let\'s do comparison test for the IGMP-ILP. Input directory is {0}'.format(input_directory))

    for path, name in get_immediate_subdirectories(input_directory):
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

            solve_restricted_median_problem(genome_files=[genome1, genome2, genome3],
                                            out_result_file=result_out_file,
                                            out_median_file=genome_out_file,
                                            gurobi_log_file=gurobi_log_file)


def test_run():
    parser = argparse.ArgumentParser(description="ILP solver for the Intermediate Genome Median Problem (IGMP)",
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

    args.log_file = os.path.join(args.paper_test, "igmp.log")
    args.gurobi_log_file = os.path.join(args.paper_test, "gurobi_igmp.log")
    enable_logging(args.log_file, overwrite=False)

    do_test_job(args.paper_test, args.gurobi_log_file)


def main():
    parser = argparse.ArgumentParser(description="ILP solver for the Intermediate Genome Median Problem (IGMP)",
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

    args.log_file = os.path.join(args.out_dir, "igmp.log")
    args.gurobi_log_file = os.path.join(args.out_dir, "gurobi_igmp.log")
    enable_logging(args.log_file, overwrite=False)

    path_to_genomes = [os.path.abspath(f) for f in args.files]
    result_out_file = os.path.join(args.out_dir, "result.txt")
    genome_out_file = os.path.join(args.out_dir, "median.gen")

    solve_restricted_median_problem(genome_files=path_to_genomes,
                                    out_result_file=result_out_file,
                                    out_median_file=genome_out_file,
                                    gurobi_log_file=args.gurobi_log_file,
                                    time_limit=args.time_limit)


if __name__ == "__main__":
    main()
