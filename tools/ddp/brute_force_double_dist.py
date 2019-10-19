import logging
import argparse
import os
from datetime import datetime

from utils.common import epilog, enable_logging, version
from double_dist import DoubleDistConf
from utils.genome import parse_genome_in_grimm_file
from impl_gurobi.common import remove_singletons_in_ord_wrt_two_dupl, remove_singletons_dupl_wrt_gene_set

logger = logging.getLogger()


def pairs_from_set(v, iterable, n, ans, used=[-1] * 200):
    if v >= n - 1 and used[v] == -1:
        return
    if used[v] > -1 and v < n - 1:
        pairs_from_set(v + 1, iterable, n, ans, used)
        return
    if v >= n - 1:
        res = set()
        for i in range(n):
            res.add(tuple(sorted([iterable[i], iterable[used[i]]])))
        ans.append(res)
    else:
        for i in range(v + 1, n):
            if used[i] < 0:
                used[i] = v
                used[v] = i
                pairs_from_set(v + 1, iterable, n, ans, used)
                used[i] = -1
                used[v] = -1


class Graph:
    def __init__(self, config):
        self.config = config
        self.vertex_num = len(config.bg_vertex_set)
        self.is_A_telomer = [False for u in range(self.vertex_num + 1)]
        for u in config.ind_bg_A_telomers:
            self.is_A_telomer[u] = True
        self.num_telomers_by_color = [0 for c in range(self.vertex_num)]
        self.colors = [-1 for u in range(self.vertex_num + 1)]
        self.edges = {u: [] for u in range(1, self.vertex_num + 1)}

    def add_edges_from_cbg_sets(self, edge_set, mask, edges_num):
        for i in range(edges_num):
            x, y = edge_set[i]
            ux, vx = self.config.equiv_map[x]
            uy, vy = self.config.equiv_map[y]
            if ((1 << i) & mask) > 0:
                self.add_edge(ux, vy)
                self.add_edge(uy, vx)
            else:
                self.add_edge(ux, uy)
                self.add_edge(vy, vx)

    def add_edge(self, u, v):
        self.edges[u].append(v)
        self.edges[v].append(u)

    def assign_colors(self, u, cur_color):
        self.colors[u] = cur_color
        if self.is_A_telomer[u]:
            self.num_telomers_by_color[cur_color] += 1
        for v in self.edges[u]:
            if self.colors[v] == -1:
                self.assign_colors(v, cur_color)

    def add_edges(self, new_edges):
        for u, v in new_edges:
            self.add_edge(u, v)

    def count_obj(self):
        max_color = 0
        for u in range(1, self.vertex_num + 1):
            if self.colors[u] == -1:
                self.assign_colors(u, max_color)
                max_color += 1
        components = max_color
        A_telomer_paths = sum(1 for i in range(max_color) if self.num_telomers_by_color[i] == 2)
        obj = components - A_telomer_paths
        return obj


def solve(config):
    max_obj = 0
    hat_r_set = config.ind_cbg_compl_R
    hat_a_set = config.ind_compl_A
    pairs_hat_r = []
    pairs_from_set(0, list(hat_r_set), len(hat_r_set), pairs_hat_r)
    if len(pairs_hat_r) == 0: pairs_hat_r = [set()]
    for hat_r_pairs in pairs_hat_r:
        x_list = list(config.ind_cbg_R_edges | hat_r_pairs)
        x_num = len(x_list)
        for x_mask in range(1 << x_num):
            pairs_hat_a = []
            pairs_from_set(0, list(hat_a_set), len(hat_a_set), pairs_hat_a)
            if len(pairs_hat_a) == 0: pairs_hat_a = [set()]
            for hat_a_pairs in pairs_hat_a:
                graph = Graph(config)
                graph.add_edges_from_cbg_sets(x_list, x_mask, x_num)
                graph.add_edges(config.ind_bg_A_edges)

                graph.add_edges(hat_a_pairs)
                max_obj = max(max_obj, graph.count_obj())
    return max_obj, config.multiplicity * len(config.s_all_genes) + len(config.ind_cbg_R_telomers) - max_obj


def write_answer(answer, result_out_file, solution_time):
    with open(result_out_file, 'w') as out:
        out.write("# Objective value\n")
        out.write(str(answer[0]) + "\n")
        out.write("# Total DCJ-Indel Distance\n")
        out.write(str(answer[1]) + "\n")
        out.write("# Total time \n")
        out.write(str(solution_time))


def solve_ddp_brute(ord_genome_file, all_dupl_genome_file, out_result_file=None, tl=50, mult=2):
    ord_genome = parse_genome_in_grimm_file(ord_genome_file)
    all_dupl_genome = parse_genome_in_grimm_file(all_dupl_genome_file)

    ord_genome = remove_singletons_in_ord_wrt_two_dupl(ord_genome, all_dupl_genome)
    all_dupl_genome = remove_singletons_dupl_wrt_gene_set(all_dupl_genome, set(ord_genome.get_gene_multiset().keys()))

    config = DoubleDistConf(ordinary_genome=ord_genome, duplicated_genome=all_dupl_genome, name="DDP",
                            log_file="brute_double_dist.log", tl=tl, mult=mult)
    start_time = datetime.now()
    answer = solve(config)
    end_time = datetime.now()
    logging.info(answer)
    if out_result_file:
        write_answer(answer, out_result_file, end_time - start_time)
    return answer


def main():
    parser = argparse.ArgumentParser(description="Brute force solver for the Double Distance Problem (DDP)",
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=epilog())

    parser.add_argument("-r", "--ordinary", dest="ordinary",
                        default="sample/R.txt", metavar="PATH", help="Ordinary genome in GRIMM format")

    parser.add_argument("-t", "--two_dupl", dest="duplicated",
                        default="sample/A.txt", metavar="PATH", help="2-duplicated genome in GRIMM format")

    parser.add_argument("-o", "--out-dir", dest="out_dir",
                        default="sample", required=False,
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

    result_out_file = os.path.join(args.out_dir, "result_brute.txt")

    solve_ddp_brute(ord_genome_file=path_to_ord_genome,
                    all_dupl_genome_file=path_to_dupl_genome,
                    out_result_file=result_out_file)


if __name__ == "__main__":
    main()
