import os

import itertools

from utils.genome import Genome, get_extremity_by_gene
import bg.breakpoint_graph


class ILPAnswer(object):
    def __init__(self, ov=0, score=0, es=0, genome=None, tm=0):
        self.obj_val = ov
        self.score = score
        self.exit_status = es
        self.genome = genome
        self.solution_time = tm

    def update_score_by_singletons(self, value):
        self.score += value

    def write_stats_file(self, result_out_file):
        with open(result_out_file, 'w') as out:
            out.write("# Objective value\n")
            out.write(str(self.obj_val) + "\n")
            out.write("# Total DCJ-Indel Distance\n")
            out.write(str(self.score) + "\n")
            out.write("# Is it solved \n")
            out.write(str(self.exit_status) + "\n")

    def write_genome_file(self, out_genome_file):
        with open(out_genome_file, 'w') as out:
            out.write(">answer\n")
            for chr_ind, chromosome in self.genome.items():
                for chr_type, blocks in chromosome:
                    string = " ".join(value if sign == "+" else (sign + value) for sign, value in blocks)
                    string += " {chr_type}".format(chr_type=chr_type)
                    out.write(string + "\n")


def get_genome_graph_from_vars(rs, r_dot, gene_set, edge_set, telomer_set, ind2vertex):
    result_graph = bg.breakpoint_graph.BreakpointGraph()

    # graph = nx.Graph()
    #
    # for gene in gene_set:
    #     graph.add_node(get_extremity_by_gene(gene, True, 1)[:-2])
    #     graph.add_node(get_extremity_by_gene(gene, False, 1)[:-2])

    name2vertices = {}
    for gene in gene_set:
        left, right = bg.vertices.TaggedBlockVertex(get_extremity_by_gene(gene, True, 1)[:-2]), \
                      bg.vertices.TaggedBlockVertex(get_extremity_by_gene(gene, False, 1)[:-2])
        left.mate_vertex = right
        right.mate_vertex = left
        name2vertices[left.name] = left
        name2vertices[right.name] = right

    # text_multiset = defaultdict(list)
    # text_set = set()
    # for u, v in edge_set:
    #     i, j = tuple(sorted([u, v]))
    #     if rs[i, j].x == 1:
    #         if ind2vertex[j][:-2] == ind2vertex[i][:-2]:
    #             print("ALLERT ALLERT")
    #         text_multiset[ind2vertex[i][:-2]].append(ind2vertex[j][:-2])
    #         text_multiset[ind2vertex[j][:-2]].append(ind2vertex[i][:-2])
    #         text_set.add(ind2vertex[i][:-2])
    #         text_set.add(ind2vertex[j][:-2])
    #         graph.add_edge(ind2vertex[i][:-2], ind2vertex[j][:-2])

    # for a, b in text_multiset.items():
    #     if len(b) != 1:
    #         print(a)
    #         print(b)

    # text_set2 = set()
    # for u, cond in telomer_set.items():
    #     if cond and r_dot[u].x == 1:
    #         text_set2.add(ind2vertex[u][:-2])
    #         if ind2vertex[u][:-2] in text_set:
    #             print("ALERT! ALERT! ALERT!")
    #         graph.add_node(ind2vertex[u][:-2])
    #
    # print("----MISSING-----")
    # for v in graph:
    #     if len(list(graph.neighbors(v))) == 0:
    #         if v not in text_set2:
    #             print("Alert!")
    #             print(v)

    for u, v in edge_set:
        i, j = tuple(sorted([u, v]))
        if rs[i, j].x == 1:
            result_graph.add_edge(name2vertices[ind2vertex[i][:-2]],
                                  name2vertices[ind2vertex[j][:-2]],
                                  bg.multicolor.Multicolor(1), merge=False)

    for u, cond in telomer_set.items():
        if cond and r_dot[u].x == 1:
            result_graph.add_edge(name2vertices[ind2vertex[u][:-2]],
                                  bg.vertices.TaggedInfinityVertex(ind2vertex[u][:-2]),
                                  bg.multicolor.Multicolor(1), merge=False)

    return result_graph.get_blocks_order()


def complete_genes_multiset(gene_set, copies):
    complete_multiset = dict()
    for gene in gene_set:
        complete_multiset[gene] = copies
    return complete_multiset


def vertex_set_from_gene_multiset(gene_multiset):
    vertex_set = set()
    for gene, copies in gene_multiset.items():
        for left in [True, False]:
            for i in range(1, copies + 1):
                vertex_set.add(get_extremity_by_gene(gene, left, i))
    return vertex_set


def observed_edges_from_gene_multiset(gene_multiset):
    obverse_edges = set()
    for gene, copies in gene_multiset.items():
        for i in range(1, copies + 1):
            obverse_edges.add((get_extremity_by_gene(gene, True, i), get_extremity_by_gene(gene, False, i)))
    return obverse_edges


def product_set(vertex_set1, vertex_set2):
    return set(tuple(sorted([u, v])) for u, v in itertools.product(vertex_set1, vertex_set2) if u != v)


def general_allowable_set(vertex_set):
    return {tuple(sorted([u, v])) for u, v in itertools.combinations(vertex_set, 2)}


def general_conditional_set(vertex_set):
    return {u: {tuple(sorted((u, v))) for v in vertex_set if u != v} for u in vertex_set}


def conserved_allowable_set(vertex_set, graph, telomers):
    allowable_edge_set = set()
    for v in vertex_set:
        for u1, u2 in graph.edges(v):
            if u1 in vertex_set and u2 in vertex_set:
                allowable_edge_set.add(tuple(sorted([u1, u2])))
    allowable_telomer_set = {v for v in vertex_set if v in telomers}
    return allowable_edge_set, allowable_telomer_set


def define_equiv_function(gene_multiset, vertex2ind, bar_vertex2ind):
    equiv = {}
    for gene, copies in gene_multiset.items():
        for left in [True, False]:
            temp = []
            for i in range(1, copies + 1):
                temp.append(bar_vertex2ind[get_extremity_by_gene(gene, left, i)])
            equiv[vertex2ind[get_extremity_by_gene(gene, left, 1)]] = temp
    return equiv


def remove_singletons_in_ord_wrt_two_dupl(ord_genome, dupl_genome):
    new_genome = Genome(ord_genome.get_name())
    dupl_gene_set = set(dupl_genome.get_gene_multiset().keys())

    matching, _ = dupl_genome.convert_to_genome_graph()
    matching = set(matching)

    def get_partitioned_gene_set(gene_multiset):
        s_1, s_2 = set(), set()
        for gene, count in gene_multiset.items():
            if count == 1:
                s_1.add(gene)
            else:
                s_2.add(gene)
        return s_1, s_2

    s_1_all_dupl, s_2_all_dupl = get_partitioned_gene_set(dupl_genome.get_gene_multiset())

    for chromosome in ord_genome:
        if chromosome.is_circular():
            s_1_chr, s_2_chr = get_partitioned_gene_set(chromosome.get_gene_multiset())

            if not len(s_2_chr) and s_1_chr.intersection(s_1_all_dupl) == s_1_chr:
                chromosome_matching = []
                chromosome.convert_to_genome_graph(edges=chromosome_matching)
                chromosome_matching, _ = chromosome.convert_to_matching()

                flag = False
                for u, v in chromosome_matching:
                    if (u, v) not in matching and (v, u) not in matching:
                        flag = True

                if flag:
                    new_genome.append(chromosome)

            elif len(chromosome.get_gene_set().intersection(dupl_gene_set)):
                new_genome.append(chromosome)
        else:
            new_genome.append(chromosome)

    return new_genome

#######################################################
#
#
# BEGINNING_INDEX = 1
#
#
# def complete_genes_multiset(gene_set, copies):
#     complete_multiset = {}
#     for gene in gene_set:
#         complete_multiset[gene] = copies
#     return complete_multiset
#
#
# def enumerate_vertex_multiset(multiset):
#     cbg_vertex2ind, cbg_ind2vertex = {}, {}
#     max_size = BEGINNING_INDEX
#
#     for gene, copies in multiset.items():
#         for left in [True, False]:
#             for i in range(1, copies + 1):
#                 cbg_vertex2ind[get_extremity_by_gene(gene, left, i)] = max_size
#                 cbg_ind2vertex[max_size] = get_extremity_by_gene(gene, left, i)
#                 max_size += 1
#
#     return cbg_vertex2ind, cbg_ind2vertex
#
#
# def indexing_obverse_edges_in_n_dupl_genome(all_genes, cbg_vertex2ind):
#     obverse_edges = set()
#     for gene, copies in all_genes.items():
#         for i in range(1, copies + 1):
#             obverse_edges.add(tuple(sorted([cbg_vertex2ind[get_extremity_by_gene(gene, True, i)],
#                                             cbg_vertex2ind[get_extremity_by_gene(gene, False, i)]])))
#     return obverse_edges
#
#
# def indexing_graph_edges(matching, cbg_vertex2ind):
#     indexing_edges = set()
#     for u, v in matching:
#         indexing_edges.add(tuple(sorted([cbg_vertex2ind[u], cbg_vertex2ind[v]])))
#     return indexing_edges
#
#
# def indexing_vertex_set(vertex_set, cbg_vertex2ind):
#     indexing_v_set = set()
#     for u in vertex_set:
#         indexing_v_set.add(cbg_vertex2ind[u])
#     return indexing_v_set
#
#
# def indexing_gene_multiset(gene_multiset, cbg_vertex2ind):
#     indexing_v_set = set()
#     for gene, copies in gene_multiset.items():
#         for i in range(1, copies + 1):
#             for left in [True, False]:
#                 indexing_v_set.add(cbg_vertex2ind[get_extremity_by_gene(gene, left, i)])
#     return indexing_v_set
#
#
# def get_biggest_constant(number_of_vertices):
#     return BEGINNING_INDEX + number_of_vertices
#
#
# def define_equiv_function(gene_set, copies, cbg_vertex2ind, bar_vertex2ind):
#     equiv = {}
#     for gene in gene_set:
#         for left in [True, False]:
#             temp = []
#             for i in range(1, copies + 1):
#                 temp.append(bar_vertex2ind[get_extremity_by_gene(gene, left, i)])
#             equiv[cbg_vertex2ind[get_extremity_by_gene(gene, left, 1)]] = temp
#     return equiv
#
# def write_stats_file(result_out_file, obj_val, score, status):
#     with open(result_out_file, 'w') as out:
#         out.write("# Number of cycles\n")
#         out.write(str(obj_val) + "\n")
#         out.write("# Total DCJ-Indel Distance\n")
#         out.write(str(score) + "\n")
#         out.write("# Is it solved \n")
#         out.write(str(status) + "\n")
#         out.write(str(1) + "\n")
#
#
# def write_genome_file(out_genome_file, genome):
#     with open(out_genome_file, 'w') as out:
#         out.write(">answer\n")
#         for chr_ind, chromosome in genome.items():
#             for chr_type, blocks in chromosome:
#                 string = " ".join(value if sign == "+" else (sign + value) for sign, value in blocks)
#                 string += " {chr_type}".format(chr_type=chr_type)
#                 out.write(string + "\n")
#
#
# def create_gene_sets_for_pair_genomes(ordinary_genome, all_dupl_genome):
#     S_B = ordinary_genome.get_gene_set()
#     S_1_A, S_2_A = all_dupl_genome.get_partitioned_gene_set()
#     S_A = S_1_A | S_2_A
#     S_R = (S_B & S_1_A) | S_2_A
#     return S_B, S_R, S_A, S_1_A, S_2_A
#
#
# def remove_singletons_in_A_wrt_set_B(genome_A, gene_set_B):
#     new_genome = Genome(genome_A.get_name())
#     for chromosome in genome_A:
#         if chromosome.is_circular():
#             if len(chromosome.get_gene_set().intersection(gene_set_B)):
#                 new_genome.append(chromosome)
#         else:
#             new_genome.append(chromosome)
#     return new_genome
#
#
# def remove_singletons_in_ord_A_wrt_2_dupl_B(ord_genome, all_dupl_genome):
#     new_genome = Genome(ord_genome.get_name())
#     all_dupl_gene_set = all_dupl_genome.get_gene_set()
#     matching, _ = all_dupl_genome.convert_to_matching()
#     matching = set(matching)
#     s_1_all_dupl, s_2_all_dupl = all_dupl_genome.get_partitioned_gene_set()
#
#     for chromosome in ord_genome:
#         if chromosome.is_circular():
#             s_1_chr, s_2_chr = chromosome.get_partitioned_gene_set()
#
#             if not len(s_2_chr) and s_1_chr.intersection(s_1_all_dupl) == s_1_chr:
#                 chromosome_matching, _ = chromosome.convert_to_matching()
#                 flag = False
#                 for u, v in chromosome_matching:
#                     if (u, v) not in matching and (v, u) not in matching:
#                         flag = True
#
#                 if flag:
#                     new_genome.append(chromosome)
#             elif len(chromosome.get_gene_set().intersection(all_dupl_gene_set)):
#                 new_genome.append(chromosome)
#         else:
#             new_genome.append(chromosome)
#
#     return new_genome
