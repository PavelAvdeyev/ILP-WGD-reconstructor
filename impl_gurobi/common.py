from utils.genome import get_extremity_by_gene
import bg.breakpoint_graph


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



