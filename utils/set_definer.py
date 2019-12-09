# -*- coding: utf-8 -*-

import itertools

from utils.genome import get_extremity_by_gene


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


def define_equiv_function(gene_multiset, vertex2ind, bar_vertex2ind):
    equiv = {}
    for gene, copies in gene_multiset.items():
        for left in [True, False]:
            temp = []
            for i in range(1, copies + 1):
                temp.append(bar_vertex2ind[get_extremity_by_gene(gene, left, i)])
            equiv[vertex2ind[get_extremity_by_gene(gene, left, 1)]] = temp
    return equiv

