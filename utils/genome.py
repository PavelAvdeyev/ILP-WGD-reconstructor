# -*- coding: utf-8 -*-
""" The description of genome and chromosome objects

This module contains functionality for reading and storing Genome and Chromosome objects.
"""

from collections import defaultdict


class Genome(object):
    """ Class providing implementation of genome and operations on it.

    Attributes:
        name (str): The name of the genome.
        chromosomes (list[Chromosome]): The list of chromosomes in the genome.
    """

    def __init__(self, name=""):
        """ Initialization of a :class:'Genome' object.

        Args:
            name (str): The name of the genome.
        """

        self.name = name
        self.chromosomes = []

    def set_name(self, name: str) -> None:
        """Set the genome name."""
        self.name = name

    def append(self, chromosome) -> None:
        """ Add a chromosome to the genome."""
        self.chromosomes.append(chromosome)

    def add_chromosome(self, is_circular, chromosome):
        """ Add a chromosome to the genome by constacting emplace."""
        self.chromosomes.append(Chromosome(is_circular, chromosome))

    def number_of_chromosomes(self) -> int:
        """Return the number of chromosomes in the genome."""
        return len(self.chromosomes)

    def get_gene_multiset(self):
        """Get multiset of genes that are present in genome"""
        result = defaultdict(int)
        for chrm in self.chromosomes:
            chrm._get_gene_multiset(result)
        return result

    def convert_to_genome_graph(self):
        appearing_genes = defaultdict(int)
        list_of_edges, list_of_telomers = [], []
        for chrm in self.chromosomes:
            chrm.convert_to_genome_graph(list_of_edges, list_of_telomers, appearing_genes)
        return list_of_edges, list_of_telomers

    def convert_to_contracted_genome_graph(self):
        list_of_edges, list_of_telomers = [], []
        for chrm in self.chromosomes:
            chrm.convert_to_contracted_genome_graph(list_of_edges, list_of_telomers)
        return list_of_edges, list_of_telomers

    def get_chromosome(self, ind):
        """Get a chromosome with specific index"""
        return self.chromosomes[ind]

    def get_name(self):
        return self.name

    def __iter__(self):
        return iter(self.chromosomes)


class Chromosome(object):
    """ Class providing implementation of chromosome and operations on it.

        Attributes:
            is_circ (bool): The indicator of chromosome circularity (True -- circular chromosome).
            blocks (list[T]): The list of genes/synteny blocks in chromosome.
    """

    def __init__(self, is_circular=False, blocks=None):
        if blocks is None:
            blocks = []
        self.is_circ = is_circular
        self.blocks = blocks

    def size(self):
        return len(self.blocks)

    def set_circular(self, is_circular):
        self.is_circ = is_circular

    def append(self, block):
        self.blocks.append(block)

    def is_circular(self):
        return self.is_circ

    def get_gene_set(self):
        result = defaultdict(int)
        self._get_gene_multiset(result)
        return set(result.keys())

    def get_gene_multiset(self):
        result = defaultdict(int)
        self._get_gene_multiset(result)
        return result

    def _get_gene_multiset(self, result):
        for block in self.blocks:
            result[abs(block)] += 1

    def convert_to_genome_graph(self, edges=None, telomers=None, appearing_genes=defaultdict(int)):
        if telomers is None:
            telomers = []

        if edges is None:
            edges = []

        appearing_genes[abs(self.blocks[0])] += 1
        start = last = appearing_genes[abs(self.blocks[0])]

        for block1, block2 in zip(self.blocks, self.blocks[1:]):
            appearing_genes[abs(block2)] += 1
            edges.append((get_extremity_by_gene(block1, False, last),
                          get_extremity_by_gene(block2, True, appearing_genes[abs(block2)])))
            last = appearing_genes[abs(block2)]

        if self.is_circ:
            edges.append((get_extremity_by_gene(self.blocks[-1], False, last),
                          get_extremity_by_gene(self.blocks[0], True, start)))
        else:
            telomers.append(get_extremity_by_gene(self.blocks[-1], False, last))
            telomers.append(get_extremity_by_gene(self.blocks[0], True, start))

    def convert_to_contracted_genome_graph(self, edges=None, telomers=None):
        if telomers is None:
            telomers = []

        if edges is None:
            edges = []

        for block1, block2 in zip(self.blocks, self.blocks[1:]):
            edges.append((get_extremity_by_gene(block1, False, 1), get_extremity_by_gene(block2, True, 1)))

        if self.is_circ:
            edges.append(
                (get_extremity_by_gene(self.blocks[-1], False, 1), get_extremity_by_gene(self.blocks[0], True, 1)))
        else:
            telomers.append(get_extremity_by_gene(self.blocks[-1], False, 1))
            telomers.append(get_extremity_by_gene(self.blocks[0], True, 1))

    def __iter__(self):
        return iter(self.blocks)


def get_extremity_by_gene(block, is_left, gene_copy=1):
    if block < 0:
        vertex = str(abs(block)) + ('h' if is_left else 't') + "_" + str(gene_copy)
    else:
        vertex = str(abs(block)) + ('t' if is_left else 'h') + "_" + str(gene_copy)
    return vertex


def parse_genome_in_grimm_file(full_name):
    """
    This function parses a genome file in grimm format.
    # - comments
    > [name] - name of a genome
    return genome object (see Genome class)
    """
    genome = Genome()

    with open(full_name, 'r') as f:
        chromosome = Chromosome()
        for line in f:
            line = line.strip(' \n\t')

            if line.find("#") != -1 or len(line) == 0:
                continue

            if line.find('>') != -1:
                genome.set_name(line[1:])
                continue

            for gene in line.split(' '):
                str_gene = gene.strip(' \n\t')
                if str_gene == '$':
                    chromosome.set_circular(False)
                    genome.append(chromosome)
                    chromosome = Chromosome()
                elif str_gene == '@':
                    chromosome.set_circular(True)
                    genome.append(chromosome)
                    chromosome = Chromosome()
                elif len(str_gene) != 0:
                    chromosome.append(int(str_gene))
    return genome


def parse_genomes_in_grimm_file(full_name):
    """
    This function parses genomes in file of grimm format.
    # - comments
    > [name] - name on genome
    return [genomes] where genome is genome object (see Genome class)
    """
    genomes = []

    with open(full_name, 'r') as f:
        genome = Genome()
        chromosome = Chromosome()

        for line in f:
            line = line.strip(' \n\t')

            if line.find("#") != -1 or len(line) == 0:
                continue

            if line.find('>') != -1:
                if chromosome.size() != 0:
                    genome.append(chromosome)

                genomes.append(genome)

                chromosome = Chromosome()
                genome = Genome()
                genome.set_name(line[1:])
                continue

            for gene in line.split(' '):
                str_gene = gene.strip(' \n\t')
                if str_gene == '$':
                    chromosome.set_circular(False)
                    genome.append(chromosome)
                    chromosome = Chromosome()
                elif str_gene == '@':
                    chromosome.set_circular(True)
                    genome.append(chromosome)
                    chromosome = Chromosome()
                elif len(str_gene) != 0:
                    chromosome.append(int(str_gene))

        if chromosome.size() != 0:
            genome.append(chromosome)

        genomes.append(genome)

    return genomes[1:]
