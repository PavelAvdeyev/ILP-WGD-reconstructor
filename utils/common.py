# -*- coding: utf-8 -*-

import logging

import os

from utils.genome import Genome

logger = logging.getLogger()


def version() -> str:
    """Returns current version of the tool"""
    return "1.1.0b"


def epilog() -> str:
    return "Gurobi solver is required! Input genomes must be in GRIMM format."


def enable_logging(log_file: str, overwrite: bool) -> None:
    """
    Turns on logging, sets debug levels and assigns a log file.

    Args:
            log_file (str): The path to log file
            overwrite (bool): The indicator for overwriting existing log file

    Returns:
         no value
    """
    console_formatter = logging.Formatter("[%(asctime)s] %(levelname)s: "
                                          "%(message)s", "%Y-%m-%d %H:%M:%S")
    console_log = logging.StreamHandler()
    console_log.setFormatter(console_formatter)
    console_log.setLevel(logging.INFO)
    logger.addHandler(console_log)

    log_formatter = logging.Formatter("[%(asctime)s] %(name)s: %(levelname)s: "
                                      "%(message)s", "%Y-%m-%d %H:%M:%S")
    if overwrite:
        open(log_file, "w").close()
    file_handler = logging.FileHandler(log_file, mode="a")
    file_handler.setFormatter(log_formatter)

    logger.setLevel(logging.INFO)
    logger.addHandler(file_handler)


def get_immediate_subdirectories(a_dir):
    """
    This function gets subdirectories of directory
    """
    return ((os.path.join(a_dir, name), name) for name in os.listdir(a_dir) if os.path.isdir(os.path.join(a_dir, name)))


def get_immediate_files(a_dir):
    """
    This function gets files in the directory
    """
    return ((os.path.join(a_dir, name), name) for name in os.listdir(a_dir) if
            os.path.isfile(os.path.join(a_dir, name)))


def remove_singletons_wrt_gene_set(genome, gene_set):
    new_genome = Genome(genome.get_name())
    number_of_singletons = 0
    for chromosome in genome:
        if chromosome.is_circular() and len(chromosome.get_gene_set().intersection(gene_set)) == 0:
            number_of_singletons += 1
            continue
        new_genome.append(chromosome)
    return new_genome, number_of_singletons


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
