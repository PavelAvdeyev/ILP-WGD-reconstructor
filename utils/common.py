import logging

import os

from utils.genome import Genome

logger = logging.getLogger()


def version():
    return "1.1.0b"


def epilog():
    return "Gurobi solver is required! Input genomes must be in GRIMM format."


def enable_logging(log_file, overwrite):
    """
    Turns on logging, sets debug levels and assigns a log file
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
    This function get subdirectories
    """
    return ((os.path.join(a_dir, name), name) for name in os.listdir(a_dir) if os.path.isdir(os.path.join(a_dir, name)))


def get_immediate_files(a_dir):
    """
    This function get files
    """
    return ((os.path.join(a_dir, name), name) for name in os.listdir(a_dir) if
            os.path.isfile(os.path.join(a_dir, name)))


def remove_singletons_wrt_gene_set(genome, gene_set):
    new_genome = Genome(genome.get_name())
    number_of_singletons = 0
    for chromosome in genome:
        if chromosome.is_circular() and len(chromosome.get_gene_set().intersection(gene_set)) == 0:
            number_of_singletons += 1
        new_genome.append(chromosome)
    return new_genome, number_of_singletons
