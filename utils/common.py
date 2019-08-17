import logging

logger = logging.getLogger()

def version():
    return "1.0.0"


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
