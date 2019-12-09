# -*- coding: utf-8 -*-


class ILPAnswer(object):
    """ Class storing the information of resulted solution from ILP solver.

    Class stores information obtained from ILP solver (for example, reconstructed genome,
    ILP solver exit status) and provides functionality to write it into files.


    """

    def __init__(self, ov=0, score=0, es=0, genome=None, tm=0):
        """ Initialization of a :class:'ILPAnswer' object.

        Args:
            ov (int): An objective value of ILP instance.
            score (int): A score of considered ancestral reconstruction problem. It is usually obtained from
            objective value.
            es (int): Exit status of ILP solver.
            genome (Genome): Reconstructed ancestral genome by ILP solver.
            tm (int): Time period that have been used for solving an instance.

        Returns:
            no value
        """
        self.obj_val = ov
        self.score = score
        self.exit_status = es
        self.genome = genome
        self.solution_time = tm

    def update_score_by_singletons(self, value : int) -> None:
        """ Update a final score by adding the number of singletons.

        Args:
            value (int): The number of singletons in the input genomes.

        Returns:
            no value
        """
        self.score += value

    def write_stats_file(self, result_out_file : str) -> None:
        """ Write statistics file about performed run

        Args:
            result_out_file (str): The path to file, which stores supplementary information.

        Returns:
            no value
        """
        with open(result_out_file, 'w') as out:
            out.write("# Objective value\n")
            out.write(str(self.obj_val) + "\n")
            out.write("# Total DCJ-Indel Distance\n")
            out.write(str(self.score) + "\n")
            out.write("# Is it solved \n")
            out.write(str(self.exit_status) + "\n")

    def write_genome_file(self, out_genome_file : str) -> None:
        """ Write file that contains gene order of reconstructed ancestral genome in  GRIMM format.

        Args:
            out_genome_file (str): The path to file, where information is recorded.

        Returns:
            no value
        """
        with open(out_genome_file, 'w') as out:
            out.write(">answer\n")
            for chr_ind, chromosome in self.genome.items():
                for chr_type, blocks in chromosome:
                    string = " ".join(value if sign == "+" else (sign + value) for sign, value in blocks)
                    string += " {chr_type}".format(chr_type=chr_type)
                    out.write(string + "\n")
