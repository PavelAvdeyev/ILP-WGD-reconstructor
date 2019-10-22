

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
