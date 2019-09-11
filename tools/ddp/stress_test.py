from double_dist import solve_double_distance_problem as ddp_ilp
from brute_force_double_dist import solve_ddp_brute as ddp_brute

import argparse
import os
from termcolor import colored
from datetime import datetime


class TestConfig():
    datasets = 10
    # Number of rearrangement events after WGD.
    E = [5, 8]
    IND = [0.10, 0.20, 0.50]
    # Percentage of deletions after WGD all events.
    deletion_p2 = 0.8
    # Percentage of insertions after WGD from all events.
    insertion_p2 = 0.2
    num_genes = [8, 10, 16]
    num_chr = [2, 3, 4, 5, 6]
    tl = 4 * 60 * 60  # in seconds


def create(all=0, config=TestConfig(), sim_path=None):
    for num_genes in config.num_genes:
        for num_chr in config.num_chr:
            if num_genes < num_chr: continue
            for e in config.E:
                if e > num_genes: continue
                for ind in config.IND:
                    for i in range(config.datasets):
                        folder_name = 'sample_{}_{}_{}_{}_{}'.format(num_genes, num_chr, e, int(100 * ind), i)
                        os.system(
                            'python2 {} -wgd="DD" -n={} -c={} -ev2={} -ip2={} -dp2={} -o={} -dl2=2 -il2=2'.format(
                                sim_path,
                                num_genes, num_chr, e, config.insertion_p2 * ind,
                                config.deletion_p2 * ind, folder_name))
                        R_file = os.path.join(folder_name, 'R.txt')
                        A_file = os.path.join(folder_name, 'A.txt')
                        ilp_result_file = os.path.join(folder_name, 'ilp_result.txt')
                        ilp_imp_result_file = os.path.join(folder_name, 'ilp_imp_result.txt')
                        start_time = datetime.now()
                        brute_ans = ddp_brute(R_file, A_file, tl=config.tl)
                        after_brute_time = datetime.now()
                        ilp_ans = ddp_ilp(R_file, A_file, ilp_result_file, time_limit=config.tl, ilp_type='basic',
                                          gurobi_log_file=os.path.join(folder_name, "log_basic"))
                        after_ilp_time = datetime.now()
                        ilp_imp_ans = ddp_ilp(R_file, A_file, ilp_imp_result_file, time_limit=config.tl,
                                              ilp_type='improved',
                                              gurobi_log_file=os.path.join(folder_name, "log_improved"))
                        end_time = datetime.now()
                        print(folder_name)
                        print(colored(
                            'Events: {},\n'
                            'ILP basic answer: {},\n'
                            'ILP improved answer: {},\n'
                            'Brute force answer: {}'.format(e, ilp_ans, ilp_imp_ans, brute_ans),
                            'green'))
                        print(colored(
                            'ILP basic time: {},\n'
                            'ILP improved time: {},\n'
                            'Brute force time: {}'.format(after_ilp_time - after_brute_time,
                                                          end_time - after_ilp_time,
                                                          after_brute_time - start_time),
                            'green'))
                        assert ilp_ans == brute_ans == ilp_imp_ans, colored('OOUPS!', 'red') + ' ' + folder_name
                        if not all:
                            return
                        os.system('rm -r ' + folder_name)


def main():
    parser = argparse.ArgumentParser(
        description="Testing tool for basic, improved and brute force solvers for the Double Distance Problem (DDP)",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-s", "--simulator",
                        default=None, metavar="PATH", help="Path to simulator.py")

    args = parser.parse_args()

    create(1, sim_path=args.simulator)


if __name__ == "__main__":
    main()
