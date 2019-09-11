from double_dist import *
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Run solver on test")
    parser.add_argument("-fn", "--folder_name", type=str, default=None,
                        help="Folder with target (R.txt) and duplicated (A.txt) genomes.")
    parser.add_argument("-tl", "--time_limit", type=int, default=5,
                        help="Time limit for gurobi in minutes.")
    parser.add_argument("-tp", "--type", type=str, default='improved',
                        help="ILP type: improved, basic or both")

    param = parser.parse_args()

    path = param.folder_name
    tl = param.time_limit
    tp = param.type
    ordinary_genome = os.path.join(path, "R.txt")
    all_dupl_genome = os.path.join(path, "A.txt")
    if tp != 'both':
        out_file = os.path.join(path, "result_" + tp + ".txt")
        solve_double_distance_problem(ord_genome_file=ordinary_genome, all_dupl_genome_file=all_dupl_genome,
                                      out_result_file=out_file, time_limit=tl,
                                      ilp_type=tp, gurobi_log_file=os.path.join(path, "gurobi_log_{}.log".format(tp)))
    else:
        out_file_imp = os.path.join(path, "result_improved.txt")
        out_file_basic = os.path.join(path, "result_basic.txt")
        solve_double_distance_problem(ord_genome_file=ordinary_genome,
                                      all_dupl_genome_file=all_dupl_genome,
                                      out_result_file=out_file_imp,
                                      time_limit=tl,
                                      ilp_type="improved",
                                      gurobi_log_file=os.path.join(path, "gurobi_log_improved.log"))
        solve_double_distance_problem(ord_genome_file=ordinary_genome,
                                      all_dupl_genome_file=all_dupl_genome,
                                      out_result_file=out_file_basic,
                                      time_limit=tl,
                                      ilp_type="basic", gurobi_log_file=os.path.join(path, "gurobi_log_basic.log"))
