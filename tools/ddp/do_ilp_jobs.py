import argparse
import logging
import os

logger = logging.getLogger()

eps = 10


def make_script(cur_path, tl, tp):
    f = open(os.path.join(cur_path, 'script_'+tp), 'w')
    f.write('#!/bin/bash\n')
    params = "#SBATCH -o {0}/dd%j.out\n#SBATCH -e {0}/dd%j.err\n\n".format(cur_path)
    f.write(params)
    command = 'python3 do_ilp_by_folder_name.py -fn="{}" -tl={} -tp="{}"'.format(cur_path, tl * 60, tp)
    f.write(command)
    f.close()
    return f


def submit(f, tl):
    os.system('chmod +x {}'.format(f.name))
    # print('sbatch -t {} {}'.format(tl*60 + eps, f.name))
    os.system('sbatch -t {} {}'.format(tl + eps, f.name))


def do(all=0, path=None, tl=10, tp='improved'):
    for cur_path, dirs, files in os.walk(path):
        if not 'A.txt' in files:
            continue
        if 'result.txt' in files:
            continue
        logging.info("Working with directory {0}".format(cur_path))
        if tp != 'both':
            f = make_script(cur_path, tl, tp)
            submit(f, tl)
        else:
            f = make_script(cur_path, tl, "improved")
            submit(f, tl)
            f = make_script(cur_path, tl, "basic")
            submit(f, tl)

        if not all:
            return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Submit jobs for ddp")
    parser.add_argument("-fn", "--folder_name", type=str, default=None,
                        help="Folder with tests.")
    parser.add_argument("-tl", "--time_limit", type=int, default=10,
                        help="Time limit for gurobi in minutes.")
    parser.add_argument("-tp", "--type", type=str, default='improved',
                        help="ILP type: improved, basic or both")

    param = parser.parse_args()

    path = param.folder_name
    tl = param.time_limit
    tp = param.type
    do(1, path, tl, tp)
