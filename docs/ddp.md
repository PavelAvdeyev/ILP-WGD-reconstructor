# Double Distance Problem
`double_dist.py` creates the ILP formulation for the Double Distance Problem
and runs Gurobi Solver.  

```
usage: double_dist.py [-h] [-r PATH] [-t PATH] -o PATH [-tl NUMBER] [-v]

ILP solver for the Double Distance Problem (DDP)

optional arguments:
  -h, --help            show this help message and exit
  -r PATH, --ordinary PATH
                        Ordinary genome in GRIMM format
  -t PATH, --two_dupl PATH
                        2-duplicated genome in GRIMM format
  -o PATH, --out-dir PATH
                        Output directory
  -tl NUMBER, --time_limit NUMBER
                        Time limit for Gurobi solver
  -v, --version         show program's version number and exit

Gurobi solver is required! Input genomes must be in GRIMM format.
```

### Examples
You can try `double_dist.py` on the provided ready-to-use examples:

```bash
 python double_dist.py --ordinary examples/ord.gen --two_dupl examples/a.gen -o test/
```

### Input
`double_dist.py` needs as input:

- Ordinary genome (i.e., each gene/synteny block is present in a single copy) **[[in GRIMM format](http://grimm.ucsd.edu/GRIMM/grimm_instr.html)]** 
- 2-duplicated genome (i.e., each gene/synteny block is present in two copies) **[[in GRIMM format](http://grimm.ucsd.edu/GRIMM/grimm_instr.html)]**
- Path to output directory

Optionally, a time limit may be set for Gurobi Solver (it is 2 hours by default). 

### Output
After running `double_dist.py`, the output directory will contain:
- `result.txt` - file with exit status (0 - stopped by reaching TL, 1 - optimal solution, 2 - other),
optimized value of the objective function, and the distance between the
pre-duplicated genome and 2-duplicated genome
- several log files