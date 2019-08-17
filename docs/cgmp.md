# Conserved Genome Median Problem
`cgmp.py` creates the ILP formulation for the Conserved Genome Median Problem 
and runs Gurobi Solver.  

```
usage: cgmp.py [-h] [-f PATH PATH PATH] -o PATH [-tl NUMBER] [-v]

ILP solver for the Conserved Genome Median Problem (CGMP)

optional arguments:
  -h, --help            show this help message and exit
  -f PATH PATH PATH, --files PATH PATH PATH
                        Three files in GRIMM format
  -o PATH, --out-dir PATH
                        Output directory
  -tl NUMBER, --time_limit NUMBER
                        Time limit for Gurobi solver
  -v, --version         show program's version number and exit

Gurobi solver is required! Input genomes must be in GRIMM format.
```

### Examples
You can try `cgmp.py` on the provided ready-to-use examples:

```bash
 python cgmp.py --files examples/S1.gen examples/S2.gen examples/S3.gen -o test/
```

### Input
`cgmp.py` needs as input:

- Three files with ordinary genomes **[[in GRIMM format](http://grimm.ucsd.edu/GRIMM/grimm_instr.html)]** <br />
 (in an ordinary genome, each gene/synteny block is present in a single copy) 
- Path to output directory

Optionally, time limit may be set for Gurobi Solver (it is 2 hours by default). 

### Output
After running `cgmp.py`, the output directory will contain:
- `median.gen` - file with ancestral (median) genome in GRIMM format 
- `result.txt` - file with exit status (0 - stopped by reaching TL, 1 - optimal solution, 2 - other), 
objective value of optimized functional, and gmp-score (the total distance from the given genomes)   
- several log files