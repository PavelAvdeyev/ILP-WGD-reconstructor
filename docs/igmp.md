# Intermediate Genome Median Problem
`igmp.py` creates the ILP formulation for the Intermediate Genome Median Problem
and runs Gurobi Solver.  

```
usage: igmp.py [-h] [-f PATH PATH PATH] -o PATH [-tl NUMBER] [-v]

ILP solver for the Intermediate Genome Median Problem (IGMP)

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
You can try `igmp.py` on the provided ready-to-use examples:

```bash
 python igmp.py --files examples/S1.gen examples/S2.gen examples/S3.gen -o test/
```

### Input
`igmp.py` needs as input:

- Three files with ordinary genomes **[[in GRIMM format](http://grimm.ucsd.edu/GRIMM/grimm_instr.html)]**.
The third genome is considered as an outgroup genome.
An intermediate genome will be found between the first and second genome.
In an ordinary genome, each gene/synteny block is present in a single copy.  
- Path to output directory

Optionally, a time limit may be set for Gurobi Solver (it is 2 hours by default).

### Output
After running `igmp.py`, the output directory will contain:
- `median.gen` - file with ancestral (intermediate) genome in GRIMM format 
- `result.txt` - file with exit status (0 - stopped by reaching TL, 1 - optimal solution, 2 - other),
optimized value of the objective function, and gmp-score (the total distance from the given genomes)   
- several log files