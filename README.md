# ILP-WGD-Reconstructor
#### A Unified ILP Framework for Core Ancestral Genome Reconstruction Problems

One of the key computational problems in comparative genomics is the reconstruction of genomes of ancestral species based on genomes of extant species. Since most dramatic changes in genomic architectures are caused by genome rearrangements, this problem is often posed as minimization of the number of genome rearrangements between extant and ancestral 
genomes. The basic case of three given genomes is known as the *genome median problem*. **Whole genome duplications** (WGDs) represent yet another type of dramatic evolutionary events and inspire 
the reconstruction of pre-duplicated ancestral genomes, referred to as the *genome halving problem*. 

**ILP-WGD-Reconstructor** contains polynomial-size integer linear programming (ILP) formulations for the aforementioned problems. Moreover, it can also produce ILP formulations for the restricted and conserved versions of the median and halving problems, which combine
the advantages of homology- and rearrangements-based methods. The full list of ILP formulations is the following:

- Genome Median Problem (GMP)
- Intermediate Genome Median Problem (IGMP)
- Conserved Genome Median Problem (CGMP)
- Double Distance Problem (DDP) 
- Guided Genome Halving Problem (GGHP)
- Restricted Guided Genome Halving Problem (RGGHP)
- Conserved Guided Genome Halving Problem (CGGHP) 

### Dependencies
**ILP-WGD-Reconstructor** requires third-party tool installations:

- The [Gurobi Optimizer software](https://www.gurobi.com).<br />The instruction how to install the Gurobi Optimizer can be found [here](https://www.gurobi.com/documentation/quickstart.html).
The Gurobi Optimizer has a [free academic license](https://www.gurobi.com/academia/academic-program-and-licenses/). 
   
- The [Gurobi Python Interface](https://www.gurobi.com/documentation/8.1/quickstart_mac/the_gurobi_python_interfac.html) 

- [Networkx Python library](http://networkx.github.io/)

- [Breakpoint Graph library](https://bg.readthedocs.io/)

### Usage
Consult to usage guide for each ILP formulation:
- [Genome Median Problem](docs/gmp.md)
- [Intermediate Genome Median Problem](docs/igmp.md)
- [Conserved Genome Median Problem](docs/cgmp.md)
- [Double Distance Problem](docs/ddp.md)
- [Guided Genome Halving Problem](docs/gghp.md)
- [Restricted Guided Genome Halving Problem](docs/rgghp.md)
- [Conserved Guided Genome Halving Problem](docs/cgghp.md)

### Publications
P. Avdeyev, N. Alexeev, Y. Rong, and M. A. Alekseyev. "A Unified ILP Framework for Core Ancestral Genome Reconstruction 
Problems". (submitted)

P. Avdeyev, N. Alexeev, Y. Rong, and M. A. Alekseyev. "A Unified ILP Framework for Genome Median, Halving, and 
Aliquoting Problems under DCJ". Lecture Notes in Computer Science, 2017. 

### Authors
- Pavel Avdeyev 
- Nikita Alexeev
- Yongwu Rong 
- Max A. Alekseyev

### Contacts
Please report any issues directly to the github issue tracker. 
Also, you can send your feedback to avdeyev@gwu.edu