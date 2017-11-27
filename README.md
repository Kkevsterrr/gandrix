# gandrix

Gandrix is a Genetic Algorithm for de novo discovery of exclusively mutated driver pathways. It searches for sets of genes in a given dataset for genes with a high coverage and exclusivity. 

Gandrix is designed to run on Python 2.7. 

## Usage

Gandrix is a highly flexible platform, that allows for a high level of control over the genetic algorithm's hyper parameters. 

```
$ python gandrix.py --help   
usage: gandrix.py [-h] [--k K] [--eval-only EVAL_ONLY]
                 [--generations GENERATIONS] [--population POPULATION]
                 [--fitness {wext,dendrix}] [--no-size-compute]
                 [--exclude EXCLUDE] [--data {gbm,test,brca}]

Genetic algorithm for de novo discovery of exclusively mutated pathways.

optional arguments:
  -h, --help            show this help message and exit
  --k K                 solution size to search for
  --eval-only EVAL_ONLY
                        only evaluate fitness for given comma-separated genes
  --generations GENERATIONS
                        number of generations to run for.
  --population POPULATION
                        size of population.
  --fitness {wext,dendrix}
                        fitness function to use.
  --no-size-compute     flag whether to disable population size computation.
  --exclude EXCLUDE     exclude certain genes
  --data {gbm,test,brca}
                        dataset to run on
```

To stop gandrix, Ctrl-C will stop execution and print the current hall of fame individuals. 

## How it works

Gandrix uses a genetic algorithm built on DEAP to evolve a "good" set of genes. The goal is to discover a set of genes that has high coverage (it "covers"/is mutated in a large number of cancer patients) with high exclusivity (there is low overlap between patients). 

Each individual in the population is a collection of genes of size *k*, and begins completely randomly initialized ($k$ random genes are uniformly sampled with replacement to create an individual *i*). When the algorith starts, an initial population  *P* (specified by the arguments) is spawned by creating |*P*| individuals. At each generation, every individual is mutated and undergoes crossover. These offspring are evaulated using a specified fitness function, and the best individuals survive. As generations continue, the individuals gradually converge to the best solution. At the conclusion of the algorithm, the "hall of fame" is printed to the user of the best individuals encountered through the life cycle of the algorithm. 

### Fitness functions
In each generation, every individual is evaluated using a fitness function. This fitness function is the codification of  which gene sets are "good" or "bad". Two fitness functions are implemented in Gandrix ("dendrix" or "wext"). The "dendrix" fitness function is a simple function: (set coverage) - (*C* * set overlap), where *C* is a constant (currently heuristically set as 3). The "wext" fitness function uses the approximation statistical test defined in the wext package, and is much slower to compute. This fitness function also rewards sets where each gene individually has exclusive high coverage (as opposed to one single high coverage gene and multiple small exclusive but low coverage genes). 

### Mutation & Crossover
Mutation and crossover are independently applied to every individual with 40% probability (so it is possible that an individual undergo both mutation and crossover). The order of these operations is left to DEAP's default. If an individual is chosen to be mutated, each gene in that individual is replaced with a different randomly sampled gene with 25% probability. If an individual is chosen for crossover, it is mated with another individual uniformly by swapping genes (each gene is swapped with 50% probability).  

### Hall of Fame
After each generation, every individual is evaluated against the best 10 individuals ever seen before. 

## Examples

By default, gandrix uses the BRCA dataset, an initial population of 5,000 individuals (of size 3), through 100 generations. 
```
$ python gandrix.py    

Loaded BRCA dataset.
Collected 202 genes.
Collected 507 patients in mutation matrix.
Searching for good gene sets of size 3 with initial population of 5000 in 100 generations.
Generation: 0 | Unique Inviduals: 2795 | Best Fitness (47,) (56/507: 11.05%) | Best ADK-MYST4(A) ENSG00000247772 3p25.1(A)
Generation: 1 | Unique Inviduals: 2442 | Best Fitness (52,) (55/507: 10.85%) | Best CPAMD8 GTSE1 RB1(D)
^C
Best individuals:
Fitness (229,) (247/507: 48.72%): GATA3 TP53 CTCF
Fitness (220,) (229/507: 45.17%): CADPS MAP3K1 TP53
Fitness (219,) (228/507: 44.97%): CHD1L MAP3K1 TP53
Fitness (217,) (229/507: 45.17%): TP53 MAP3K1 CBFB
Fitness (215,) (224/507: 44.18%): CDH1 ATN1 TP53
Fitness (214,) (229/507: 45.17%): MYOM3 PIK3CA PTEN(D)
Fitness (210,) (228/507: 44.97%): TP53 ENSG00000245029 MAP3K1
Fitness (207,) (225/507: 44.38%): FAM179A MAP3K1 TP53
Fitness (202,) (259/507: 51.08%): PIK3CA PTEN(D) MIR21(A)
Fitness (199,) (199/507: 39.25%): NINL DALRD3 TP53
```

The GBM dataset is also bundled with gandrix, as well as a small toy dataset for demonstration purposes. 
```
$ python gandrix.py --data test    
Loaded TEST dataset.
Collected 9 genes.
Collected 10 patients in mutation matrix.
     1 2 3 4 5 6 7 8 9
P1:  # . . . . # # . .
P2:  . # . . . # # . .
P3:  . # . . . # . # .
P4:  . # . . # . . # .
P5:  . . # . # . . # .
P6:  . . # . # . . # .
P7:  . . # . # . . . #
P8:  . . # . # . . . #
P9:  . . # # . . . . #
P10: . . # # . . . . #
Searching for good gene sets of size 3 with initial population of 10000 in 100 generations.
Generation: 0 | Unique Inviduals: 159 | Best Fitness (-3,) (9/10: 90.00%) | Best gene5 gene7 gene3
Generation: 1 | Unique Inviduals: 142 | Best Fitness (8,) (8/10: 80.00%) | Best gene9 gene1 gene2
Generation: 2 | Unique Inviduals: 137 | Best Fitness (9,) (9/10: 90.00%) | Best gene1 gene8 gene9
^C
Best individuals:
Fitness (10,) (10/10: 100.00%): gene9 gene7 gene8
Fitness (10,) (10/10: 100.00%): gene4 gene5 gene6
Fitness (10,) (10/10: 100.00%): gene3 gene2 gene1
Fitness (9,) (9/10: 90.00%): gene5 gene7 gene4
Fitness (9,) (9/10: 90.00%): gene9 gene1 gene8
Fitness (8,) (8/10: 80.00%): gene8 gene4 gene7
Fitness (8,) (8/10: 80.00%): gene4 gene5 gene1
Fitness (8,) (8/10: 80.00%): gene2 gene1 gene9
Fitness (7,) (10/10: 100.00%): gene9 gene6 gene8
Fitness (7,) (10/10: 100.00%): gene2 gene3 gene7
```

### Adding new datasets
