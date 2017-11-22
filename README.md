# gandrix

Gandrix is a Genetic Algorithm for de novo discovery of exclusively mutated driver pathways. It encorporates 

## Usage

Gandrix is a highly flexible platform, that allows for a high level of control over the genetic algorithm's hyper parameters. 

```
$ python solver.py --help                                                                                                                             [18:27:00]
usage: solver.py [-h] [--k K] [--eval-only EVAL_ONLY]
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

### Adding new datasets
