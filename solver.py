# pylint: disable=superfluous-parens,no-member
import argparse
import random
import sys
import warnings
import math
import os

from saddlepoint import saddlepoint
from deap import creator, base, tools, algorithms

# Compare output to http://compbio-research.cs.brown.edu/projects/multi-dendrix/results/Multi/BRCA/
def get_args():
    """
    Sets up argparse and collects arguements.
    """

    parser = argparse.ArgumentParser(description='Genetic algorithm for BRCA pathways.')
    parser.add_argument('--k', type=int, action='store', default=3, help='solution size to search for')
    parser.add_argument('--eval-only', action='store', default=None, help='only evaluate fitness for given comma-separated genes')
    parser.add_argument('--generations', type=int, action='store', default=100, help="number of generations to run for.")
    parser.add_argument('--population', type=int, action='store', default=10000, help="size of population.")
    parser.add_argument('--fitness', choices=('wext', 'dendrix'), action='store', default='dendrix', help="fitness function to use.")
    parser.add_argument('--no-size-compute', action='store_true', default=False, help="flag whether to disable population size computation.")
    parser.add_argument('--exclude', action='store', default=None, help="exclude certain genes")
    parser.add_argument('--data', action='store', choices=('gbm', 'test', 'brca'), default="brca", help="dataset to run on")
    args = parser.parse_args()
    return args


def observed_values(individual, num, gene_to_cases):
    """
    Critical function for computing metrics & table for row-exclusivity test.
    """
    # Construct the table
    k = len(individual)
    mutated_patients = set(patient for gene in individual for patient in gene_to_cases[gene])
    tbl = [0] * 2**k
    exclusive, overlap = 0, 0
    for patient in mutated_patients:
        # Determine the binary index representing which genes
        # have mutations in this sample
        i = 0
        genes_mutated = 0
        for idx, gene in enumerate(individual):
            if patient in gene_to_cases[gene]:
                genes_mutated += 1
                i = i | (1 << idx)
        tbl[i] += 1

        # Increment the Z/T count if this mutation is exclusive or
        # co-occurring
        overlap += int(genes_mutated > 1)
        exclusive += int(genes_mutated == 1)

    # Finish the table and compute X
    tbl[0] = num - len(mutated_patients)
    coverage = [len(gene_to_cases[gene]) for gene in individual]

    return coverage, exclusive, overlap, tbl


def re_test(exclusive, coverage, tbl):
    """
    Computes row-exclusivity statistical test.
    """
    total = sum(tbl)
    norm = [[float(x_i)/total] * total for x_i in coverage]
    # Ignore warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        p_value = saddlepoint(exclusive, coverage, norm)

    return p_value


def get_indices(num):
    """
    Method used by deap to create each aspect of an individual. 
    Simply returns a random int between 0 and given num - 1.
    """

    return random.randint(0, num - 1)


def fitness_from_genes(full_gene_list, patient_matrix, fitness_function, test_genes):
    """
    Given an array of gene names, compute the fitness of this population.
    Example: ['TP53', 'GATA3']
    """
    individual = []
    # Fitness evaluation is done by index for speed, so first get each index
    for gene in test_genes:
        try:
            idx = full_gene_list.index(gene)
            individual.append(idx)
        except ValueError:
            print("ERROR: Could not find %s" % gene)
            continue
    fitness = fitness_function(full_gene_list, patient_matrix, individual)
    coverage = get_coverage_exclusivity(patient_matrix, individual)[0]
    print("Fitness %s (%d/%d: %0.2f%%)" % (str(fitness), coverage, len(patient_matrix), round(100.0*float(coverage)/float(len(patient_matrix)), 2)))


def statistical_fitness(patient_matrix, individual):
    """
    Compute fitness using wext's row-exclusivity test (re_test). 

    Expected input data for re_test does not fully match up with how data is stored, 
    so massage data into better format.
    """

    num = len(patient_matrix)
    gene_to_cases = {}
    for gene in individual:
        for idx in range(0, num):
            if patient_matrix[idx][gene]:
                if gene not in gene_to_cases:
                    gene_to_cases[gene] = set()
                gene_to_cases[gene].add(idx)
    coverage, exclusive, overlap, tbl = observed_values(individual, num, gene_to_cases)
    if overlap >= exclusive:
        p_val = 1.0
    else:
        p_val = re_test(exclusive, coverage, tbl)
        if p_val == 0:
            p_val = 0.0000000000000001
    
    return -1 * math.log(p_val, 10),


def dendrix_fitness(patient_matrix, individual):
    """
    Use dendrix-like fitness measure of rewarding coverage and punishing overlap.
    """
    coverage, overlap = get_coverage_exclusivity(patient_matrix, individual)
    # Comma at the end of this return is *required* by deap
    return coverage - 3 * overlap,


def get_coverage_exclusivity(patient_matrix, individual):
    """
    Computes simple coverage & exclusivity measures for a given individual for 
    the fitness function to use.
    """
    coverage = {}
    exclusivity = 0
    for idx in individual:
        for j in range(0, len(patient_matrix)):
            if patient_matrix[j][idx]:
                if j in coverage:
                    exclusivity += 1
                coverage[j] = True
    return len(coverage), (exclusivity),


def parse_genes(filename, exclude):
    """
    Parses a given gene_file, and excludes any genes in the excludes list.
    """

    full_list = []
    with open(filename, "r") as genes_file:
        genes = genes_file.readlines()
        for gene in genes:
            gene = gene.strip()
            if gene not in exclude:
                full_list.append(gene)
    print("Collected %d genes." % len(full_list))
    return full_list


def parse_patients(filename, genes, exclude):
    """
    Parses a given patients file, and excludes any genes from the exclude list.
    """

    matrix = []
    with open(filename, "r") as patient_file:
        i = 0
        for line in patient_file.readlines():
            matrix.append([])
            mutated = line.strip().split("\t")
            mutated.pop(0)
            for gene in genes:
                if gene not in exclude:
                    matrix[i].append((gene in mutated))
            i += 1
    print("Collected %d patients in mutation matrix." % len(matrix))
    return matrix


def pretty_print(matrix):
    """
    Pretty print method exclusively for use with the test dataset for easy visualization
    of what's going on. This method does not work well with larger datasets (as the screen
    is generally not big enough to fit the whole matrix without running over a line), so 
    it's strongly recommended to leave it as such for automatic invocation with test data. 
    """

    print(" " * 5),
    for idx in range(1, 10):
        sys.stdout.write("%d " % idx)
    print("")
    for idx, genes in enumerate(matrix):
        print("P%d:%s" % ((idx + 1), (" " * (2 if (idx < 9) else 1)))),
        for gene in genes:
            if gene:
                char = "# "
            else: 
                char = ". "
            sys.stdout.write(char)
        print("")


def get_unique_population_size(population, no_size_compute):
    """
    Computes number of unique invidiuals in a population.
    Since deap does not recognize that [gene1, gene2] is equivalent to [gene2, gene1],
    sorting each individual is required. 
    """

    if no_size_compute:
        return 0
    
    uniques = {}
    for i in population:
        uniques[str(sorted(list(i)))] = True
    return len(uniques)


def add_to_hof(hof, population):
    fits = {}
    for best in hof:
        fits[str(sorted(list(best)))] = (best, best.fitness.values)
    for ind in population:
        fits[str(sorted(list(ind)))] = (ind, ind.fitness.values)
    sorted_fits = sorted(fits.iteritems(), reverse=True, key=lambda (k,v): (v[1],k)) 
    hof = []
    for best in sorted_fits[:10]:
        hof.append(best[1][0])
    return hof


def genetic_solve(options, genes, matrix, fitness_function):
    """
    Perform actual genetic solving with given options on a given gene and patient matrix,
    using given fitness function.
    """

    toolbox = base.Toolbox()
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

    toolbox.register("attr", get_indices, len(genes))
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr, n=options["solution_size"])
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", fitness_function, matrix)
    toolbox.register("mate", tools.cxUniform, indpb=0.5)
    toolbox.register("mutate", tools.mutUniformInt, low=0, up=len(genes)-1, indpb=0.25)
    toolbox.register("select", tools.selTournament, tournsize=2)
    
    hall = [] 
    population = toolbox.population(n=options["population_size"])
    
    try:
        offspring = []
        for gen in range(options["num_generations"]):    
            offspring = algorithms.varAnd(population, toolbox, cxpb=0.4, mutpb=0.4)
            fits = toolbox.map(toolbox.evaluate, offspring)
            for fit, ind in zip(fits, offspring):
                ind.fitness.values = fit
            
            population = toolbox.select(offspring, k=len(population))
            best_one = toolbox.select(offspring, k=1)[0]
            best_fit = fitness_function(matrix, best_one)
            best_coverage = get_coverage_exclusivity(matrix, best_one)[0]
            #hall.update(population)
            hall = add_to_hof(hall, population)
            print("Generation: %d | Unique Inviduals: %d | Best Fitness %s (%d/%d: %0.2f%%) | Best %s %s\r" % (gen, get_unique_population_size(population, options["no_size_compute"]), str(best_fit), best_coverage, len(matrix), round(100.0*float(best_coverage)/float(len(matrix)), 2), get_names(genes, best_one).strip(), " "*20))
    except KeyboardInterrupt:
        print("\nBest individuals:")

    return hall


def get_names(genes, individual):
    """ 
    Build string of space separated gene names from an individual.
    """
    name = ""
    for idx in individual:
        name += genes[idx] + " "
    return name


def print_results(genes, matrix, offspring, fitness_function):
    """
    Print final results from offspring.
    """
    num_p = len(matrix)
    for ind in offspring:
        coverage = float(get_coverage_exclusivity(matrix, ind)[0])
        sys.stdout.write("Fitness %s (%d/%d: %0.2f%%): " % (fitness_function(matrix, ind), coverage, num_p, round(100.0*coverage/float(num_p), 2)))
        for idx in ind:
            sys.stdout.write(genes[idx]+" ")
        print("")


def validate_dataset(dataset):
    """
    Simple error handling to help users correct comment errors given dataset.
    """

    if not os.path.exists(os.path.join("data", dataset)):
        print("ERROR: Can't find folder %s in folder data/. Is that folder present?" % dataset)
        sys.exit()

    if not os.path.exists(os.path.join("data", dataset, "gene_list.txt")):
        print("ERROR: Can't find file \"gene_list.txt\" in data/%s. Is the genes list present?" % dataset)
        sys.exit()
    
    if not os.path.exists(os.path.join("data", dataset, "patients.txt")):
        print("ERROR: Can't find file \"patientstxt\" in data/%s. Is the genes list present?" % dataset)
        sys.exit()
    
    return os.path.join("data", dataset, "gene_list.txt"), os.path.join("data", dataset, "patients.txt")


def driver():
    """
    Main workflow driver for the solver. Parses flags and input data, and initiates solving.
    """

    args = get_args() 
    genes_list_file, patients_file = validate_dataset(args.data) 

    without = args.exclude
    if without:
        print("Excluding %s from analysis" % without)
        exclude = [g.strip() for g in without.split(",")]
    else:
        exclude = []

    if args.fitness == "wext":
        fitness_function = statistical_fitness
    elif args.fitness == "dendrix":
        fitness_function = dendrix_fitness
    else:
        print("ERROR: Unknown fitness function.")
        sys.exit()
    
    eval_only = args.eval_only
    full_gene_list = parse_genes(genes_list_file, exclude)
    patient_matrix = parse_patients(patients_file, full_gene_list, exclude)

    if args.data == "test":
        pretty_print(patient_matrix)

    if eval_only:
        print("Computing fitness for M = %s" % eval_only)
        ind = []
        for gene in [x.strip() for x in eval_only.split(",")]:
            ind.append(full_gene_list.index(gene))
        print_results(full_gene_list, patient_matrix, [ind], fitness_function)
        sys.exit()

    
    options = {}
    options["population_size"] = args.population
    options["solution_size"] = args.k
    options["num_generations"] = args.generations
    options["no_size_compute"] = args.no_size_compute
    hall_of_fame = genetic_solve(options, full_gene_list, patient_matrix, fitness_function) 
    print_results(full_gene_list, patient_matrix, hall_of_fame, fitness_function)


if __name__ == "__main__":
    driver()


