import random
from easme.parse import parse_proteins, determine_mutability
from easme.fitness import calculate_fitness_from_consensus
from easme.selection import (
    get_living_population,
    prune,
    fitness_proportionate_selection,
    truncation,
)

# Import protein data.
PROTEINS = parse_proteins("data/json/sequences.json", 774)
MUTABILITY = determine_mutability(PROTEINS)

# Set evolutionary algorithm parameters.
MAX_POPULATION_SIZE = 100
MU_NUM_PARENTS = 10
LAMBDA_NUM_CHILDREN = 10
NUM_GENERATIONS = 100


def main():
    # Initialize population.
    population = []
    for protein in PROTEINS:
        if protein.location_of_origin == "Alabama":
            population.append(protein)

    # Remove duplicates.
    for protein in population:
        for other_protein in population:
            if other_protein.dna_sequence_coding == protein.dna_sequence_coding:
                population.remove(other_protein)

    # Calculate fitness.
    for protein in population:
        protein.fitness = calculate_fitness_from_consensus(protein, MUTABILITY)

    # Run evolutionary algorithm.
    for i in range(NUM_GENERATIONS):
        print(f"== COMMENCING GENERATION {i+1} ==")

        # Get living population as a flat list.
        living_population = get_living_population(population)

        # Select parents.
        parents = fitness_proportionate_selection(living_population, MU_NUM_PARENTS)

        # Generate offspring.
        for j in range(LAMBDA_NUM_CHILDREN):
            # Select two random (unique) parents.
            parents = random.sample(parents, 2)
            child = parents[0].reproduce(parents[1], MUTABILITY)
            if child:
                print(
                    f"Generated child GEN{i+1}SER{child.serial_number} with fitness {child.fitness}."
                )
                child.mutate(MUTABILITY)

        # Select survivors.
        living_population_with_new_children = get_living_population(population)
        survivors = truncation(living_population_with_new_children, MAX_POPULATION_SIZE)
        population = prune(population, survivors)

    # Save population to JSON.
    for protein in population:
        protein.save_to_json(f"data/out/{protein.name}.json")


if __name__ == "__main__":
    main()
