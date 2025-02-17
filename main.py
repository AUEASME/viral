import json
import tensorflow as tf
from easme.parse import parse_proteins, determine_mutability
from easme.selection import (
    fitness_proportionate_selection,
    truncation,
)

# Import protein data.
PROTEINS = parse_proteins("dat/json/sequences.json", 774)
MUTABILITY = determine_mutability(PROTEINS)

# Set evolutionary algorithm parameters.
MAX_POPULATION_SIZE = 100
MU_NUM_PARENTS = 2
LAMBDA_NUM_CHILDREN = 10
NUM_GENERATIONS = 100000

# Load spam filter models.
VALIDITY_MODEL = tf.keras.models.load_model("dat/models/validity.keras")
AGGREGATION_MODEL = tf.keras.models.load_model("dat/models/agg.keras")


class Position:
    def __init__(self, position):
        self.position = position
        self.amino_acids = []


def find_conserved_subsequences(proteins):
    unique_amino_acids_at_position = []
    for i in range(len(proteins[0].amino_acid_sequence)):
        position = Position(i)
        for protein in proteins:
            if protein.amino_acid_sequence[i] not in position.amino_acids:
                position.amino_acids.append(protein.amino_acid_sequence[i])
        unique_amino_acids_at_position.append(position)

    # Remove all cases where len(position.amino_acids) > 1.
    unique_amino_acids_at_position = [
        position
        for position in unique_amino_acids_at_position
        if len(position.amino_acids) == 1
    ]

    return unique_amino_acids_at_position


CONSERVED_SUBSEQUENCES = find_conserved_subsequences(PROTEINS)


def test_for_conserved_subsequences(protein):
    # Check if the amino acid sequence of this protein violates any conserved subsequences.
    for position in CONSERVED_SUBSEQUENCES:
        if protein.amino_acid_sequence[position.position] not in position.amino_acids:
            protein.fitness = 0.0
            return False

    return True


def calculate_fitness_from_spam_filter(protein):
    if not test_for_conserved_subsequences(protein):
        return 0.0

    validity_score = VALIDITY_MODEL.predict(protein.hydro_list)[0][0]
    aggregation_score = AGGREGATION_MODEL.predict(protein.hydro_list)[0][0]
    return validity_score * aggregation_score


def main():
    # Initialize population.
    population = []
    for protein in PROTEINS:
        if protein.date_of_discovery:
            population.append(protein)

    # Remove duplicates.
    for protein in population:
        for other_protein in population:
            if other_protein.dna_sequence_coding == protein.dna_sequence_coding:
                population.remove(other_protein)

    # Calculate fitness.
    for protein in population:
        protein.fitness = calculate_fitness_from_spam_filter(protein)

    # Run evolutionary algorithm.
    max_fitness = 0.0
    for i in range(NUM_GENERATIONS):
        print(f"== GENERATION {i+1} ==")

        # Generate offspring.
        for _ in range(LAMBDA_NUM_CHILDREN):
            # Choose parents.
            parents = fitness_proportionate_selection(population, MU_NUM_PARENTS)
            # Crossover.
            child = parents[0].reproduce(parents[1])
            # Mutate.
            child.mutate()
            # Test fitness.
            child.fitness = calculate_fitness_from_spam_filter(child)
            # Add child to population.
            population.append(child)

        # Select survivors.
        population = truncation(population, MAX_POPULATION_SIZE)
        max_fitness = max([protein.fitness for protein in population])
        print(f"Max fitness: {max_fitness}")
        output = f"dat/out2/gen_{i+1}.json"
        json_output = [protein.dump_to_json() for protein in population]
        with open(output, "w+") as file:
            json.dump(json_output, file, indent=2)

    # Save population to JSON.
    for protein in population:
        protein.save_to_json(f"dat/out2/{protein.name}.json")


if __name__ == "__main__":
    main()
