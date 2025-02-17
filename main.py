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
NUM_GENERATIONS = 100000


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


def create_unique_child(generation, parents, mutability):
    # Select two random (unique) parents.
    parents = random.sample(parents, 2)

    # Get all unique amino acid sequences in descendants of both parents.
    # (BEFORE generating child.)
    unique_amino_acid_sequences = set()
    for parent in [parents[0], parents[1]]:
        for descendant in parent.descendants:
            unique_amino_acid_sequences.add(descendant.amino_acid_sequence)

    # Generate child.
    child = parents[0].reproduce(parents[1], mutability)
    if child:
        print(
            f"Generated child GEN{generation+1}SER{child.serial_number} with fitness {child.fitness}."
        )
        # Generate amino acid sequence of new child.
        child.transcribe()
        child.translate()
        test_for_conserved_subsequences(child)
        

        # Check if child has a unique amino acid sequence.
        if (
            (parents[0].amino_acid_sequence != child.amino_acid_sequence)
            and (parents[1].amino_acid_sequence != child.amino_acid_sequence)
            and (child.amino_acid_sequence not in unique_amino_acid_sequences)
        ):
            print(f"Child has a unique amino acid sequence. Keeping…")
            # Mutate child.
            child.mutate(mutability)
        else:
            # If child has a duplicate amino acid sequence, discard it.
            print(f"Child has a duplicate amino acid sequence. Discarding…")
            parents[0].descendants.remove(child)
            parents[1].descendants.remove(child)


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
        print(f"== GENERATION {i+1} ==")

        # Get living population as a flat list.
        living_population = get_living_population(population)

        # Select parents.
        parents = fitness_proportionate_selection(living_population, MU_NUM_PARENTS)

        # Generate offspring.
        for _ in range(LAMBDA_NUM_CHILDREN):
            create_unique_child(i, parents, MUTABILITY)

        # Select survivors.
        living_population_with_new_children = get_living_population(population)
        survivors = truncation(living_population_with_new_children, MAX_POPULATION_SIZE)
        population = prune(population, survivors)

    # Save population to JSON.
    for protein in population:
        protein.save_to_json(f"data/out/{protein.name}.json")


if __name__ == "__main__":
    main()
