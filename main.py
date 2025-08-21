import json
from difflib import SequenceMatcher
from copy import deepcopy
import random
from easme.evolve.protein import Protein, CODONS
from easme.evolve.align import determine_conservation_likelihoods
from easme.evolve.selection import k_tournament_without_replacement
from easme.triage.filter import test_protein


# PARSING PARAMETERS
STANDARD_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
# EVOLUTION PARAMETERS
NUM_GENERATIONS = 15000
NUM_PARENTS = 2
LAMBDA_NUM_CHILDREN = 10
MAX_POPULATION_SIZE = 640
MUTATION_RATE = [1]


def similarity(seq1: str, seq2: str) -> float:
    """Calculate the similarity between two sequences using SequenceMatcher."""
    return SequenceMatcher(None, seq1, seq2).ratio()


def determine_mutation_likelihoods(
    protein: Protein, comparison_proteins: list[Protein]
) -> list[float]:
    """
    Determine the mutation likelihoods for a given protein.

    Args:
        protein (Protein): The protein to determine mutation likelihoods for.
        comparison_proteins (list[Protein]): A list of proteins to compare against.

    Returns:
        list[float]: A list of mutation likelihoods for each amino acid in the protein.
    """
    # Randomly shuffle the comparison proteins to ensure diversity in the selection.
    random.shuffle(comparison_proteins)

    # Prepend the START protein to the list of comparison proteins.
    comparison_proteins = [protein] + comparison_proteins[:30]

    # Determine conservation likelihoods for the protein.
    conservation_likelihoods = determine_conservation_likelihoods(comparison_proteins)

    # Determine mutation likelihoods for the protein.
    mutation_likelihoods = [1.0 - likelihood for likelihood in conservation_likelihoods]

    return mutation_likelihoods


def directed_mutation(
    protein: Protein, mutation_likelihoods: list[float], count: int
) -> Protein:
    """
    Perform directed mutation on a protein based on mutation likelihoods.

    Args:
        protein (Protein): The protein to mutate.
        mutation_likelihoods (list[float]): The mutation likelihoods for each amino acid in the protein.
        count (int): The number of mutations to perform.

    Returns:
        Protein: The mutated protein.
    """
    # Duplicate the protein to avoid modifying the original.
    new_protein = deepcopy(protein)

    # Weighted sampling of indices to mutate based on mutation likelihoods.
    unique_indices = []
    while len(unique_indices) < count:
        # Sample indices based on mutation likelihoods.
        indices_to_mutate = random.choices(
            range(len(new_protein.amino_acid_sequence)),
            weights=mutation_likelihoods,
            k=count,
        )
        # Ensure unique indices.
        unique_indices = list(set(indices_to_mutate))
        # If we have enough unique indices, break the loop.
        if len(unique_indices) >= count:
            break

    # Mutate the selected indices in the protein's amino acid sequence.
    amino_acid_list = list(new_protein.amino_acid_sequence)
    for index in unique_indices:
        amino_acid_list[index] = random.choice(
            [aa for aa in STANDARD_AMINO_ACIDS if aa != amino_acid_list[index]]
        )
    new_protein.amino_acid_sequence = "".join(amino_acid_list)
    new_protein.detrans()

    return new_protein


def optimize_codons(protein: Protein) -> Protein:
    """
    Optimize the codons of a protein to ensure they are in the most efficient form.

    Args:
        protein (Protein): The protein to optimize.

    Returns:
        Protein: The protein with optimized codons.
    """
    # Use the freq_ecoli field of the amino acid dictionaries to optimize codons.
    optimized_sequence = []
    for aa in protein.amino_acid_sequence:
        # Get all codons where the amino_acid is the same as the current one.
        codons = [codon for codon in CODONS if codon["amino_acid"] == aa]
        # Use a weighted random choice to select a codon based on its frequency in E. coli.
        codon_choice = random.choices(
            codons,
            weights=[codon["freq_ecoli"] for codon in codons],
            k=1,
        )[0]
        optimized_sequence.append(codon_choice["codon"])
    # Join the codons to form the optimized DNA sequence.
    optimized_dna_sequence = (
        "".join(optimized_sequence).replace("U", "T").replace("u", "t")
    )
    # Create a new Protein object with the optimized DNA sequence.
    optimized_protein = deepcopy(protein)
    optimized_protein.dna_sequence_coding = optimized_dna_sequence
    optimized_protein.trans()

    return optimized_protein


def evolution(
    start_protein: list[Protein],
    mutation_likelihoods: list[float],
    generations: int,
    count: int = 1,
) -> list[Protein]:
    """
    Perform evolution on a protein over a specified number of generations.

    Args:
        start_protein (Protein): The starting protein for evolution.
        mutation_likelihoods (list[float]): The mutation likelihoods for each amino acid in the protein.
        other_fps (list[Protein]): A list of other fluorescent proteins to compare against.
        generations (int): The number of generations to evolve the protein.

    Returns:
        list[Protein]: A list of evolved proteins.
    """
    # Initialize the population with the starting protein.
    population = []
    for protein in start_protein:
        # Set the fitness of the starting protein.
        protein.fitness = test_protein(protein=protein)
        population.append(protein)

    for i in range(generations):
        # Print the current generation number.
        print(f"Generation {i + 1} running…")

        # Create children through directed mutation.
        unique_children = []
        while len(unique_children) < LAMBDA_NUM_CHILDREN:
            # Perform directed mutation on the parent to create a child.
            child = directed_mutation(
                protein=start_protein,
                mutation_likelihoods=mutation_likelihoods,
                count=count,
            )

            child.fitness = test_protein(protein=child)

            # Add the child to the population.
            if (
                child.amino_acid_sequence
                not in [p.amino_acid_sequence for p in unique_children]
                and child.amino_acid_sequence
                not in [p.amino_acid_sequence for p in population]
                and child.amino_acid_sequence != start_protein.amino_acid_sequence
            ):
                unique_children.append(child)

        # Add the unique children to the population.
        population.extend(unique_children)

        # Limit the population size to MAX_POPULATION_SIZE.
        if len(population) > MAX_POPULATION_SIZE:
            population = k_tournament_without_replacement(
                population=population,
                n=MAX_POPULATION_SIZE,
                k=NUM_PARENTS,
                elitism=True,
            )

        # Sort the population by fitness in descending order.
        population.sort(key=lambda x: x.fitness, reverse=True)

        # Optimize the codons of the surviving proteins.
        for protein in population:
            protein = optimize_codons(protein)

        # Print the best protein in the population.
        print(
            f"Best fitness: {population[0].fitness} (start fitness: {start_protein.fitness})"
        )

        # Save the population.
        with open(f"dat/out/round1/{str(count)}.json", "w", encoding="utf-8") as f:
            list_of_dicts = [p.dump_to_dict() for p in population]
            json.dump(list_of_dicts, f, indent=4)

    return population


if __name__ == "__main__":
    # Load start proteins from dat/in/formatted/sequences.json.
    with open("dat/in/formatted/sequences.json", "r", encoding="utf-8") as f:
        raw_proteins = json.load(f)

    # Create a Protein object for each start protein.
    start_proteins = []
    for protein in raw_proteins:
        start_proteins.append(
            Protein(
                name=protein["name"],
                dna_sequence_coding=protein["dna_sequence_coding"],
            )
        )

    """
    # Calculate the mutation likelihoods for pmagenta based on the red fluorescent proteins.
    mutation_likelihoods = determine_mutation_likelihoods(pmagenta, FP_RED)
    """
