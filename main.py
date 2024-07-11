from easme.parse import parse_proteins, determine_mutability
from easme.fitness import calculate_fitness_from_consensus

PROTEINS = parse_proteins("data/json/sequences.json", 774)
MUTABILITY = determine_mutability(PROTEINS)


def main():
    initial_population = []
    for protein in PROTEINS:
        if protein.location_of_origin == "Alabama":
            initial_population.append(protein)

    for protein in initial_population:
        for other_protein in initial_population:
            if other_protein.dna_sequence_coding == protein.dna_sequence_coding:
                # Remove other_protein from the list.
                initial_population.remove(other_protein)

    for protein in initial_population:
        protein.fitness = calculate_fitness_from_consensus(protein, MUTABILITY)

    for _ in range(10000):
        for protein in initial_population:
            new_individual = protein.mutate(MUTABILITY)
            if new_individual:
                print("New individual added!")

    for protein in initial_population:
        protein.save_to_json(f"data/out/{protein.name}.json")


if __name__ == "__main__":
    main()
