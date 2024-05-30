from fitness import calculate_fitness
from parse import parse, parse_proteins
from protein import PhylogeneticTree

SEQUENCES = "data/json/sequences.json"
MUTABILITY = "data/json/mutability.json"


def main():
    proteins = parse_proteins(SEQUENCES)
    mutability = parse(MUTABILITY)

    initial_population = []
    for protein in proteins:
        if protein.location == "Alabama":
            initial_population.append(protein)

    for protein in initial_population:
        for other_protein in initial_population:
            if other_protein.sequence == protein.sequence:
                # Remove other_protein from the list.
                initial_population.remove(other_protein)


if __name__ == "__main__":
    main()
