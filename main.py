from parse import parse_proteins

PROTEINS = parse_proteins("data/json/sequences.json")


def main():
    initial_population = []
    for protein in PROTEINS:
        if protein.location == "Alabama":
            initial_population.append(protein)

    for protein in initial_population:
        for other_protein in initial_population:
            if other_protein.dna == protein.dna:
                # Remove other_protein from the list.
                initial_population.remove(other_protein)

    for _ in range(10000):
        for protein in initial_population:
            protein.mutate()

    for protein in initial_population:
        protein.save_to_json(f"data/out/{protein.name}.json")


if __name__ == "__main__":
    main()
