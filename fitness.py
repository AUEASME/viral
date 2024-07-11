# Improved fitness function which simply adds 1 to summed_fitness if the amino acid is in the mutabibility list and has the highest count out of all amino acids at that position.
def calculate_fitness(protein, mutability):
    summed_fitness = 0.0
    length = len(protein.dna)
    for i in range(length):
        base_pair = protein.dna[i]
        # Check if mutability[i] has an object for this amino acid.
        base_pair_obj = next(
            (x for x in mutability[i] if x["base_pair"] == base_pair), None)
        if base_pair_obj and base_pair_obj["count"] == max([x["count"] for x in mutability[i]]):
            summed_fitness += 1.0

    return summed_fitness / length
