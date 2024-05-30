SEQUENCE_LENGTH = 358


'''
def calculate_fitness(protein, mutability):
    summed_fitness = 0.0
    for i in range(len(protein.sequence)):
        amino_acid = protein.sequence[i]
        # Check if mutability[i] has an object for this amino acid.
        amino_acid_obj = next(
            (x for x in mutability[i] if x["amino_acid"] == amino_acid), None)
        if amino_acid_obj:
            summed_fitness += (amino_acid_obj["count"] / SEQUENCE_LENGTH)

    return summed_fitness / SEQUENCE_LENGTH
'''


# Improved fitness function which simply adds 1 to summed_fitness if the amino acid is in the mutabibility list and has the highest count out of all amino acids at that position.
def calculate_fitness(protein, mutability):
    summed_fitness = 0.0
    for i in range(len(protein.sequence)):
        amino_acid = protein.sequence[i]
        # Check if mutability[i] has an object for this amino acid.
        amino_acid_obj = next(
            (x for x in mutability[i] if x["amino_acid"] == amino_acid), None)
        if amino_acid_obj and amino_acid_obj["count"] == max([x["count"] for x in mutability[i]]):
            summed_fitness += 1.0

    return summed_fitness / SEQUENCE_LENGTH
