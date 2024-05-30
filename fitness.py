SEQUENCE_LENGTH = 358


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
