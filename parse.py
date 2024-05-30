import json


class Protein:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence

    def __str__(self):
        return f"{self.name}: {self.sequence}"

    def __len__(self):
        return len(self.sequence)


def parse(input_file):
    sequences = []
    with open(input_file) as file:
        sequences = file.readlines()

    proteins = []
    for i in range(int(len(sequences) / 2)):
        # print(sequences[2 * i], sequences[2 * i + 1])
        protein = Protein(
            sequences[2 * i][1:].replace("\n", ""), sequences[2 * i + 1].strip())
        proteins.append(protein)

    # Dump to JSON.
    filename = input_file.split('.')[0]
    with open(f'{filename}.json', 'w+') as json_file:
        json.dump(proteins, json_file, default=lambda x: x.__dict__, indent=4)

    return proteins


def get_prototype(proteins):
    most_common_amino_acids = []
    for i in range(len(proteins[0].sequence)):
        # Iterate over each protein, getting the most common amino acid at each position.
        amino_acids = {}
        for protein in proteins:
            try:
                amino_acid = protein.sequence[i]
                if amino_acid in amino_acids:
                    amino_acids[amino_acid] += 1
                else:
                    amino_acids[amino_acid] = 1
            except IndexError:
                pass

        # Get the most common amino acid.
        most_common_amino_acid = max(amino_acids, key=amino_acids.get)
        most_common_amino_acids.append(most_common_amino_acid)

    with open("prototype.txt", "w+") as file:
        file.write("".join(most_common_amino_acids))


def determine_mutability(proteins):
    # Remove proteins that are shorter than the longest variant (we'll figure out how to handle them later).
    max_length = max([len(protein) for protein in proteins])
    proteins = [protein for protein in proteins if len(protein) == max_length]

    # For each index, we want to count how many unique amino acids have been seen.
    mutability = []
    for i in range(len(proteins[0].sequence)):
        amino_acids = set()
        for protein in proteins:
            try:
                amino_acid = protein.sequence[i]
                amino_acids.add(amino_acid)
            except IndexError:
                pass

        mutability.append(list(amino_acids))

    with open("mutability.json", "w+") as file:
        json.dump(mutability, file, indent=4)


if __name__ == "__main__":
    proteins = parse("sequences.txt")
    get_prototype(proteins)
    determine_mutability(proteins)
