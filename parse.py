import json
from protein import Protein


def parse_proteins(input_file):
    data = None

    # Parse JSON file.
    with open(input_file, "r") as file:
        data = json.load(file)

    # Create Protein objects.
    proteins = []
    for protein in data:
        new_protein = Protein(protein["name"], protein["sequence"])
        new_protein.location = protein["location"]
        proteins.append(new_protein)

    return proteins


def parse(input_file):
    data = None

    with open(input_file, "r") as file:
        data = json.load(file)

    return data


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

    with open("data/prototype.txt", "w+") as file:
        file.write("".join(most_common_amino_acids))


def determine_mutability(proteins):
    # Remove proteins that are shorter than the longest variant (we'll figure out how to handle them later).
    max_length = max([len(protein) for protein in proteins])
    proteins = [protein for protein in proteins if len(protein) == max_length]

    # For each index, we want to count how many unique amino acids have been seen, and how many times they've been seen.
    mutability = []
    for i in range(len(proteins[0].sequence)):
        amino_acids = []
        for protein in proteins:
            amino_acid = protein.sequence[i]
            # Check if an object for this amino acid already exists.
            amino_acid_obj = next(
                (x for x in amino_acids if x["amino_acid"] == amino_acid), None)
            if amino_acid_obj:
                amino_acid_obj["count"] += 1
            else:
                amino_acids.append({"amino_acid": amino_acid, "count": 1})

        # Sort the amino acids by count.
        amino_acids.sort(key=lambda x: x["count"], reverse=True)
        mutability.append(amino_acids)

    with open("data/json/mutability.json", "w+") as file:
        json.dump(mutability, file, indent=4)


def add_locations(proteins):
    # Load the location data.
    with open("data/json/locations.json", "r") as file:
        locations = json.load(file)

        for location in locations:
            for protein in proteins:
                if protein.name in location["accessions"]:
                    protein.location = location["location"]

    with open("data/json/sequences.json", "w+") as file:
        json.dump([{"name": protein.name, "sequence": protein.sequence,
                  "location": protein.location} for protein in proteins], file, indent=2)


def add_alabamas(proteins):
    for protein in proteins:
        if protein.location == None:
            protein.location = "Alabama"

    with open("data/json/sequences.json", "w+") as file:
        json.dump([{"name": protein.name, "sequence": protein.sequence,
                  "location": protein.location} for protein in proteins], file, indent=2)


if __name__ == "__main__":
    proteins = parse_proteins("data/json/sequences.json")
    # get_prototype(proteins)
    # determine_mutability(proteins)
    # add_locations(proteins)
    add_alabamas(proteins)
