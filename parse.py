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
        new_protein = Protein(protein["name"], protein["dna"], protein["sequence"])
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
    # Determine mode length of proteins.
    mode_length = max(
        set([len(protein) for protein in proteins]),
        key=[len(protein) for protein in proteins].count,
    )
    # Remove proteins that are not the mode length.
    proteins = [protein for protein in proteins if len(protein) == mode_length]

    # For each index, we want to count how many unique base pairs/amino acids have been seen, and how many times they've been seen.
    mutability = []
    for i in range(len(proteins[0].dna)):
        amino_acids = []
        for protein in proteins:
            base_pair = protein.dna[i]
            # Check if an object for this amino acid already exists.
            amino_acid_obj = next(
                (x for x in amino_acids if x["base_pair"] == base_pair), None
            )
            if amino_acid_obj:
                amino_acid_obj["count"] += 1
            else:
                amino_acids.append({"base_pair": base_pair, "count": 1})

        # Sort the amino acids by count.
        amino_acids.sort(key=lambda x: x["count"], reverse=True)
        mutability.append(amino_acids)

    with open("data/json/mutability.json", "w+") as file:
        json.dump(mutability, file, indent=2)


def add_locations(proteins):
    # Load the location data.
    with open("data/json/locations.json", "r") as file:
        locations = json.load(file)

        for location in locations:
            for protein in proteins:
                if protein.name in location["accessions"]:
                    protein.location = location["location"]

    with open("data/json/sequences.json", "w+") as file:
        json.dump(
            [
                {
                    "name": protein.name,
                    "sequence": protein.sequence,
                    "location": protein.location,
                }
                for protein in proteins
            ],
            file,
            indent=2,
        )


def add_alabamas(proteins):
    for protein in proteins:
        if protein.location == None:
            protein.location = "Alabama"

    with open("data/json/sequences.json", "w+") as file:
        json.dump(
            [
                {
                    "name": protein.name,
                    "sequence": protein.sequence,
                    "location": protein.location,
                }
                for protein in proteins
            ],
            file,
            indent=2,
        )


def fix_dna():
    # Load existing protein data.
    proteins = json.load(open("data/json/sequences.json", "r"))

    # Add a dna field to each protein.
    for protein in proteins:
        protein["dna"] = ""

    # Load the DNA data.
    dna = None
    with open("data/json/dna.txt", "r") as file:
        dna = file.read()

    # Split the DNA data.
    dna = dna.split()
    print(dna)

    target = dna[0].replace(">", "")
    for i in range(len(dna)):
        if dna[i].startswith(">"):
            target = dna[i].replace(">", "")
        else:
            # Find protein with same name as target.
            for protein in proteins:
                if protein["name"] == target:
                    protein["dna"] += dna[i]

    with open("data/json/sequences.json", "w+") as file:
        json.dump(proteins, file, indent=2)


if __name__ == "__main__":
    print("Parsing data...")
    proteins = parse_proteins("data/json/sequences.json")

    print("Determining mutability...")
    determine_mutability(proteins)
    # add_locations(proteins)
    # add_alabamas(proteins)
