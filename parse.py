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
        if len(protein["dna"]) == 774:
            new_protein = Protein(protein["name"], protein["dna"], protein["sequence"])
            new_protein.location = protein["location"]
            proteins.append(new_protein)

    return proteins


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
        base_pairs = []
        for protein in proteins:
            base_pair = protein.dna[i]
            # Check if an object for this amino acid already exists.
            amino_acid_obj = next(
                (x for x in base_pairs if x["base_pair"] == base_pair), None
            )
            if amino_acid_obj:
                amino_acid_obj["count"] += 1
            else:
                base_pairs.append({"base_pair": base_pair, "count": 1})

        # Sort the amino acids by count.
        base_pairs.sort(key=lambda x: x["count"], reverse=True)
        mutability.append(base_pairs)

    with open("data/json/mutability.json", "w+") as file:
        json.dump(mutability, file, indent=2)


if __name__ == "__main__":
    print("Parsing data...")
    proteins = parse_proteins("data/json/sequences.json")

    print("Determining mutability...")
    determine_mutability(proteins)
