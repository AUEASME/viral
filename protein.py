import json
import random
from fitness import calculate_fitness
from uuid import uuid4


BASE_PAIRS = "ACGT"
MUTABILITY = json.load(open("data/json/mutability.json"))


class Protein:
    def __init__(self, name, dna, sequence):
        self.name = name
        self.dna = dna
        self.sequence = sequence
        self.location = None
        self.fitness = calculate_fitness(self, MUTABILITY)
        self.descendants = []

    def __str__(self):
        return f"{self.name}: {self.dna}"

    def __len__(self):
        return len(self.dna)

    def mutate_dna(self):
        # Select a random amino acid in the sequence.
        index = random.randint(0, len(self.dna) - 1)
        dna = list(self.dna)
        # Change the amino acid.
        dna[index] = random.choice(BASE_PAIRS)
        self.dna = "".join(dna)
        # Update the fitness.
        self.fitness = calculate_fitness(self, MUTABILITY)

    def select_random_node(self):
        if not self.descendants:
            return self
        else:
            return random.choice(
                [self]
                + [descendant.select_random_node() for descendant in self.descendants]
            )

    def calculate_all_fitness(self):
        self.fitness = calculate_fitness(self, MUTABILITY)
        for descendant in self.descendants:
            descendant.calculate_all_fitness()

    def mutate(self):
        random_node = self.select_random_node()
        duplicate = Protein(str(uuid4()), random_node.sequence)
        duplicate.mutate_dna()
        if (duplicate.fitness > random_node.fitness) and (
            duplicate.sequence
            not in [descendant.sequence for descendant in random_node.descendants]
        ):
            random_node.descendants.append(duplicate)
            print("New descendant added.")

    def dump_to_json(self):
        return {
            "name": self.name,
            "dna": self.dna,
            "sequence": self.sequence,
            "location": self.location,
            "fitness": self.fitness,
            "descendants": [
                descendant.dump_to_json() for descendant in self.descendants
            ],
        }

    def save_to_json(self, output_file):
        with open(output_file, "w+") as file:
            json.dump(self.dump_to_json(), file, indent=2)
