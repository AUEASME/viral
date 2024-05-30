import random
from fitness import calculate_fitness


AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


class Protein:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.location = None
        self.fitness = calculate_fitness(sequence)
        self.descendants = []

    def __str__(self):
        return f"{self.name}: {self.sequence}"
    
    def __len__(self):
        return len(self.sequence)
    
    def mutate_sequence(self):
        # Select a random amino acid in the sequence.
        index = random.randint(0, len(self.sequence) - 1)
        sequence = list(self.sequence)
        # Change the amino acid.
        sequence[index] = random.choice(AMINO_ACIDS)
        self.sequence = "".join(sequence)
        # Update the fitness.
        self.fitness = calculate_fitness(self.sequence)

    def select_random_node(self):
        if not self.descendants:
            return self
        else:
            return random.choice([self] + [descendant.select_random_node() for descendant in self.descendants])
        
    def calculate_all_fitness(self):
        self.fitness = calculate_fitness(self.sequence)
        for descendant in self.descendants:
            descendant.calculate_all_fitness()

    def mutate(self):
        random_node = self.select_random_node()
        duplicate = Protein(random_node.name, random_node.sequence)
        duplicate.mutate_sequence()
        if (duplicate.fitness > random_node.fitness):
            self.descendants.append(duplicate)
            print("New descendant added.")