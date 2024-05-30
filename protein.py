class Protein:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.location = None

    def __str__(self):
        return f"{self.name}: {self.sequence}"

    def __len__(self):
        return len(self.sequence)