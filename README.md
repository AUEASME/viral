# EASME TSWV Project

This is a project by the EASME team to analyze and predict the evolution of the [tomato spotted wilt virus nucleocapsid protein](https://www.uniprot.org/uniprotkb/P25999/entry) [in Alabama](https://www.ncbi.nlm.nih.gov/nucleotide/PQ369628.1).

Make sure you have the EASME core library installed in your Python `Lib` folder!

## How It Works

The core evolutionary loop begins with an initial population of TSWV individuals found in Alabama, then grows new individuals from those. Reproduction is done through single-point crossover between two individuals, and new children are mutated using a simple (but biologically-accurate) substitution operator. New individuals that have a higher fitness are appended to the phylogenetic trees of their parents. Fitness is judged by similarity to a consensus sequence generated from the most common base pairs for each sequence position in the input data.

## Misc. Notes

### BLAST Types

https://sequenceserver.com/blog/choosing-blast-algorithms/

- **BLASTN:** nucleotide BLAST (DNA)
- **BLASTP:** protein BLAST (amino acids)
- **BLASTX:** translates from nucleotides to amino acids, then searches
- **TBLASTX:** "This translates both query and database sequences into proteins in all 6 frames, then aligns them."
