from Bio import SeqIO
from Bio.Blast import NCBIXML

sequences = []
with open('data/psiblast_search.xml') as f:
    psiblast_rounds = NCBIXML.parse(f)
    # Iterate over Psiblast rounds
    for psiblast_round in psiblast_rounds:
        for alignment in psiblast_round.alignments:
            name = alignment.title.split()[0]
            if name not in sequences:
                for hsp in alignment.hsps:
                    sequences.append(name)

sequences = list(set(sequences))  # Remove duplicates
print("Number of hits retrieved usign psi-blast: {}".format(len(sequences)))