from Bio import SeqIO
import pandas as pd
import os

# original dataset
PATH_ORIGINAL_DB = '../../datasets/original.txt'
# datasate with human sequences in swissprot
PATH_REFERENCE_DB = '../../data/SwissProt_humans_reference_all.fasta'
# map from pdb chains to uniprot entries
PATH_PDB_UNIPROT_REL = '../../data/pdb_chain_uniprot.tsv'

# output dataset: save pdb dataset
OUTPUT_PATH = '../../datasets/pdb.cvs'

def main():

    list_human = getListHumans(PATH_REFERENCE_DB)

    # Import original dataset
    original_proteins = []
    with open(PATH_ORIGINAL_DB) as file:
        for line in file:
            original_proteins.append(line.strip())

    # Import pdb - uniprot relation file
    pdb_rel = pd.read_csv(PATH_PDB_UNIPROT_REL, sep = '\t', header = 1)
    pdb_rel.columns = list(map(lambda x: x.lower(), pdb_rel.columns.values))
    
    pdb_dataset = pdb_rel.loc[pdb_rel.sp_primary.isin(original_proteins),['pdb','sp_primary','chain','sp_beg','sp_end']].copy()

    # Create final version of the PDB dataset: original proteins in PDB + human proteins with same chain as the originals
    frames = [pdb_dataset, pdb_rel[(pdb_rel.pdb.isin(pdb_dataset.pdb)) & 
              ~(pdb_rel.sp_primary.isin(pdb_dataset.sp_primary))
              & (pdb_rel.sp_primary.isin(list_human))]
              ]
    pdb_dataset = pd.concat(frames, sort = False)[
        ['pdb','sp_primary','chain','sp_beg','sp_end']
        ].reset_index(drop=True)

    pdb_dataset.to_csv(OUTPUT_PATH, index=False)
    print("Dataset saved as {}".format(OUTPUT_PATH))


def getListHumans(path):
    """
    Return the list of human sequences in swissprot
    """
    human = SeqIO.parse(path, 'fasta')
    list_human = []
    for sequence in human:
        name = sequence.id # name is in the form sp|P46108|CRK_HUMAN
        list_human.append(name.split('|')[1])

    print('There are {} human proteins in SwissProt'.format(len(list_human)))
    return list_human


if __name__ == "__main__":
    # set working directory
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
    
    main()