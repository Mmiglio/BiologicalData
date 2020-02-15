from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import subprocess
import numpy as np
import re
import os

PDBS_PATH = "../../datasets/covering_pdbs"
TMALIGN_PATH = "./TMalign/TMalign"
SCORE = 'tmscore' # tmscore or rmsd

def main():
    # read list of PDBs that will be used
    with open(PDBS_PATH+"/list_pdbs.txt", 'r') as f:
        pdbs = f.read().split('\n')[:-1]
        print("Read {} pdbs".format(len(pdbs)))

    alignments_score = np.zeros((len(pdbs), len(pdbs)))
    # performe pairwise alignment (upper triangular matrix)
    for i, pdb1 in enumerate(pdbs):
        if i%10==0:
            print("row {} / {}".format(i+1, len(pdbs)))
        for j, pdb2 in enumerate(pdbs):
            if j<i:
                continue
            alignments_score[i][j] = TMalign(pdb1, pdb2, score=SCORE)

    if SCORE == "tmscore":
        # set 0 as minimum distance instead of 1
        iu = np.triu_indices(alignments_score.shape[0])
        alignments_score[iu] = np.abs(alignments_score[iu] - 1)

    fig = plt.figure(figsize=(12,6))
    labelList = [pdb.split('.')[0] for pdb in pdbs]
    dendrogram(linkage(alignments_score, method='average'), labels=labelList, leaf_rotation = 90)
    plt.title("Hierarchical clustering based on {}".format("TM-score" if SCORE == "tmscore" else "RMSD"))
    plt.ylabel('Distance')
    plt.tight_layout()
    plt.savefig("../../figures/dendrogram_{}.pdf".format(SCORE))
    plt.show()


def TMalign(pdb1, pdb2, score='rmsd'):
    """
    Perform alignment usign TMalign and parse output
    """
    cmd = TMALIGN_PATH +" {0}/{1} {0}/{2}".format(PDBS_PATH, pdb1, pdb2)
    results = subprocess.run(
        cmd, shell=True, universal_newlines=True, 
        stdout=subprocess.PIPE
        )
    res = results.stdout
    # parse TMalign output
    if score == 'rmsd':
        return float(re.findall(r'RMSD=[\s]*([\d.]*\d+)', res)[0])
    elif score == 'tmscore':
        return float(re.findall(r'TM-score=[\s]*([\d.]*\d+)', res)[0])
    else:
        print("Invalid score! use either 'rmsd' or 'tmscore'")

if __name__ == "__main__":
    # set working directory
    os.chdir(os.path.abspath(os.path.dirname(__file__)))  
    main()