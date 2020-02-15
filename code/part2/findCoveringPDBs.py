import os
# set working directory
os.chdir(os.path.abspath(os.path.dirname(__file__)))  
# use functions of part 1
import sys
sys.path.insert(0,'../part1')

from parsers import parsePsiBlastOutput
from utils import createPositionSet

import pandas as pd
import urllib

PSIBLAST_PATH = '../../results/psiblast_search.txt'
PDB_PATH =  '../../datasets/pdb.cvs'
# Path to save the downloaded pdbs
COVERING_PDBS_PATH = "../../datasets/covering_pdbs"

# number of pdbs per domain position found by our model
NUM_PDBS = 1

# Threshold used to define if domains are overlappings
THRESHOLD = 0.5

# set to true to re-download everything
DOWNLOAD_PDB = False

def main():
    psiblast_sh2_positions = parsePsiBlastOutput(PSIBLAST_PATH)

    # Load pdb dataset
    pdb = pd.read_csv(PDB_PATH) 

    found_pdb_score = {}
    for seq in psiblast_sh2_positions:
        # find all pdbs for this sequence
        pdb_seq = pdb.loc[pdb['sp_primary']==seq]
        
        # sanity check
        if pdb_seq.shape[0] == 0:
            continue
        
        count_seq = 0
        for pos in psiblast_sh2_positions[seq]:
            res_df = find_overlap(pdb_seq, pos, threshold=THRESHOLD)
            if res_df.shape[0] == 0:
                # no overlap
                continue
            keys = res_df[:NUM_PDBS].pdb.values
            scores = res_df[:NUM_PDBS].overlap.values
            
            add_pdb_dict(found_pdb_score, keys, scores)
    
    print("Retrieved {} pdbs".format(len(found_pdb_score.values())))

    if DOWNLOAD_PDB:
        status = []
        for i, pdb_id in enumerate(covering_pdbs):
            print("Downloading {}/{}".format(i+1, len(covering_pdbs)))
            status.append(get_pdb(pdb_id))
            
        # retry failed download
        status = [e for e in status if e!=0]
        for i, pdb_id in enumerate(status):
            print("Downloading {}/{}".format(i+1, len(status)))
            status.append(get_pdb(pdb_id))

def find_overlap(pdb_seq, domain, threshold=0.5):
    """
    Add a column measuring the overlap score and remove 
    non overlapping pdbs
    """
    pdb_seq = pdb_seq.copy()
    pdb_seq['overlap'] = pdb_seq.apply(
        lambda row: check_overlap(row['sp_beg'], row['sp_end'], [domain]), axis=1
    )
    pdb_seq = pdb_seq.loc[pdb_seq['overlap']>threshold].sort_values('overlap', ascending=False)
    return pdb_seq


def check_overlap(pdb_start, pdb_end, domain):
    """
    Check if domain overlaps with pdb positions
    """   
    domain_pos = createPositionSet(domain)     
    pdb_pos = set(i for i in range(pdb_start,pdb_end+1))
    
    overlap = pdb_pos.intersection(domain_pos)
    
    if overlap:
        true_positive = len(overlap)
        false_negative = len(domain_pos) - len(overlap)
        false_positive = len(pdb_pos) - len(overlap)

        precision = true_positive / (true_positive + false_positive)
        sensitivity = true_positive / (true_positive + false_negative)
        f1 = 2 * (precision * sensitivity) / (precision + sensitivity)
    
        return f1 #len(overlap)/len(pdb_pos)
    else:
        return 0

def add_pdb_dict(d, pdbs, scores):
    """
    Update the score of a certain pdb if it has already
    been inserted in the dict
    """
    for i in range(len(pdbs)):
        key = pdbs[i]
        score = scores[i]
        if key in d:
            if score > d[key]:
                d[key] = score              
        else:
            d[key] = score

def get_pdb(pdb_id):
    url = "https://files.rcsb.org/download/" + pdb_id + '.pdb'  
    try: 
        response =  urllib.request.urlopen(url).read()
        # write to file
        with open("covering_pdbs/{}.pdb".format(pdb_id), 'w') as f:
            f.write(response.decode())      
    except:
        return 'Error for id '+ pdb_id

if __name__ == "__main__":   
    main()