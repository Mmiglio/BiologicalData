from models import create_profile, search_psiblast
from references import getHumanSH2, getPositionReference, countSequences
from utils import evaluatePositionsSH2, evaluateSequencesSH2

# Profile + psi-blast parameters
PROFILE_PATH = "../models/profile.pssm"
INPUT_DATA_PATH = "../data/BLAST_uniref90.fasta"
MSA_PATH = "../data/PF00017_seed.fasta"

SEARCH_RESULT_PATH = "../results/psiblast_search.txt"
SEARCH_DB = "../data/SwissProt_humans_reference_all.fasta"

# References dataset path
PATH_SEQUENCES_SH2 = "../data/SwissProt_humans_reference.fasta"
PATH_SEQUENCES_HUMAN = "../data/SwissProt_humans_reference_all.fasta"
PATH_POSITION_REFERENCE = "../data/interpo-PF00017.json"


def main():
    # Create pssm
    create_profile(
        profile_path = PROFILE_PATH,
        data_path = INPUT_DATA_PATH,
        msa_path = MSA_PATH
    )

    # use pssm to search using psiblast
    search_psiblast(
        result_path = SEARCH_RESULT_PATH,
        pssm_path = PROFILE_PATH,
        db_path = SEARCH_DB
    )

    # parse search result
    psiblast_sh2_positions = parsePsiBlastOutput(SEARCH_RESULT_PATH)

    # get references 
    reference_human_sh2 = getHumanSH2(PATH_SEQUENCES_SH2) # human sequences with sh2
    num_human_sequences = countSequences(PATH_SEQUENCES_HUMAN) # number of sequences in search dataset
    reference_sh2_positions = getPositionReference(PATH_POSITION_REFERENCE) # sh2 position infos

    # evaluate ability of retrieving sequences containing SH2
    print("\nAbility of retreiving sequences containing SH2:")
    evaluateSequencesSH2(psiblast_sh2_positions, reference_human_sh2, num_human_sequences)

    # Evaluate ability of matching domains
    print("\nAbility of matching domains")
    evaluatePositionsSH2(psiblast_sh2_positions, reference_sh2_positions)


def parsePsiBlastOutput(output_path=SEARCH_RESULT_PATH):
    """
    return dict of the type {'P16885':[{'start': 532, 'end': 617}], ...}
    """
    psiblast_sh2_positions = {}
    with open(output_path, 'r') as f:
        for line in f:
            if len(line)>1:
                qseqid, sseqid, pident, length, mismatch, gapopen, \
                qstart, qend, sstart, send, evalue, bitscore = line.strip().split()
                if sseqid not in psiblast_sh2_positions:
                    psiblast_sh2_positions[sseqid] = [{'start':int(sstart), 'end':int(send)}]
                else:
                    pos = {'start':int(sstart), 'end':int(send)}
                    # check if the position has alredy been inserted
                    # otherwise insert it in the dictionary
                    if pos in psiblast_sh2_positions[sseqid]:
                        # position already inserted
                        continue
                    else:
                        psiblast_sh2_positions[sseqid].append(pos)
            else:
                break
                
    print("{} sequences found with psi-blast".format(len(psiblast_sh2_positions.keys())))
    return psiblast_sh2_positions

if __name__ == "__main__":
    main()
