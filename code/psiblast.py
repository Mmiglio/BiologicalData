from models import create_profile, search_psiblast
from parsers import parsePsiBlastOutput
from references import getHumanSH2, getPositionReference, countSequences
from utils import evaluatePositionsSH2, evaluateSequencesSH2

#### Profile + psi-blast ####

# Path where generated profile will be stored
PROFILE_PATH = "../models/profile.pssm"
# Path of the dataset used to create the profile
INPUT_DATA_PATH = "../data/BLAST_uniref90.fasta"
# Path of the multiple sequence alignment
MSA_PATH = "../data/PF00017_seed.fasta"

# Path where psi-blast output will be saved
SEARCH_RESULT_PATH = "../results/psiblast_search.txt"
# Database used to perform the search
SEARCH_DB = "../data/SwissProt_humans_reference_all.fasta"

#### Reference datasets ####

# File containing reviewed human sequences with SH2 domain
PATH_SEQUENCES_HUMAN_SH2 = "../data/SwissProt_humans_reference.fasta"
# File containing all eviewed human sequences
PATH_SEQUENCES_HUMAN = "../data/SwissProt_humans_reference_all.fasta"
# Json containing information about the position of the domains
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
    reference_human_sh2 = getHumanSH2(PATH_SEQUENCES_HUMAN_SH2) # human sequences with sh2
    num_human_sequences = countSequences(PATH_SEQUENCES_HUMAN) # number of sequences in search dataset
    reference_sh2_positions = getPositionReference(PATH_POSITION_REFERENCE) # sh2 position infos

    # evaluate ability of retrieving sequences containing SH2
    print("\nAbility of retreiving sequences containing SH2:")
    evaluateSequencesSH2(psiblast_sh2_positions, reference_human_sh2, num_human_sequences)

    # Evaluate ability of matching domains
    print("\nAbility of matching domains")
    evaluatePositionsSH2(psiblast_sh2_positions, reference_sh2_positions)

if __name__ == "__main__":
    main()