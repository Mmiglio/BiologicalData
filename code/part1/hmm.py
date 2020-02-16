from models import create_hmm, search_hmm
from references import getHumanSH2, getPositionReference, countSequences
from utils import evaluatePositionsSH2, evaluateSequencesSH2
from parsers import parseHmmerOutput
import os

#### Hmmer model  #####

# Path where generated model will be stored
MODEL_PATH = "../../models/hmm_model.hmm"
# Path of the multiple sequence alignment
MSA_PATH = "../../data/msa_edited.fasta"

# Path where hmmer output will be saved
SEARCH_RESULT_PATH = "../../results/hmmsearch.hmmer_domtblout"
# Database used to perform the search
SEARCH_DB = "../../data/SwissProt_humans_reference_all.fasta"

# save retrieved hits on the "original dataset"
SAVE_HITS = False

#### Reference datasets ####
# File containing reviewed human sequences with SH2 domain
PATH_SEQUENCES_HUMAN_SH2 = "../../data/SwissProt_humans_reference.fasta"
# File containing all eviewed human sequences
PATH_SEQUENCES_HUMAN = "../../data/SwissProt_humans_reference_all.fasta"
# Json containing information about the position of the domains
PATH_POSITION_REFERENCE = "../../data/interpo-PF00017.json"


def main():
    # Create hmm
    create_hmm(
        model_path = MODEL_PATH,
        msa_path = MSA_PATH
    )

    # search with hmmer
    search_hmm(
        result_path = SEARCH_RESULT_PATH,
        model_path = MODEL_PATH,
        db_path = SEARCH_DB
    )

    # parse search result
    hmm_sh2_positions = parseHmmerOutput(SEARCH_RESULT_PATH)

    # get references 
    reference_human_sh2 = getHumanSH2(PATH_SEQUENCES_HUMAN_SH2) # human sequences with sh2
    num_human_sequences = countSequences(PATH_SEQUENCES_HUMAN) # number of sequences in search dataset
    reference_sh2_positions = getPositionReference(PATH_POSITION_REFERENCE) # sh2 position infos

    # evaluate ability of retrieving sequences containing SH2
    print("\nAbility of retreiving sequences containing SH2:")
    resulting_metrics = evaluateSequencesSH2(
        hmm_sh2_positions, reference_human_sh2, num_human_sequences
        )

    # Evaluate ability of matching domains
    print("\nAbility of matching domains")
    evaluatePositionsSH2(hmm_sh2_positions, reference_sh2_positions)

    # Create "original dataset" with retrieved hits
    if SAVE_HITS:
        print("Saving hits on the original dataset")
        with open('../datasets/original.txt', 'w') as fout:
            for seqId in hmm_sh2_positions.keys():
                fout.write("{}\n".format(seqId))


if __name__ == "__main__":
    # set working directory
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

    main()