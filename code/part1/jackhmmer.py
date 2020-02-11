from models import search_jackhmmer
from references import getHumanSH2, getPositionReference, countSequences
from utils import evaluatePositionsSH2, evaluateSequencesSH2
from parsers import parseHmmerOutput
import os

#### jackhmmer #####

# Path of the multiple sequence alignment
INPUT_SEQUENCE = "../../data/sequenceP23615.fasta "

# Path where hmmer output will be saved
SEARCH_RESULT_PATH = "../../results/jackhmmsearch.hmmer_domtblout"
# Database used to perform the search
SEARCH_DB = "../../data/SwissProt_humans_reference_all.fasta"

# number of iterations
ITERATIONS = 4
# evalue for the search
EVALUE = 0.01

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

    # search with jackhmmer
    search_jackhmmer(
        result_path = SEARCH_RESULT_PATH,
        sequence_path = INPUT_SEQUENCE,
        db_path = SEARCH_DB,
        iterations = ITERATIONS,
        evalue = EVALUE
    )

    # parse search result
    jackhmm_sh2_positions = parseHmmerOutput(SEARCH_RESULT_PATH)

    # get references 
    reference_human_sh2 = getHumanSH2(PATH_SEQUENCES_HUMAN_SH2) # human sequences with sh2
    num_human_sequences = countSequences(PATH_SEQUENCES_HUMAN) # number of sequences in search dataset
    reference_sh2_positions = getPositionReference(PATH_POSITION_REFERENCE) # sh2 position infos

    # evaluate ability of retrieving sequences containing SH2
    print("\nAbility of retreiving sequences containing SH2:")
    evaluateSequencesSH2(jackhmm_sh2_positions, reference_human_sh2, num_human_sequences)

    # Evaluate ability of matching domains
    print("\nAbility of matching domains")
    evaluatePositionsSH2(jackhmm_sh2_positions, reference_sh2_positions)

    # Create "original dataset" with retrieved hits
    if SAVE_HITS:
        print("Saving hits on the original dataset")
        with open('../../datasets/original.txt', 'w') as fout:
            for seqId in jackhmm_sh2_positions.keys():
                fout.write("{}\n".format(seqId))


if __name__ == "__main__":
    # set working directory
    os.chdir(os.path.abspath(os.path.dirname(__file__)))

    main()