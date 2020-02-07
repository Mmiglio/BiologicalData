from models import create_hmm, search_hmm
from references import getHumanSH2, getPositionReference, countSequences
from utils import evaluatePositionsSH2, evaluateSequencesSH2

# Hmmer model parameters
MODEL_PATH = "../models/hmm_model.hmm"
MSA_PATH = "../data/PF00017_seed.fasta"

SEARCH_RESULT_PATH = "../results/hmmsearch.hmmer_domtblout"
SEARCH_DB = "../data/SwissProt_humans_reference_all.fasta"

# References dataset path
PATH_SEQUENCES_SH2 = "../data/SwissProt_humans_reference.fasta"
PATH_SEQUENCES_HUMAN = "../data/SwissProt_humans_reference_all.fasta"
PATH_POSITION_REFERENCE = "../data/interpo-PF00017.json"


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
    reference_human_sh2 = getHumanSH2(PATH_SEQUENCES_SH2) # human sequences with sh2
    num_human_sequences = countSequences(PATH_SEQUENCES_HUMAN) # number of sequences in search dataset
    reference_sh2_positions = getPositionReference(PATH_POSITION_REFERENCE) # sh2 position infos

    # evaluate ability of retrieving sequences containing SH2
    print("\nAbility of retreiving sequences containing SH2:")
    evaluateSequencesSH2(hmm_sh2_positions, reference_human_sh2, num_human_sequences)

    # Evaluate ability of matching domains
    print("\nAbility of matching domains")
    evaluatePositionsSH2(hmm_sh2_positions, reference_sh2_positions)


def parseHmmerOutput(output_path=SEARCH_RESULT_PATH):
    """
    return dict of the type {'P16885':[{'start': 532, 'end': 617}], ...}
    """
    hmm_sh2_positions = {}
    with open(SEARCH_RESULT_PATH) as f:
        for line in f:
            if line[0] != "#":
                target, tacc, tlen, query, qacc, qlen, \
                fevalue, fscore, fbias, _, _, dcevalue, dievalue, dscore, dbias, _, _ , \
                alifrom, alito, evfrom, envto, accuracy, desc = line.strip().split()[:23]
                tacc = target.split("|")[1]

                if tacc not in hmm_sh2_positions:
                    hmm_sh2_positions[tacc] = [{'start':int(alifrom), 'end':int(alito)}]
                else:
                    pos = {'start':int(alifrom), 'end':int(alito)}
                    # check if the position has alredy been inserted
                    # otherwise insert it in the dictionary
                    if pos in hmm_sh2_positions[tacc]:
                        continue
                    else:
                        hmm_sh2_positions[tacc].append(pos)
    print("{} sequences found with hmm-search".format(len(hmm_sh2_positions.keys())))
    return hmm_sh2_positions

if __name__ == "__main__":
    main()
