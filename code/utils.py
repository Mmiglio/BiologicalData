import math

def computeMetrics(true_positive, true_negative, false_positive, false_negative):
    """
    Compute a set of metrics starting from true positives, true negatives etc.
    """
    accuracy = (true_positive + true_negative) / (true_positive + true_negative + false_negative + false_positive)
    precision = true_positive / (true_positive + false_positive)
    sensitivity = true_positive / (true_positive + false_negative)
    specificity = true_negative / (true_negative + false_positive)
    f1 = 2 * (precision * sensitivity) / (precision + sensitivity)
    mcc = (true_positive*true_negative - false_positive*false_negative) / (math.sqrt(
        (true_positive + false_positive) * (true_positive + false_negative) * 
        (true_negative + false_positive) * (true_negative + false_negative)
    ))
    
    print("\t*Accuracy: {:.2f}".format(accuracy))
    print("\t*Precision: {:.2f}".format(precision))
    print("\t*Sensitivity: {:.2f}".format(sensitivity))
    print("\t*Specificity: {:.2f}".format(specificity))
    print("\t*F1 score: {:.2f}".format(f1))
    print("\t*MCC: {:.2f}".format(mcc)) 
    print("")

def createPositionSet(positions_list):
    """
    position_list = [{'start': '532', 'end': '617'}, {'start': '646', 'end': '720'}]
    return a list {532, 533, ...., 617, 646, ..., 720}
    """
    positions = set()
    for pos in positions_list:
        p = [i for i in range(pos['start'], pos['end']+1)]
        positions = positions.union(p)
    return positions

def evaluateSequencesSH2(retrieved_sequences, ground_truth, num_sequences):
    """
    retrieved_sequences = dict of the type {'P16885':[{'start': 532, 'end': 617}], ...}
    ground_truth = set of human sequences containing SH2 (set_human_sh2 in this notebook)
    """
    retrieved_sequences = set(retrieved_sequences.keys())
    
    true_positive = len(set(retrieved_sequences).intersection(ground_truth))
    # If there aren't true postive later there will be a division by 0
    if true_positive == 0:
        print("No true positive!")
        return
    
    false_negative = len(ground_truth) - true_positive
    false_positive = len(set(retrieved_sequences)) - true_positive
    true_negative = num_sequences - true_positive - false_positive - false_negative
    
    computeMetrics(
        true_positive,
        true_negative,
        false_positive,
        false_negative
    )  

def evaluatePositionsSH2(predicted_positions_sh2, reference_positions_sh2):
    """
    predicted_positions_sh2: dict of the type {'P16885':[{'start': 532, 'end': 617}], ...}
    reference_positions_sh2 {'P16885':{'length': 1265, 'positions': [{'start': 532, 'end': 617}]}, ...}
    """
    list_tp = []
    list_fp = []
    list_fn = []
    list_tn = []

    for seqid in predicted_positions_sh2:

        if seqid in reference_positions_sh2:
            sequence_length = reference_positions_sh2[seqid]['length']
            reference_positions = createPositionSet(reference_positions_sh2[seqid]['positions'])
            identified_positions = createPositionSet(predicted_positions_sh2[seqid])    

            overlap  = reference_positions.intersection(identified_positions)

            true_positive = len(overlap)
            false_positive = len(identified_positions.difference(overlap))
            false_negative = len(reference_positions.difference(overlap))
            true_negative  = sequence_length - true_positive - false_positive - false_negative

            list_tp.append(true_positive)
            list_fp.append(false_positive)
            list_fn.append(false_negative)
            list_tn.append(true_negative)
            
    computeMetrics(
        true_positive=sum(list_tp),
        true_negative=sum(list_tn),
        false_positive=sum(list_fp),
        false_negative=sum(list_fn)
    )