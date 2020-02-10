
def parsePsiBlastOutput(output_path):
    """
    Parse psi-blast output files psiblast_search.txt
    Return dict of the type {'P16885':[{'start': 532, 'end': 617}], ...}
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


def parseHmmerOutput(output_path):
    """
    Parse hmmserach output file hmmsearch.hmmer_domtblout
    Return dict of the type {'P16885':[{'start': 532, 'end': 617}], ...}
    """
    hmm_sh2_positions = {}
    with open(output_path) as f:
        for line in f:
            if line[0] != "#" and line[0] != "-":
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