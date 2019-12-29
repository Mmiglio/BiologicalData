sequences_psiblast = []
with open("data/psiblast_search.txt") as f:
    for line in f:
        if line:
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore = line.strip().split()
            if sseqid not in sequences_psiblast:
                sequences_psiblast.append(sseqid)
print("Number of hits retrieved usign psi-blast: {}".format(len(sequences_psiblast)))

sequences_hmm = []
with open("data/hmmsearch.hmmer_domtblout") as f:
    for line in f:
        if line[0] != "#":
            target, tacc, tlen, query, qacc, qlen, \
            fevalue, fscore, fbias, _, _, dcevalue, dievalue, dscore, dbias, _, _ , \
            alifrom, alito, evfrom, envto, accuracy, desc = line.strip().split()[:23]
            tacc = target.split("|")[1]
            sequences_hmm.append(tacc)
print("Number of hits retrieved usign hmmsearch: {}".format(len(sequences_hmm)))
print