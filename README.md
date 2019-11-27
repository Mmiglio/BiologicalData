# BiologicalData

short description

## Steps

### 1: Retrieve homologous proteins in UniProt

Perform a [BLAST search](https://www.uniprot.org/blast/) with our input sequence 

>YYFPFNGRQAEDYLRSKERGEFVIRQSSRGDDHLVITWKLDKDLFQHIDIQELEKENPLALGKVLIVDNQKYNDLDQIIVEY

on UniProt (Should we use a different database? e.g. UniProtKB/humans, UniprotKB/swiss-prot, UniRef90, ...). I've used `BLOSUM-62` substitution matrix, allowed gaps in the comparison and returned 250 hits.

Output is saved in `data/BLAST_UniProt_250.fasta`

### 2: Generate a multiple sequence alignment (MSA)

We can align sequences using [CLUSTAL Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Set output to `Pearson/FASTA`. The result of the alignment is saved in `data/MSA_clustalomega.fasta`. We should understand if we need to edit rows and columns.

### 3: Build a PSSM model starting from the MSA using BLAST

We can do it with 

```
psiblast -subject data/BLAST_UniProt_250.fasta -in_msa data/MSA_clustalomega.fasta -out_pssm data/profile.pssm
```

Does it make any sense? I'm using a blast result as subject... 

### 4: Build a HMM model starting from the MSA using HMMER

If previous steps are correct it should be something like
```
hmmbuild models/hmm_model.hmm data/MSA_clustalomega.fasta
```
