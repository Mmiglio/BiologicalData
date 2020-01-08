# BiologicalData

Short Description

## Domain models

### 1: Retrieve homologous proteins in UniProt

Perform a [BLAST search](https://www.uniprot.org/blast/) with our input sequence 

>YYFPFNGRQAEDYLRSKERGEFVIRQSSRGDDHLVITWKLDKDLFQHIDIQELEKENPLALGKVLIVDNQKYNDLDQIIVEY

on UniRef90 (suggested by prof. Piovesan) using `BLOSUM-62` substitution matrix, allowed gaps in the comparison and returned at most 500 hits.

Output is saved in `data/BLAST_uniref90.fasta`

### 2: Generate a multiple sequence alignment (MSA)

We can align sequences using [CLUSTAL Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Remember to set output to `Pearson/FASTA`. 

Result of the alignment is saved in `data/MSA_clustalomega.fasta`. 

### Edit 1
File edited usign Jalview is saved in `data/MSA_clustalomega_edited.fasta`.

### Edit 2
File `data/MSA_clustalomega_soloseq.fasta` obtained by keeping only nucleotides aligned with the starting sequence, 118 (counting gaps) with respect to 81 starting nucleotides (as done in Pfam).

### 3: Build a PSSM model starting from the MSA using BLAST

We can do it with 

```
psiblast -subject data/BLAST_uniref90.fasta -in_msa data/MSA_clustalomega_soloseq.fasta -out_pssm models/profile.pssm
```

Output profile saved in `models/profile.pssm`.

### 4: Build a HMM model starting from the MSA using HMMER

If previous steps are correct it should be something like
```
hmmbuild models/hmm_model.hmm data/MSA_clustalomega_soloseq.fasta
```

HMM model saved in `models/hmm_model.hmm`.

### 5: Evaluate the model against human proteins available in SwissProt

Use metrics such as accuracy, precision, sensitivity, specificity and MCC (Matthews Correlation Coefficient). 
The evaluation is composed by the steps described below.

#### a: Define you ground truth/reference 
The reference database is saved in `data/SwissProt_reference.fasta`. It can be obtained with the following query on UniProt 
>database:(type:pfam pf00017) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"

We can use it to create a database (Blast indexes, phr+pin+psq files) that will be used by HMM-SEARCH and PSI-BLAST to retrieve proteins:
>makeblastdb -dbtype prot -in data/SwissProt_humans_reference.fasta -parse_seqids

You can test if it is working correctly by searching a sequence in the database:
>blastdbcmd -entry "CRK_HUMAN" -db data/SwissProt_humans_reference.fasta

#### b: Find significant hits using HMM-SEARCH and PSI-BLAST respectively for the HMM and PSSM model

We can now search in the database usign PSI-BLAST with the pssm generated in step 3:

>psiblast -in_pssm models/profile.pssm -db data/SwissProt_humans_reference.fasta -outfmt 6 -num_iterations 3 -evalue 0.001 > results/psiblast_search.txt

The output is saved in `results/psiblast_search.txt`. If you want to print the results on screen remove from the command `> results/psiblast_search.txt`.

The search with HMMER is performed with the following command:
> hmmsearch --domtblout results/hmmsearch.hmmer_domtblout models/hmm_model.hmm data/SwissProt_humans_reference.fasta > results/hmmsearch_results.hmmer_align

and it generates two outputs `results/hmmsearch.hmmer_domtblout` and `results/hmmsearch_results.hmmer_align`.

#### c: Evaluate the ability of retrieving proteins with that domain.

#### d: Evaluate the ability of matching the domain position, i.e. the alignment position of the model in the retrieved proteins (Pfam reference position is available in InterPro).



## Notes

* Reference database used in the first part of the project `SwissProt_reference.fasta`. It can be obtained with the following query on UniProt `database:(type:pfam pf00017) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"`
* To visualize the generate PSSM use the following command `psiblast -subject data/BLAST_uniprot_human.fasta -in_msa data/MSA_clustalomega.fasta -out_ascii_pssm test`
