# BiologicalData

Short Description

## Domain models

### 1: Retrieve homologous proteins in UniProt

Perform a [BLAST search](https://www.uniprot.org/blast/) with our input sequence 

>YYFPFNGRQAEDYLRSKERGEFVIRQSSRGDDHLVITWKLDKDLFQHIDIQELEKENPLALGKVLIVDNQKYNDLDQIIVEY

on UniProtKB/humans using `BLOSUM-62` substitution matrix, allowed gaps in the comparison and returned at most 250 hits.

Output is saved in `data/BLAST_uniprot_human.fasta`

### 2: Generate a multiple sequence alignment (MSA)

We can align sequences using [CLUSTAL Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). Remember to set output to `Pearson/FASTA`. 

Result of the alignment is saved in `data/MSA_clustalomega.fasta`. 

We should understand if we need to edit rows and columns.

### 3: Build a PSSM model starting from the MSA using BLAST

We can do it with 

```
psiblast -subject data/BLAST_uniprot_human.fasta -in_msa data/MSA_clustalomega.fasta -out_pssm data/profile.pssm
```

Output profile saved in `data/profile.pssm`.

Does it make any sense? I'm using a blast result as subject... 

### 4: Build a HMM model starting from the MSA using HMMER

If previous steps are correct it should be something like
```
hmmbuild models/hmm_model.hmm data/MSA_clustalomega.fasta
```

HMM model saved in `models/hmm_model.hmm`.

### 5: Evaluate the model against human proteins available in SwissProt

The reference database is saved in `data/SwissProt_reference.fasta`. It can be obtained with the following query on UniProt 
>database:(type:pfam pf00017) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"



## Notes

* Reference database used in the first part of the project `SwissProt_reference.fasta`. It can be obtained with the following query on UniProt `database:(type:pfam pf00017) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]"`