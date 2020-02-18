# Biological Data
 This goal of this project is to investigate the functional and structural properties of the SH2 domain starting from a sample sequence. Our sequence is identified by UniProt ID `P23615(1258-1339)` and can be found as a fasta file `data/sequenceP23615.fasta`.

 ## Requirements
 The project have been developed using Python 3 with the packages cotained in `requirements.txt`.
 It is possible to install the dependencies using `pip`:
 
```
 pip install -r requirements.txt
```

 For the structural alignment the softwar TMalign is needed. Source code and instruction on how to compile it can be found in `code/part2/TMalign`. The executable present in the directory has been compiled on OSX. 

 We used jalview to edit the multiple sequence alignments.

# Part 1: Domain Models

The goal of this first part is to build a PSSM and HMM models representing the asigned domain.

### Step 1
The first step is retrieve homologous sequences from UniProt. To do this we used [Blast](https://www.uniprot.org/blast/) on `UniRef90` and `500` hits. These sequences are save in a fasta file `data/BLAST_uniref90.fasta`.

### Step 2
The retrieved hits have been used to generate a multiple sequence alignment usign [Clustal Omega](https://www.ebi.ac.uk/Tools/msa/clustalo/). The aligment is saved in `data/msa_clustalw.fasta`. 
We then edited the MSA usign JalView. The result is saved in `data/msa_edited.fasta`

### Step 3 
In this step we created and evaluated the models. Models are evaluated against  Pfam annotation of the domain, `PF00017`, in the human organism.  

To do this we need to get from UniProt all the human sequences available in SwissProt containing our domain. This can be done with the query on UniProt

> organism:"Homo sapiens (Human) [9606] AND reviewed:yes  

Results have been downloaded and saved in `data/SwissProt_humans_reference_all.fasta`. The list of true positive, i.e. sequences containing our domain, can be retrieved with a similar query

```
database:(type:pfam pf00017) AND reviewed:yes AND organism:"Homo sapiens (Human) [9606]
```

This has been saved as `data/SwissProt_humans_reference.fasta`, even if it is not necessary.

These files can be used to crete a BLAST database (Blast indexes, phr+pin+psq files) that will be used by HMM-SEARCH and PSI-BLAST to retrieve proteins:

```
makeblastdb -dbtype prot -in data/SwissProt_humans_reference_all.fasta -parse_seqids
```

You can test if it is working correctly by searching a sequence in the database:

```
blastdbcmd -entry "CRK_HUMAN" -db data/SwissProt_humans_reference_all.fasta
```

Creation and evaluation of the models are described on the notebook `code/part1/ModelsEvaluation.ipynb`. Here you can find also a comparison between them. Scripts performing all the steps automatically are in `code/part1/{psiblast.py, hmm.py, jackhmmer.py}`. Other scripts in `code/part1/` contain utility functions.
To run one of them use 

```
python code/part1/psiblast.py
```

Created models are saved in the directory `models` and search output on `results`. 
Parameters for this scripts are contained in the headers. 

# Part 2: Domain family characterization

Once the list of human sequences matching the model is defined, we can look at functional and structural properties.

The first step consists in the creation of four different datasets that will be enriched (gene ontology and disease ontology). 

The first one is called `Original dataset` and is composed of the protein retrieved by our model against human proteins in SwissProt. It can be found in `datasets/original.txt` and it is created by the scripts described in part 1.   

The second is called `Architectures datasets`. This dataset is composed by multiple datasets, where each of them contains proteins of the original dataset with the same combination of Pfam domains. The script used to create this dataset is `code/part2/architectures_dataset.py` and it saves the dataset in `datasets/architectures_datasets.json`. Domains are used as key of the json.

The third is called `PDB_network` and it contains all the sequences of the original dataset with a PDB plus other human proteins found as chains in the same PDB.
The dataset is created with the script `code/part2/pdb_network.py` and saved in `datasets/pdb.csv`.
To find the PDBs the script utilize the file `data/pdb_chain_uniprot.tsv.gz` that can be downloaded from [SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
Together with the sequence and PDB id we saved the position of the PDB in the sequences. This will be used later on to find PDBs covering our domain.

The last dataset is called `STRING network`, composed by all the proteins in the original dataset plus all direct interactors found in the `STRING database`. This dataset is created with the script `code/part2/string_network.py` and save in `datasets/string.txt`.

In the notebooks `code/part2/{GeneOntology.ipynb, DiseaseOntology.ipynb}` the background datasets are identified and the enrichment is performed. 
