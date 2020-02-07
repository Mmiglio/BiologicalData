import subprocess
from subprocess import DEVNULL

def create_profile(profile_path, data_path, msa_path):
    cmd = """
    rm {0}
    psiblast -subject {1} \
             -in_msa  {2} \
             -out_pssm {0} 
    """.format(
        profile_path,
        data_path,
        msa_path
    )
    
    results = subprocess.run(
        cmd, shell=True, universal_newlines=True, check=True, stdout=DEVNULL)

    if results.stderr:
        print("An error eccured")
    
    print("\nCreated PSSM! Saved as {}\n".format(profile_path))


def create_hmm(model_path, msa_path):
    cmd = """
    rm {0}
    hmmbuild {0} {1}
    """.format(
        model_path,
        msa_path
    )

    results = subprocess.run(
    cmd, shell=True, universal_newlines=True, check=True, stdout=DEVNULL)
    
    if results.stderr:
        print("An error eccured")
    
    print("\nCreate HMM! Saved as {}\n".format(model_path))


def search_psiblast(result_path, pssm_path, db_path):
    cmd = """
    rm {0}
    psiblast -in_pssm {1} -db {2} \
         -outfmt 6 -num_iterations 3 -evalue 0.001 > {0}
    """.format(
        result_path,
        pssm_path,
        db_path
    )

    results = subprocess.run(
    cmd, shell=True, universal_newlines=True, check=True)
    
    if results.stderr:
        print("An error eccured")
    
    print("\nSearch completed! Results saved as {}\n".format(result_path))


def search_hmm(result_path, model_path, db_path):
    cmd = """
    rm {0}
    hmmsearch --domtblout {0} {1} {2}
    """.format(
        result_path,
        model_path,
        db_path
    ) 

    results = subprocess.run(
    cmd, shell=True, universal_newlines=True, check=True, stdout=DEVNULL)
    
    if results.stderr:
        print("An error eccured")

    print("\nSearch completed! Results saved as {}\n".format(result_path))


if __name__ == "__main__":
    #tests 
    create_profile(
        profile_path = "../models/profile.pssm",
        data_path = "../data/BLAST_uniref90.fasta",
        msa_path = "../data/PF00017_seed.fasta"
        )

    create_hmm(
        model_path = "../models/hmm_model.hmm",
        msa_path = "../data/PF00017_seed.fasta"
    )

    search_psiblast(
        result_path = "../results/psiblast_search.txt",
        pssm_path = "../models/profile.pssm",
        db_path = "../data/SwissProt_humans_reference_all.fasta"
    )

    search_hmm(
        result_path = "../results/hmmsearch.hmmer_domtblout",
        model_path = "../models/hmm_model.hmm",
        db_path = "../data/SwissProt_humans_reference_all.fasta"
    )