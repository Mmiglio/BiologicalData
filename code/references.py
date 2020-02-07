from Bio import SeqIO
import json

def getHumanSH2(fasta_path):
    """
    Return a set of sequences {'P16885', ...}
    """
    human_sh2 = SeqIO.parse(fasta_path,'fasta')
    list_human_sh2 = []
    for sequence in human_sh2:
        name = sequence.id # name is in the form sp|P46108|CRK_HUMAN
        list_human_sh2.append(name.split('|')[1])

    set_human_sh2 = set(list_human_sh2)
    print("There are {} human proteins containing SH2 domain in SwissProt".format(len(list_human_sh2)))
    return set_human_sh2

def countSequences(fasta_path):
    seq = SeqIO.parse(fasta_path, 'fasta')
    list_seq = []
    for sequence in seq:
        name = sequence.id # name is in the form sp|P46108|CRK_HUMAN
        list_seq.append(name.split('|')[1])
    print("There are {} human proteins in SwissProt".format(len(list_seq)))
    return len(list_seq)

def getPositionReference(json_path):
    """
    return dict of the type
    {'P16885' :{'length': 1265,
               'positions': [{'start': 532, 'end': 617}, {'start': 646, 'end': 720}]}}
    """
    with open(json_path, "r") as file:
        interpro_data = json.load(file)       
    reference_sh2_positions = {}
    for el in interpro_data:
        name = el['metadata']['accession']
        length = el['metadata']['length']
        if name in reference_sh2_positions:
            print("{} already in reference dict!".format(name))
            continue
        else:
            reference_sh2_positions[name] = {'length':length, 'positions':[]}
        entries = el['entries']
        for entry in entries:
            if entry['accession']=='PF00017':
                domain_locations = entry['entry_protein_locations']
                for location in domain_locations:
                    if len(location['fragments'])>1:
                        print("len fragments > 1")
                    position_dict = {
                        'start': location['fragments'][0]['start'],
                        'end': location['fragments'][0]['end']
                    }
                    reference_sh2_positions[name]['positions'].append(position_dict)

    return reference_sh2_positions



    

