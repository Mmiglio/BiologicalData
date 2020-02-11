import urllib.request
import urllib.parse
import json
import re

ORIGINAL_DATASET = "../../datasets/original.txt"
ARCHITECTURES_DATASET = "../../datasets/architectures_datasets.json"

def main():
    # read sequences from the original dataset
    with open(ORIGINAL_DATASET, "r") as file:
        sequences = []
        for seq in file:
            sequences.append(seq.strip())

    response = make_query(query_sequences=" ".join(sequences))
    rows = response.split('\n') # split string into rows
    header = rows.pop(0).split('\t') #remove header

    architecture_dataset = {}

    for row in rows: 
        if len(row) > 0:
            parsed_row = row.strip().split('\t')
            protein_id = parsed_row[0]
            # split by tabs, pfams are the second column (pfam1;pfam2;)
            pfams = parsed_row[1]
            # [:-1] is needed because the last element is an empty string (last ;)
            set_pfams = set(pfams.split(';')[:-1]) 
            # notes for domains
            names = re.findall(r'/note="([ \w]+)"', parsed_row[2])
            
            # set of domains is the key of the architecture_dataset
            key = frozenset(set_pfams)
            
            # insert protein in the dataset
            architecture_dataset.setdefault(key, []).append(protein_id)

    # key is of the type frozenset(PF07525,PF12610,PF00017)
    # transform it into a string
    for key in list(architecture_dataset.keys()):
        new_key = ",".join(key)
        architecture_dataset[new_key] = architecture_dataset.pop(key)

    with open(ARCHITECTURES_DATASET, 'w') as fp:
        json.dump(architecture_dataset, fp, indent=4)

    print("Dataset saved as {}".format(ARCHITECTURES_DATASET))

def make_query(query_sequences):
    """
    Query uniprot and find pfams for each sequence
    """
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': 'ACC',
        'to': 'ACC',
        'format': 'tab',
        'query': query_sequences,
        'columns': ','.join(['id','database(Pfam)','feature(DOMAIN EXTENT)'])
    }

    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)

    print("Quering UniprotKB")

    with urllib.request.urlopen(req) as f:
        response = f.read()

    print("Done!")

    return response.decode('utf-8')

if __name__ == "__main__":
    main()