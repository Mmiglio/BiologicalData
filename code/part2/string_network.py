import urllib.parse
import urllib.request
import pandas as pd
import os

URL_UNIPROT = 'https://www.uniprot.org/uploadlists/'
URL_STRING = 'http://string-db.org/api/tsv/interactorsList?'
PATH_ORIGINAL_DB = '../../datasets/original.txt'

def main():

    # Import original dataset: convert to a unique query string
    original_proteins = ''
    with open(PATH_ORIGINAL_DB) as file:
        for line in file:
            original_proteins += line[:-1] + ' '

    response = make_query(
        URL = URL_UNIPROT,
        params = {
            'from': 'ACC',
            'to': 'STRING_ID',
            'format': 'tab',
            'query': original_proteins[:-1],
            'columns': 'id'
            }
        )

    # Create dataframe with Swissprot - STRING mapping for original dataset
    string_ids = pd.DataFrame(
                list(map(lambda x: x.split(), response.split('\n')[:-1]))[1:], 
                columns = ['id','string'] 
            )

    # retrieve interactors from string
    # with this procedure we are quering STRING for each
    # protein. This may take some time (~1min)
    interactors = []
    for protein in string_ids.string.values:
        # Retrive interactors for the query protein (STRING ID)
        string_request = make_query(
            URL = URL_STRING,
            params = {
                'identifier': protein
            }
        )
        # Header and query protein excluded; last empty string excluded
        string_request = string_request.split('\n')[2:-1]
        interactors.extend(string_request)   

    # Convert string identifiers to uniprot IDs
    map_string_uniprot = make_query(
        URL = URL_UNIPROT,
        params = {
            'from': 'STRING_ID',
            'to': 'ACC',
            'format': 'tab',
            'query': ','.join(set(interactors)),
            'columns': 'id'
        }
    )

    # Create dataframe with Swissprot - STRING mapping for original dataset + interactors
    frames = [string_ids, pd.DataFrame( list(map(
                            lambda x: x.split(), 
                            map_string_uniprot.split('\n')[:-1]
                            ))[1:],
                            columns = ['id','string'] 
                            )
             ]
    string_dataset = pd.concat(frames, sort = False)
    
    # Save ids to file
    print("Saving string network")
    with open('../../datasets/string.txt', 'w') as fout:
        for uniprot_id in set(string_dataset['id'].values):
            fout.write("{}\n".format(uniprot_id))

def make_query(URL, params):
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(URL, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()

    return response.decode('utf-8')

if __name__ == "__main__":
    # set working directory
    os.chdir(os.path.abspath(os.path.dirname(__file__)))
    
    main()