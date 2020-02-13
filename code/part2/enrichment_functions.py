from scipy.stats import fisher_exact
import pandas as pd
import copy
import json
import gzip


def parse_gene_ontology(ontology):

    # parents is a dict { term : list_of_parent_terms }
    parents = get_parents(ontology)

    # labels is a dict { term : definition }
    labels = get_labels(ontology)
    # from labels it is possible to retrieve nodes
    # in the following way:  nodes = list(labels.keys())

    # compute the dictionary of ancestors for each node
    ancestors = get_ancestors(
        nodes = list(labels.keys()),
        parents = parents
    )

    # minimum distance from the root for each term
    min_depth = get_min_depth(
        nodes = list(labels.keys()),
        parents = parents 
    )

    return labels, ancestors, min_depth


def get_parents(ontology):
    """
    Retrieve the GO parents, { term : list_of_parent_terms }.
    """
    parents = {} # { term : list_of_parent_terms }
    for edge in ontology["graphs"][0]["edges"]:
        # select only is_a edges
        if edge["pred"] == "is_a":
            parents.setdefault(edge["sub"].split("_")[1], []).append(edge["obj"].split("_")[1])
    return parents


def get_labels(ontology):
    """
    Get labels, i.e a dictionary {term (GO_id): definition}
    """
    labels = {}  # { term : definition }
    for node in ontology["graphs"][0]["nodes"]:
        # exclude obsolete terms
        if "GO_" in node["id"] and "deprecated" not in node["meta"]:
            labels[node["id"].split("_")[1]] = node["lbl"]  
    return labels


def get_ancestors(nodes, parents):
    """
    Compute the list of ancestors for each node
    """
    ancestors = {}  # { term : list_of_ancestor_terms }
    for node in nodes:
        node_ancestors = []
        node_parents = parents.get(node)
        # Loop parent levels until no more parents
        while node_parents:
            node_ancestors.extend(node_parents)
            # Get the parents of current parents (1 level up)
            node_parents = [term for parent in node_parents for term in parents.get(parent, [])]
        ancestors[node] = node_ancestors
    return ancestors

def get_min_depth(nodes, parents):   
    """
    Calculate the minimum depth (distance from the root) of each term
    """
    depth = {}  # { term : min_depth }
    # get list of roots
    roots = set(nodes) - set(parents.keys())
    for node in nodes:
        c = 0  # Depth level
        node_parents = parents.get(node)
        while node_parents:
            c += 1
            if roots.intersection(set(node_parents)):  # break the loop if the root is among parents
                break
            # Get the parents of current parents (1 level up)
            node_parents = [term for parent in node_parents for term in parents.get(parent, [])]
        depth[node] = c
    return depth

def map_protein_to_go(map_file):
    """
    Map proteins ID to GO term
    """
    protein_to_go = {}  # { protein_id : (GO terms) }
    with gzip.open(map_file) as f:
        for acc, annotations in gen_block(f):
            protein_to_go[acc] = annotations
    return protein_to_go

def gen_block(f):
    """
    Parse and split the input.
    The input must be sorted by target name, second column.
    """
    name, old_name = None, None
    chunk = []
    for line in f:
        line = line.decode()
        if line and line[0] != "!":
            _, name, _, _, term, _, ec, _, namespace, protein_name = line.split("\t")[:10]
            term = term[3:]  # remove "GO:" from the term ID
            if name != old_name and old_name:
                # return a set as there can be repetitions, i.e. the same term with different evidence codes
                yield (old_name, set(chunk))  
                chunk = []
            old_name = name
            chunk.append(term)
    # Last line
    if old_name:
        yield (old_name, set(chunk))

def count_ancestors(protein_list, ancestors, protein_to_go):
    """
    Count the ancestors for a list of proteins
    """
    counts = {}

    # the intersection between protein_list and protein_to_go.keys
    # is needed since not all human proteins are annotated
    for protein in set(protein_list).intersection(set(protein_to_go.keys())):
        annotations = protein_to_go[protein]

        terms_ancestors = copy.copy(annotations)  # annotations + ancestor terms
        for term in annotations:  # directly annotated terms
            terms_ancestors.update(set(ancestors.get(term, [])))  # add ancestors
        for term in terms_ancestors:
            counts.setdefault(term, 0)
            counts[term] += 1
        
    return counts


def fisher_test(d_count, bg_count, depth, l):
    """
    Perform a fisher exact test 
    """

    # Init result dict
    results = {}
    
    # Get the tot number of counts
    tot_d = sum(list(d_count.values()))
    tot_bg = sum(list(bg_count.values()))
    
    key_intersection = set(d_count.keys()).intersection(set(bg_count.keys()))
    
    for key in key_intersection:
        # Number of occurrences of the specific GO term in d_count   
        a = d_count[key]
        # Number of occurrences of the specific GO term in bg_count
        b = bg_count[key]
        # Number of GO terms that are different from the specific one in d_count
        not_a = tot_d - a
        # Number of GO terms that are different from the specific one in bg_count
        not_b = tot_bg - b
        
        # Perform Fisher Exact Test
        fisher_results = fisher_exact([[a, b],[not_a, not_b]])
        
        # Create dataframe with results
        results.setdefault(key, {'OddRatio': fisher_results[0], 'p-value': fisher_results[1]})
    
    # Return the DataFrame
    return pd.DataFrame(results).transpose()


def add_depth_description(df, min_depth, labels):
    """
    Add depth and description columns to the dataframe
    """
    labels_list, min_depth_list = [], []
    for key in df.index:
        labels_list.append(labels[key])
        min_depth_list.append(min_depth[key])
    
    df['depth'] = min_depth_list
    df['label'] = labels_list

    return df 