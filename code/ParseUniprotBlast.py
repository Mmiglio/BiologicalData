"""
Parse XML file downloaded from https://www.uniprot.org/blast/
and convert it into fasta file
Use https://jsonformatter.org/xml-viewer to inspect the XML
"""

import xml.etree.ElementTree as ET

tree = ET.parse('../data/BLAST_uniref90.xml')
root = tree.getroot()

NS = {'ebi': 'http://www.ebi.ac.uk/schema' }

hits = root.findall('ebi:SequenceSimilaritySearchResult/ebi:hits/ebi:hit', NS)

output_file = open('../data/BLAST_uniref90_cut.fasta', 'a')

for hit in hits:
    fasta_entry = ""
    
    # add header
    fasta_entry += ">{} {}\n".format(hit.attrib['id'],hit.attrib['description'])
    
    # add sequence
    sequence = hit.find('ebi:alignments/ebi:alignment/ebi:matchSeq',NS).text
    fasta_entry += "\n".join([ x[i:i+60] for i in range(0, len(sequence), 60)])
    
    # new line
    fasta_entry += "\n"
    
    # write to file
    output_file.write(fasta_entry)

output_file.close()