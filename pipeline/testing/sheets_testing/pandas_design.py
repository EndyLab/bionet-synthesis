import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from difflib import ndiff
import re


import json
import io
import sys
import math

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Read in genes
design_genes = pd.read_csv(sys.stdin)


# Config fragments
BsaI_prefix = "GGAGGTCTCA".lower()
BsaI_suffix = "CGAGACCGCT".lower()


# Add preix / suffix to each gene
for index, row in design_genes.iterrows():
    part_type = (row['Part_Type'])
    with open("overhang-map.json", "r") as json_file:
        data_type = json.load(json_file)
        case = (data_type[part_type]["type"]["special_case"])
        prefix = (data_type[part_type]["prefix"])
        suffix = (data_type[part_type]["suffix"])

        if case == "true":
            ## Work on this part
            design_genes['Sequence'] = design_genes['Sequence'] 
        else:
            design_genes['Sequence'] = BsaI_prefix + prefix + design_genes['Sequence'] + suffix + BsaI_suffix

print (design_genes)

# What size are our things?
eprint("Gene Number")
for i in np.arange(1,5)*1500:
    eprint(">", i, (design_genes['Sequence'].str.len() >= i).sum())
eprint()


small_prefix = "GCGATGcatgcttgca".lower()
small_suffix = "cctctgaataCATCGC".lower()
frag_prefix = "GAAGACAT".lower()
frag_suffix = "TCGTCTTC".lower()


def fragment_gene(seq):
    max_len = 1500
    num_frags = len(seq) // max_len + 1
    frag_len = len(seq) // num_frags

    frags = []
    if len(seq) < 300:
        frag = small_prefix + seq + small_suffix
        frags.append(frag)
    else:
        for i in range(num_frags):
            frag = seq[max(0, i * frag_len - 2):min((i+1) * frag_len + 2,len(seq))]
            frag = frag_prefix + frag + frag_suffix
            frags.append(frag)

    return frags

fragments = pd.DataFrame({
            'Gene': [],
            'Fragment': [],
            'Sequence': []
        })

# Fragment big (>300 bp) genes
small_genes = pd.DataFrame({
            'Gene': [],
            'Fragment': [],
            'Sequence': []
        })
for i, gene in design_genes.iterrows():
    # Skip small genes for assembly into multi-gene constructs
    if len(gene['Sequence']) < 300:
        frags = fragment_gene(gene['Sequence'])
        small_genes = small_genes.append(
                pd.DataFrame({
            'Gene': [gene['Gene']] * len(frags),
            'Fragment': range(1,len(frags)+1),
            'Sequence': frags
            }), ignore_index=True)
        continue

    frags = fragment_gene(gene['Sequence'])
    fragments = fragments.append(
        pd.DataFrame({
            'Gene': [gene['Gene']] * len(frags),
            'Fragment': range(1,len(frags)+1),
            'Sequence': frags
        }), ignore_index=True)


# Add fragment to gene name
for index, row in fragments.iterrows():
    name_with_fragment= (row['Gene']) + "_" + str(int((row['Fragment'])))
    fragments.at[index, 'Gene'] = name_with_fragment

# Sort dataframe by length
new_index = fragments.Sequence.str.len().sort_values().index
fragments = (fragments.reindex(new_index))
fragments = fragments.reset_index(drop=True)

for index, row in small_genes.iterrows():
    sequence = (row['Sequence'])
    name = (fragments.get_value(index, 'Gene')) + "_link_" + (row['Gene'])
    sequence = (fragments.get_value(index, 'Sequence')) + (row['Gene']) 
    fragments['Gene'].replace(
            to_replace=(fragments.get_value(index, 'Gene')),
            value= name,
            inplace=True
            )
    fragments['Sequence'].replace(
            to_replace=(fragments.get_value(index, 'Sequence')),
            value= sequence,
            inplace=True
            )



eprint(small_genes)
eprint(fragments)


# Print some stats
eprint("We will synthesize {} bp of DNA in {} fragments".format(fragments['Sequence'].str.len().sum(), fragments.shape[0]))


# Output DNA in Twist order format
twist_dna = pd.DataFrame({
        'Sequence ID': fragments['Gene'], 
        'Sequence Description': [''] * fragments.shape[0],
        'Insert Sequence': fragments['Sequence'],
        'Output Format': 'linear',
        'Vessel Type': '96 well'
        }, columns=['Sequence ID', 'Sequence Description', 'Insert Sequence', 'Output Format', 'Vessel Type'])
twist_dna.to_csv(sys.stdout)

