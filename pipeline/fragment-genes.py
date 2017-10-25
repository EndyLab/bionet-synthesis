import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from difflib import ndiff
import re

import io
import sys
import math

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Read in genes
design_genes = pd.read_csv(sys.stdin)

# Add assembly flanking sequences

# Prefix is BbsI cut site, followed by one arbitrary nucleotide
# Suffix is the same in reverse complement
# Our synthesized DNA will have 25bp flanking sequences added by Twist, so we don't need to
# extend out further.

prefix = "GAAGACAT".lower()
suffix = "TCGTCTTC".lower()

# Add an A to each end of each sequence to create the standard MoClo sticky ends
# Note that we assume the stop codons on incoming genes are TGA
design_genes['Sequence'] = "a" + design_genes['Sequence'] + "a"

# What size are our things?
eprint("Gene Number")
for i in np.arange(1,5)*1780:
    eprint(i, (design_genes['Sequence'].str.len() >= i).sum
eprint()

# Break genes into fragments
# We're going to divide the genes that are bigger than 1800 bp in half, and those bigger than 3600 into thirds. We'll store the new sequences in a dataframe of actual fragments to synthesize. Nice future project: generalize this for arbitrary sizes.
# Note: We use sizes slightly smaller than 1800/3600 to compensate for the addition of 2 bp of overlap on each fragment and the size of the prefix/suffix.

def fragment_gene(seq):
    max_len = 1800-2-len(prefix)-len(suffix)
    num_frags = len(seq) // max_len + 1
    frag_len = len(seq) // num_frags

    # eprint(max_len, num_frags, frag_len)

    frags = []
    for i in range(num_frags):
        frag = seq[max(0, i * frag_len - 2):min((i+1) * frag_len + 2,len(seq))]
        frag = prefix + frag + suffix
        frags.append(frag)

    return frags

fragments = pd.DataFrame({
            'Gene': [],
            'Fragment': [],
            'Sequence': []
        })

# Fragment big (>300 bp) genes
small_genes = []
for i, gene in design_genes.iterrows():
    # Skip small genes for assembly into multi-gene constructs
    if len(gene['Sequence']) < 300:
        small_genes.append(gene)
        continue

    frags = fragment_gene(gene['Sequence'])
    fragments = fragments.append(
        pd.DataFrame({
            'Gene': [gene['Gene']] * len(frags),
            'Fragment': range(1,len(frags)+1),
            'Sequence': frags
        }), ignore_index=True)

# Link small (<300) bp genes
small_prefix = "GAAGACAT"
small_linker = "TCGTCTTCGCGATGcatgcttgca"
small_suffix = "cctctgaataCATCGC"

# Pair biggest with smallest, etc
small_genes.sort(key=lambda x: len(x['Sequence']))
small_pairs = zip(small_genes[:math.ceil(len(small_genes)/2)], small_genes[:math.floor(len(small_genes)/2)-1:-1])

eprint('Small Pairs')
for a,b in small_pairs:
    eprint(a['Gene'], b['Gene'])
    construct = small_prefix + a['Sequence'] + small_linker + b['Sequence'] + small_suffix

    fragments = fragments.append(
            pd.DataFrame({
            'Gene': [a['Gene'] + "_link_" + b['Gene']],
            'Fragment': [1],
            'Sequence': [construct]
        }), ignore_index=True)
eprint()

# Print some stats
eprint("We will synthesize {} bp of DNA in {} fragments".format(fragments['Sequence'].str.len().sum(), fragments.shape[0]))

# Output DNA in Twist order format
twist_dna = pd.DataFrame({
        'Sequence ID': fragments['Gene'] + "_" + fragments['Fragment'].astype(np.int).astype(np.str),
        'Sequence Description': [''] * fragments.shape[0],
        'Insert Sequence': fragments['Sequence'],
        'Output Format': 'linear',
        'Vessel Type': '96 well'
        }, columns=['Sequence ID', 'Sequence Description', 'Insert Sequence', 'Output Format', 'Vessel Type'])

twist_dna.to_csv(sys.stdout)
