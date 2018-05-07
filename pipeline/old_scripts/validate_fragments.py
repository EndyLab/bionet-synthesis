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

validation = []
twist = pd.read_csv(sys.stdin)
fragments = pd.DataFrame({
    'Gene': twist['Sequence ID'].str.rsplit('_', n=1).str[0],
    'Fragment': twist['Sequence ID'].str.rsplit('_', n=1).str[1],
    'Sequence': twist['Insert Sequence']
})

# Process single gene fragments
for i, gene in fragments[~fragments['Gene'].str.contains('_link_')].groupby('Gene'):
    frags = []
    for seq in gene['Sequence'].str.upper():
        # Virtual digest using BbsI
        seq = seq[seq.find('GAAGAC')+8:]
        seq = seq[:seq.find('GTCTTC')-2]
        frags.append(seq)

    # Now we have a list of 'digested' fragments; assemble them based on overlaps.
    while len(frags) > 1:
        frag = frags.pop(0)
        for other in frags:
            if other[-4:] == frag[:4]:
                frag = other[:-4] + frag
                frags.remove(other)
                break
        frags.append(frag)

    # Add assembled full sequence to the validation database
    # Remove the first and last nucleotides to get the CDS instead of MoClo adapters
    validation.append({'Gene': gene['Gene'].iloc[0], 'Sequence': frags[0][1:-1]})

# Process multigene fragments
for i, gene in fragments[fragments['Gene'].str.contains('_link_')].iterrows():
    names = gene['Gene'].split('_link_')

    # Virtual digest using BbsI
    seq = gene['Sequence'].upper()
    seq = seq[seq.find('GAAGAC')+8:]
    seq = seq[:seq.find('GTCTTC')-2]
    validation.append({'Gene': names[0], 'Sequence': seq[1:-1]})

    # We have a multi-gene fragment, digest with BtgZI
    seq = gene['Sequence'].upper()
    seq = seq[seq.find('GCGATG')+16:]
    seq = seq[:seq.find('CATCGC')-10]
    validation.append({'Gene': names[1], 'Sequence': seq[1:-1]})

validation_dna = pd.DataFrame(validation)
validation_dna.to_csv(sys.stdout)
