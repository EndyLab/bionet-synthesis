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

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Read in genes
design_genes = pd.read_csv(sys.stdin)

# Basic validation tests
def is_seq(seq):
    return seq.replace(r'\s +', '').strip() != ""

def is_triplet(seq):
    return len(seq) % 3 == 0

def translate(seq):
    return str(Seq(seq, generic_dna).translate(table=11))

def dna_matches_aa(seq, aa):
    prot = translate(seq)
    if prot[-1] == "*":
        prot = prot[:-1]

    return prot == aa

def has_start(seq):
    return seq[:3] == "ATG"

def has_stop(seq):
    return translate(seq)[-1] == "*"

def has_no_internal_stops(seq):
    return not "*" in translate(seq)[:-1]

basic_tests = [is_seq, is_triplet, has_start, has_stop, has_no_internal_stops]

for i, gene in design_genes.iterrows():
    for test in basic_tests:
        if not test(gene['Sequence']):
            eprint("{} failed on {}".format(gene['Gene'], test.__name__))
            eprint(gene['Sequence'])
            eprint(translate(gene['Sequence']))
            eprint()


# Force stop codons to TGA
design_genes['Sequence'] = design_genes['Sequence'].str[:-3] + "TGA"

# Optimization

# Rough and ready codon optimisation code
ec_codon_usage_data = '''
Triplet AA Fraction Frequency Number
TTT F 0.58 22.1 ( 80995)  TCT S 0.17 10.4 ( 38027)  TAT Y 0.59 17.5 ( 63937)  TGT C 0.46  5.2 ( 19138)
TTC F 0.42 16.0 ( 58774)  TCC S 0.15  9.1 ( 33430)  TAC Y 0.41 12.2 ( 44631)  TGC C 0.54  6.1 ( 22188)
TTA L 0.14 14.3 ( 52382)  TCA S 0.14  8.9 ( 32715)  TAA * 0.61  2.0 (  7356)  TGA * 0.30  1.0 (  3623)
TTG L 0.13 13.0 ( 47500)  TCG S 0.14  8.5 ( 31146)  TAG * 0.09  0.3 (   989)  TGG W 1.00 13.9 ( 50991)
CTT L 0.12 11.9 ( 43449)  CCT P 0.18  7.5 ( 27340)  CAT H 0.57 12.5 ( 45879)  CGT R 0.36 20.0 ( 73197)
CTC L 0.10 10.2 ( 37347)  CCC P 0.13  5.4 ( 19666)  CAC H 0.43  9.3 ( 34078)  CGC R 0.36 19.7 ( 72212)
CTA L 0.04  4.2 ( 15409)  CCA P 0.20  8.6 ( 31534)  CAA Q 0.34 14.6 ( 53394)  CGA R 0.07  3.8 ( 13844)
CTG L 0.47 48.4 (177210)  CCG P 0.49 20.9 ( 76644)  CAG Q 0.66 28.4 (104171)  CGG R 0.11  5.9 ( 21552)
ATT I 0.49 29.8 (109072)  ACT T 0.19 10.3 ( 37842)  AAT N 0.49 20.6 ( 75436)  AGT S 0.16  9.9 ( 36097)
ATC I 0.39 23.7 ( 86796)  ACC T 0.40 22.0 ( 80547)  AAC N 0.51 21.4 ( 78443)  AGC S 0.25 15.2 ( 55551)
ATA I 0.11  6.8 ( 24984)  ACA T 0.17  9.3 ( 33910)  AAA K 0.74 35.3 (129137)  AGA R 0.07  3.6 ( 13152)
ATG M 1.00 26.4 ( 96695)  ACG T 0.25 13.7 ( 50269)  AAG K 0.26 12.4 ( 45459)  AGG R 0.04  2.1 (  7607)
GTT V 0.28 19.8 ( 72584)  GCT A 0.18 17.1 ( 62479)  GAT D 0.63 32.7 (119939)  GGT G 0.35 25.5 ( 93325)
GTC V 0.20 14.3 ( 52439)  GCC A 0.26 24.2 ( 88721)  GAC D 0.37 19.2 ( 70394)  GGC G 0.37 27.1 ( 99390)
GTA V 0.17 11.6 ( 42420)  GCA A 0.23 21.2 ( 77547)  GAA E 0.68 39.1 (143353)  GGA G 0.13  9.5 ( 34799)
GTG V 0.35 24.4 ( 89265)  GCG A 0.33 30.1 (110308)  GAG E 0.32 18.7 ( 68609)  GGG G 0.15 11.3 ( 41277)
'''
# From http://www.genscript.com/tools/codon-frequency-table

ec_codon_usage_data = re.sub(r'\( *', '', ec_codon_usage_data)
ec_codon_usage_data = re.sub(r'\) *', '\n', ec_codon_usage_data)

ec_codon_usage = pd.read_table(io.StringIO(ec_codon_usage_data), sep=r' +')
ec_codon_usage.set_index(['AA', 'Triplet'], inplace=True)
ec_codon_usage.sort_index(inplace=True)

ec_codon_usage_10plus = ec_codon_usage.ix[ec_codon_usage.Fraction >= 0.1]
ec_codon_usage_10plus = ec_codon_usage_10plus.groupby(level=0).transform(lambda x: x / x.sum())

ec_codons = ec_codon_usage.reset_index(level=0).ix[:,'AA']

def pick_codon(aa):
    return ec_codon_usage_10plus.ix[x].iloc[((ec_codon_usage_10plus.ix[x].Fraction).cumsum() < np.random.rand()).sum()].name

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

def recode_sequence(seq, rep):
    pos = seq.find(rep)

    if pos < 0:
        return seq

    pos -= pos % 3

    for i in range(pos, pos + (len(rep) // 3 + 1) * 3, 3):
        codon = seq[i:i+3]
        choices = ec_codon_usage_10plus.ix[ec_codons.ix[codon]]
        choices = choices[choices.index != codon]

        if choices.shape[0] > 0:
            newcodon = choices.iloc[(choices.Fraction.cumsum() / choices.Fraction.cumsum().max() < np.random.rand()).sum()].name # Stochastically allocate codon
#             newcodon = choices[choices.index != codon].Fraction.idxmax() # Deterministically allocate codons by using the most frequence one
            break

    eprint("{} -> {}".format(codon, newcodon))

    return seq[:i] + newcodon + seq[i+3:]

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)

def gc_in_range(seq):
    return (gc_content(seq) > 0.3) & (gc_content(seq) <= 0.65)

def size_in_range(seq):
    return (len(seq) >= 300) & (len(seq) <= 1800)

def max_homopolymer(seq):
    prev = ''
    count = 0
    max_count = 0
    loc = None

    for i, char in enumerate(seq):
        if char == prev:
            count += 1
        else:
            count = 0

        if count > max_count:
            max_count = count
            loc = i - count

        prev = char

    return max_count, loc

def homopolymer_in_range(seq):
    count, loc =  max_homopolymer(seq)
    return max_homopolymer(seq) < 6

cut_sites = [
    ("BfuAI", "ACCTGC"),
    ("AarI", "CACCTGC"),
    ("BtgZI", "GCGATG"),
    ("BbsI", "GAAGAC"),
    ("BsmBI", "CGTCTC"),
    ("SapI", "GCTCTTC"),
    ("BsaI", "GGTCTC")]

def remove_cutsites(name, seq):
    for enzyme, cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
        while cut in seq:
            eprint("{} cuts {} ({})".format(enzyme, name, cut))
            seq = recode_sequence(seq, cut)

    return seq

for i, gene in design_genes.iterrows():
    eprint("Optimizing gene {}".format(gene['Gene']))
    for j in range(1000): # try 1000 times to get a working sequence
        changed = False

        oldgene = gene['Sequence']
        gene['Sequence'] = remove_cutsites(gene['Gene'], gene['Sequence'])
        if gene['Sequence'] != oldgene:
            changed = True

        homopolymer_length, homopolymer_pos = max_homopolymer(gene['Sequence'])
        while homopolymer_length >= 6:
            changed = True
            gene['Sequence'] = recode_sequence(gene['Sequence'], gene['Sequence'][homopolymer_pos:homopolymer_pos+homopolymer_length])
            homopolymer_length, homopolymer_pos = max_homopolymer(gene['Sequence'])

        if not gc_in_range(gene['Sequence']):
            eprint ("GC out of range")

        if not changed:
            break
    else:
        eprint("Warning: {} could not be optimized".format(gene['Gene']))

    eprint("")

# Force final codons to our standard set
force_codons = {
    'M': 'ATG',
    'W': 'TGG',
    'F': 'TTT',
    'L': 'CTG',
    'I': 'ATT',
    'V': 'GTG',
    'S': 'TCC',
    'P': 'CCA',
    'T': 'ACC',
    'A': 'GCC',
    'Y': 'TAC',
    'H': 'CAT',
    'Q': 'CAG',
    'N': 'AAC',
    'K': 'AAG',
    'D': 'GAT',
    'E': 'GAG',
    'C': 'TGC',
    'R': 'CGC',
    'G': 'GGC'
}

design_genes['Sequence'] = design_genes['Sequence'].str[:-6] + design_genes['Sequence'].str[-6:-3].apply(lambda x: force_codons[ec_codons.ix[x]]) + design_genes['Sequence'].str[-3:]

# Make sure we didn't introduce an RE site
# Make sure we didn't introduce any restriction sites
for enzyme, cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
    if design_genes['Sequence'].str.contains(cut).any():
        eprint("Cannot normalize final codon without introducing restriction site {} ({})".format(enzyme, cut))
        eprint(design_genes.ix[design_genes['Sequence'].str.contains(cut), ['Gene','Sequence']])
        eprint()
        sys.exit(1)

#         seq = gene['Sequence']
#         while gc_content(seq) < 0.35:
#             # Too low, let's bump it up

#             # Pick a random codon
#             seq = gene['Sequence']
#             index = np.random.randint(len(seq))
#             index -= index % 3
#             codon = seq[index:index+3]

#             # Recode it for higher GC by replacing with the highest GC codon we can use
#             # Yup, this is a seriously inefficient way to do this
#             codons = ec_codon_usage.ix[ec_codons[codon]].index.values.tolist()
#             codons.sort(key=gc_content)
#             newcodon = codons[-1]
#             seq = seq[:index] + newcodon + seq[index+3:]

#            eprint('{} -> {}'.format(codon, newcodon))
#             break

#    eprint("")

design_genes.to_csv(sys.stdout)
