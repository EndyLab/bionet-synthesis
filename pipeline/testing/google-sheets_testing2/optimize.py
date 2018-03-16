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
import codon

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def optimize_gene(sequence,taxid):
    gene = [{'Gene': 'gene', 'Sequence': sequence}]
    design_genes = pd.DataFrame(gene)

    # Read in genes
    #design_genes = pd.read_csv(sys.stdin)

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

                # REVIEW: Changed from original script
                # If no stop codon is present, add one
                if test == has_stop:
                    eprint("Added stop codon")
                    gene['Sequence'] += "TGA"
                eprint()


    # Force stop codons to TGA
    design_genes['Sequence'] = design_genes['Sequence'].str[:-3] + "TGA"

    ### CLEAN UP THIS CODE:
    codon_table = codon.load_codon_table(taxonomy_id=taxid)
    codon_10plus = codon.codon_table_10plus(codon_table)
    ec_codon_usage_10plus = codon_10plus
    ec_codon_usage = codon_table

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
        changes = 0

        for enzyme, cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
            while cut in seq:
                eprint("{} cuts {} ({})".format(enzyme, name, cut))
                changes += 1
                seq = recode_sequence(seq, cut)

        return seq, changes

    for i, gene in design_genes.iterrows():
        eprint("Optimizing gene {}".format(gene['Gene']))

        sequence = gene['Sequence']
        for j in range(1000): # try 1000 times to get a working sequence
            sequence, changed = remove_cutsites(gene['Gene'], sequence)

            homopolymer_length, homopolymer_pos = max_homopolymer(sequence)
            while homopolymer_length >= 6:
                changed = True
                sequence = recode_sequence(sequence, sequence[homopolymer_pos:homopolymer_pos+homopolymer_length])
                homopolymer_length, homopolymer_pos = max_homopolymer(sequence)

            if not gc_in_range(sequence):
                eprint ("GC out of range")

            if not changed:
                break
        else:
            eprint("Warning: {} could not be optimized".format(gene['Gene']))

        design_genes.ix[i, 'Sequence'] = sequence

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

    # design_genes['Sequence'] = design_genes['Sequence'].str[:-6] + design_genes['Sequence'].str[-6:-3].apply(lambda x: force_codons[ec_codons.ix[x]]) + design_genes['Sequence'].str[-3:]

    # Make sure we didn't introduce an RE site
    # Make sure we didn't introduce any restriction sites
    error = False
    for enzyme, cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
        if design_genes['Sequence'].str.contains(cut).any():
            error = True
            eprint("Cannot normalize final codon without introducing restriction site {} ({})".format(enzyme, cut))

            for i, (g, s) in design_genes.ix[design_genes['Sequence'].str.contains(cut), ['Gene','Sequence']].iterrows():
                eprint(g, s[s.find(cut):])
            eprint()

    if error:
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

    return(design_genes.iloc[0]['Sequence'])
    #design_genes.to_csv(sys.stdout)
