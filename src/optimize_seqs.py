import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import re
import pandas as pd
from config import *
from db_config import *
import ot_functions as ot
from synbiolib import codon
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from numpy.random import choice,randint
from itertools import groupby

## Import data from the database
session,engine = connect_db()
parts_query = "Select * FROM parts"
parts_df = pd.read_sql_query(parts_query,con=engine)

## Initalize starting conditions
cut_sites = [
        ("BbsI", "GAAGAC"),
        ("BtgZI", "GCGATG"),
        ("BsaI", "GGTCTC"),
]

end_codons = {
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
    'G': 'GGC',
    '*': 'TGA'
}

default_table = codon.load_codon_table(taxonomy_id="custom_1", custom=True)
default_table = default_table.reset_index()

test_seq = parts_df.loc[0].original_seq
test_seq = test_seq[:43]+'GAAGAC'+test_seq[43:90]+'GCGATG'+test_seq[90:142]+'AAAAAAAAAAAA'+test_seq[142:]
test_seq = 'ATATAGGAGAAAAAAAAATGGAGATCACA'

# cutsite_check(test_seq)
used_codons = generate_codon_dict(test_seq)
used_codons

seq = test_seq

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

def generate_codon_dict(seq):
    '''Creates a dictionary to store the codons that have been attempted thus far'''
    coding_dna = Seq(seq, generic_dna)
    translation = coding_dna.translate()
    return dict(enumerate(zip(translation,[[seq[n:n+3]] for n in range(0, len(seq), 3)])))

def confirm_translation(original,optimized):
    if Seq(original, generic_dna).translate() == Seq(original, generic_dna).translate():
        print('Matched')
    else:
        print('failed')

def check_cutsites(table,seq,cut_sites,used_codons):
    for enzyme,cut in cut_sites + [(e, reverse_complement(c)) for e, c in cut_sites]:
        start = seq.find(cut)
        print(enzyme,start)
        if start < 0:
            continue
        start = int((start - (start % 3))/3)
        seq = replace_codons(table,seq,range(start,start + (len(cut) // 3 + 1)),used_codons)
        return check_cutsites(table,seq,cut_sites,used_codons)
    print('no cutsites')
    return seq

def check_homopolymer(table,seq,used_codons):
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
    print(max_count)
    if max_count + 1 >= 6:
        start = int((loc - (loc % 3))/3)
        seq = replace_codons(table,seq,range(start,start + (max_count // 3 + 1)),used_codons)
        return check_homopolymer(table,seq,used_codons)
    else:
        print('fixed homopolymer')
        return seq

def replace_codons(table,seq,codon_indexes,used_codons):
    print(codon_indexes)
    if len(codon_indexes) == 0:
        return seq
    section = table[table.AA == used_codons[codon_indexes[0]][0]]
    available = section[~section.Triplet.isin(used_codons[codon_indexes[0]][1])]

    if len(available) == 0:
        print(codon_indexes[0],"can't change")
        return replace_codon(table,seq,codon_indexes[1:],used_codons)
    options = available.Fraction.tolist()
    norm = [float(i)/sum(options) for i in options]
    codons = available.Triplet.tolist()
    new_codon = codons[int(choice(len(options),1,p=norm))]
    used_codons[codon_indexes[0]][1].append(new_codon)
    seq = ''.join((seq[:codon_indexes[0]*3],new_codon,seq[(codon_indexes[0]*3)+3:]))
    return seq

print(test_seq)

seq = check_cutsites(default_table,test_seq,cut_sites,used_codons)
seq = check_homopolymer(default_table,test_seq,used_codons)

print(seq)


size = randint(150,5000)
size -= size % 3
print(type(size))
print(choice(['A','T','C','G']))
test_seq = ''.join(['ATG']+[choice(['A','T','C','G']) for _ in range(size)]+['TGA'])

print(len(test_seq))
print(str(test_seq))
