import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math
import datetime
from datetime import datetime

from Bio import pairwise2
from Bio import AlignIO
from Bio import Align
from Bio import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import SeqIO

gene_names = []
wells = []
counter = 0

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv("../../../pipeline/testing/data_testing/10K_CDS.csv")
dictionary = dict(zip(data['gene_name'], data['idnum']))

dist_map = pd.read_csv('MMSYN1_plate1.csv')

backbone_seq = SeqIO.read('../../../pipeline/popen_v1-1_backbone.fasta', 'fasta')
backbone_seq = str(backbone_seq.seq)

old_insert = "ATGCGGTCTTCCGCATCGCCTGGCACGACAGGTTTCCCGACTGGAAAGCGGGCAGTGAGCGCAACGCAATTAATGTGAGTTAGCTCACTCATTAGGCACCCCAGGCTTTACACTTTATGCTTCCGGCTCGTATGTTGTGTGGAATTGTGAGCGGATAACAATTTCACACATACTAGAGAAAGAGGAGAAATACTAGATGGCTTCCTCCGAAGATGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGTTCCGTTAACGGTCACGAGTTCGAAATCGAAGGTGAAGGTGAAGGTCGTCCGTACGAAGGTACCCAGACCGCTAAACTGAAAGTTACCAAAGGTGGTCCGCTGCCGTTCGCTTGGGACATCCTGTCCCCGCAGTTCCAGTACGGTTCCAAAGCTTACGTTAAACACCCGGCTGACATCCCGGACTACCTGAAACTGTCCTTCCCGGAAGGTTTCAAATGGGAACGTGTTATGAACTTCGAGGACGGTGGTGTTGTTACCGTTACCCAGGACTCCTCCCTGCAAGACGGTGAGTTCATCTACAAAGTTAAACTGCGTGGTACCAACTTCCCGTCCGACGGTCCGGTTATGCAGAAAAAAACCATGGGTTGGGAAGCTTCCACCGAACGTATGTACCCGGAGGACGGTGCTCTGAAAGGTGAAATCAAAATGCGTCTGAAACTGAAAGACGGTGGTCACTACGACGCTGAAGTTAAAACCACCTACATGGCTAAAAAACCGGTTCAGTTACCGGGTGCTTACAAAACCGACATCAAACTGGACATCACCTCCCACAACGAGGACTACACCATCGTTGAACAGTACGAACGTGCTGAAGGTCGTCACTCCACCGGTGCTTAAGCGATGTTGAAGACCATGA"

for index,row in dist_map.iterrows():

    gene_name = row['Gene Name']
    if gene_name in dictionary.keys():
        id_num = dictionary[gene_name]
    else:
        id_num = "Empty"
        print('Empty')
        continue
    # Uses the ID# to open the correct .json file and then adds the location
    for file in glob.glob("../../../data/{}/{}.json".format(id_num,id_num)):
        print(file)
        # Open the json file
        with open(file,"r") as json_file:
            data = json.load(json_file)
    seq = data["sequence"]["optimized_sequence"]
    print(id_num,seq)
    new_seq = backbone_seq.replace(old_insert,seq)
    print(new_seq)

    # Generate the fasta file with the gene sequence
    fasta = open("./plasmid_maps/{}_popen_v1-1.fasta".format(gene_name),"w+")
    fasta.write(">{}_popen_v1-1".format(gene_name))
    fasta.write("\n")
    fasta.write(new_seq)
    fasta.close()



#
