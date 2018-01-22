## FOR GIT REPOSITORY -- Sorting Sequencing files
## Read in a sequencing data and copy the files to the gene directory that they are associated with

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re

import shutil

from Bio import pairwise2
from Bio import AlignIO
from Bio import Align
from Bio import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_dna

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv("./testing/data_testing/10K_CDS.csv")
dictionary = dict(zip(data['gene_name'], data['idnum']))

#State the primers
forward_primer = "M13-Forward---20-"
reverse_primer = "M13-Reverse"

# Track genes that are not in the dictionary
not_in_dict = []
count = 0
already_sorted = False

# Iterate through every sequencing file in the sequencing_files directory
for seqfile in glob.glob("../sequencing_files/*/*.ab1"):

    count += 1
    if count == 10:
        break

    # If statement to accomodate for the occasional "-" in front of the name
    if seqfile[55] == "-":
        initials, order_number, plate_number, well_number, hyphen, sample_name, primer_name, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_(-)([A-Za-z0-9_-]+)_([A-Za-z0-9-]+)_([A-H][0-9]{2}).ab1',
            seqfile).groups()
        revfile = "{}/{}_{}-{}{}_{}{}_{}_{}.ab1".format(os.path.dirname(seqfile), initials, order_number, (int(plate_number)+1), well_number, hyphen, sample_name, reverse_primer, well_address)
    else:
        initials, order_number, plate_number, well_number, sample_name, primer_name, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_([A-Za-z0-9_-]+)_([A-Za-z0-9-]+)_([A-H][0-9]{2}).ab1',
            seqfile).groups()
        revfile = "{}/{}_{}-{}{}_{}_{}_{}.ab1".format(os.path.dirname(seqfile), initials, order_number, (int(plate_number)+1), well_number, sample_name, reverse_primer, well_address)

    print("sample name", sample_name)

    # ONLY HANDLES FIRST SEQUENCE: Checks if it is a linked gene and if it is only take the first one
    if 'link' in sample_name:
        print('link')
        gene, gene2 = sample_name.split('_link_')

    # Checks if the name is already the ID# and if not it converts it
    if sample_name[:6] != "BBF10K":
        if sample_name[:-2] in dictionary.keys():
            idnum = dictionary[sample_name[:-2]]
        else:
            print("{} not in dictionary".format(sample_name[:-2]))
            not_in_dict.append(sample_name)
    else:
        print("{} already in proper format".format(sample_name))
        idnum = sample_name
    print("idnum", idnum)

    # Copies the sequencing file to the corresponding directory
    shutil.copy(seqfile, '../data/{}'.format(idnum))

# Prints all of the sequences that could not be sorted
print("Not in dictionary or database:")
print(not_in_dict)
