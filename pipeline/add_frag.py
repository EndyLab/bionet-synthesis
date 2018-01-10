## FOR GIT REPOSITORY: Takes in a file containing all of the fragment sequences and then appends the information to the json files

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re
import math

import shutil

counter = 0

no_frag = []
linked_genes = []
not_in_dict = []
genes = []

# Create a dictionary to link the gene name to the corresponding id number
ref = pd.read_csv("./testing/data_testing/10K_CDS.csv")
dictionary = dict(zip(ref['gene_name'], ref['idnum']))

# Imports the csv containing all of the fragment sequences
data = pd.read_csv("./testing/data_testing/fragments1-5.csv")

# Iterates through all of the items in the spreadsheet
for index, row in data.iterrows():

    # Takes in the name and sequence of the fragment and decides if it has a fragment or not
    name = row["Fragment"]
    seq = row["Sequence"]
    if name[-2] != "_":
        no_frag.append(name)
        continue

    # Splits the fragment name into gene name and fragment number
    gene = name[:-2]
    print(gene)
    fragment = name[-1]

    # Checks if it is a linked gene, right now it only takes the first gene
    link_count = 0
    if 'link' in gene:
        print('link')

        # Creates an array with all of the linked genes inside
        genes = gene.split("_link_")

        # Iterates through the each gene to assign the fragment name and sequence
        for linked in genes:
            idnum = dictionary[linked]
            linked_genes.append(idnum)
            id_frag = idnum + "_" + fragment
            for file in glob.glob("../data/{}/{}.json".format(idnum,idnum)):
                print(file)

                # Open the json file
                with open(file,"r") as json_file:
                    data = json.load(json_file)
                data["location"]["fragments"][id_frag] = ""
                data["sequence"]["fragment_sequences"][id_frag] = seq
                with open(file,'w') as json_file:
                    json.dump(data,json_file,indent=2)

                # Uses link count to break out of two loops
                link_count = 1
    if link_count == 1:
        continue

    # Checks if the gene is in the dictionary and if not record it and move on
    if gene not in dictionary:
        print("{} not in dictionary".format(gene))
        not_in_dict.append(gene)
        continue

    # Set the ID# based on the dictionary specification
    print("name", name)
    idnum = dictionary[gene]
    id_frag = idnum + "_" + fragment
    print("id_frag", id_frag)

    # Look up the specific json file associated with the current fragment
    for file in glob.glob("../data/{}/{}.json".format(idnum,idnum)):
        print(file)
        # Open the json file
        with open(file,"r") as json_file:
            data = json.load(json_file)
        data["location"]["fragments"][id_frag] = ""
        data["sequence"]["fragment_sequences"][id_frag] = seq
        with open(file,'w') as json_file:
                json.dump(data,json_file,indent=2)
    print("cycle")
    print()


print(no_frag)
print(not_in_dict)
print()
print("linked genes: ", linked_genes)
