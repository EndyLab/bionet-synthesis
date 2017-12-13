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

# Import the Gene -> ID dictionary to allow the use of the gene name
for file in glob.glob("./gene_id_dict.json"):
        print("Opening dictionary file: ", file)
        with open(file,"r") as json_file_dict:
            #print(json_file_dict)
            dictionary = json.load(json_file_dict)

# Imports the csv containing all of our data
data = pd.read_csv("./testing/data_testing/fragments1-5.csv")

# Iterates through all of the items in the spreadsheet
for index, row in data.iterrows():

    # Stops it from going through every file in the database
    #counter = counter + 1
    #if counter == 10:
    #    break

    # Takes in the name of the fragment and decides if it has a fragment or not
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
    if 'link' in gene:
        print('link')
        linked_genes.append(gene)
        gene =  gene[:11]

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
        #print(data["location"]["fragments"][name])
        with open(file,'w') as json_file:
                json.dump(data,json_file,indent=2)
    print("cycle")
    print()


print(no_frag)
print(not_in_dict)
