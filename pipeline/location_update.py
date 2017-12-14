## FOR GIT REPOSITORY -- ADD LOCATIONS AND STATE BUILDABILITY
## Read in a plate_map.csv, add locations, and assesses if it is buildable by iterating through the list

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

from Bio import pairwise2
from Bio import AlignIO
from Bio import Align
from Bio import Seq
from Bio import SeqRecord
from Bio.Alphabet import generic_dna

linked_genes = []
not_in_dict = []
names = []

# Import the Gene -> ID dictionary to allow the use of the gene name
for file in glob.glob("./gene_id_dict.json"):
        print(file)
        with open(file,"r") as json_file_dict:
            print(json_file_dict)
            dictionary = json.load(json_file_dict)

# Iterate through all of the plate maps in the directory
for plate in glob.glob("../plate_maps/*.csv"):
    plates = pd.read_csv(plate)

    # Modify the customer line item ids to extract the gene name/ID from them
    plates['customer_line_item_id'] = plates['customer_line_item_id'].str.strip()
    ids = plates['customer_line_item_id'].str.rsplit('_', n=1)
    plates['Gene'] = ids.str[0]
    plates['Fragment'] = ids.str[1]
    plates['Location'] = plates['Plate'] + "_" + plates["Well"]

    # Iterate through all of the items listed in the plate map
    for index, row in plates.iterrows():
        name = row['Gene']
        print("name", name)

        # Currently only dealing with the BbsI flanked fragments but tracks the others
        link_count = 0
        if 'link' in name:
            print('link')
            names = name.split("_link_")
            for linked in names:
                gene = dictionary[linked]
                linked_genes.append(gene)
                fragment_num = str(row['Fragment'])
                fragment = gene + "_" + fragment_num
                # Search the data base for the corresponding json file
                for file in glob.glob("../data/{}/{}.json".format(gene,gene)):
                    # Acts as a counter to determine if all of the fragments have been received
                    num_locations = 0
                    print("file", file)

                    # Open the json file
                    with open(file,"r") as json_file:
                        data = json.load(json_file)

                    # Iterate through all of the fragments
                    for frag in data['location']['fragments']:

                        # Searches through the fragments to find the one matching the line item in the plate map
                        if frag == fragment:
                            print("match", frag)
                            print("old location", data['location']['fragments'][frag])

                            # Updates the data with the new location
                            data['location']['fragments'][frag] = row['Location']
                            print("new location", data['location']['fragments'][frag])
                        else:
                            print("no match")

                        # Adds to the counter as it iterates through the fragments
                        if data['location']['fragments'][frag] != "":
                            print("has location")
                            num_locations = num_locations + 1
                        else:
                            print("has no location")

                    # Determines if the gene is build ready
                    if num_locations == len(data['location']['fragments']):
                        print("buildable")
                        data['status']["build_ready"] = "TRUE"
                    else:
                        print("not buildable")

                    # Writes a new json with the modified data
                    with open(file,'w') as json_file:
                            json.dump(data,json_file,indent=2)
                    print("cycle done")
                    print("")
                link_count = 1
        if link_count == 1:
            continue

        # FOR INITIAL PLATES: checks if it is in the database
        if name not in dictionary:
            not_in_dict.append(name)
            continue

        # FOR INITIAL PLATES: Remake the fragment name with the correct gene name
        gene = dictionary[name]
        fragment_num = str(row['Fragment'])
        fragment = gene + "_" + fragment_num
        print("fragment", fragment)

        # Search the data base for the corresponding json file
        for file in glob.glob("../data/{}/{}.json".format(gene,gene)):
            # Acts as a counter to determine if all of the fragments have been received
            num_locations = 0
            print("file", file)

            # Open the json file
            with open(file,"r") as json_file:
                data = json.load(json_file)

            # Iterate through all of the fragments
            for frag in data['location']['fragments']:

                # Searches through the fragments to find the one matching the line item in the plate map
                if frag == fragment:
                    print("match", frag)
                    print("old location", data['location']['fragments'][frag])

                    # Updates the data with the new location
                    data['location']['fragments'][frag] = row['Location']
                    print("new location", data['location']['fragments'][frag])
                else:
                    print("no match")

                # Adds to the counter as it iterates through the fragments
                if data['location']['fragments'][frag] != "":
                    print("has location")
                    num_locations = num_locations + 1
                else:
                    print("has no location")

            # Determines if the gene is build ready
            if num_locations == len(data['location']['fragments']):
                print("buildable")
                data['status']["build_ready"] = "TRUE"
            else:
                print("not buildable")

            # Writes a new json with the modified data
            with open(file,'w') as json_file:
                    json.dump(data,json_file,indent=2)
            print("cycle done")
            print("")
    print("plate {} done".format(plate))
print("Genes not in the dictionary: ")
print(not_in_dict)
print("Genes with 'link':")
print(linked_genes)
