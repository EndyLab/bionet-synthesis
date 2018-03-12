## FOR GIT REPOSITORY: Takes in a file containing all of the fragment sequences and then appends the information to the json files

import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

import shutil
from config import *

PIPELINE_PATH = BASE_PATH + "/pipeline"
BUILDS_PATH = BASE_PATH + "/builds"
DATA_PATH = BASE_PATH + "/data"

BACKBONE_PATH = BASE_PATH + "/sequencing_files/popen_v1-1_backbone.fasta"
DICTIONARY_PATH = PIPELINE_PATH + "/testing/data_testing/10K_CDS.csv"

counter = 0

no_frag = []
linked_genes = []
not_in_dict = []
genes = []

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv(DICTIONARY_PATH)
dictionary = dict(zip(data['gene_name'], data['idnum']))

counter = 0
plate_map_number = []
plate = []
well = []
name = []
seq = []
not_in_dict = []
seq_disc = []
all_maps = []
no_frag_seq = []




## Take in all current plate maps
## ============================================
for file in glob.glob("{}/plate_maps/*.csv".format(BASE_PATH)):
    plate_map = pd.read_csv(file)
    all_maps.append(plate_map)

    
    order_num,attempt_num,date = re.match(
        r'.+\/O-([0-9]{3})_A-([0-9]{3})_([0-9.]+)\.csv',
        file).groups()
    print(file,order_num,attempt_num)
all_maps = pd.concat(all_maps)
input("stop")

## Generates the necessary dataframe
## ============================================
for index, row in all_maps.iterrows():
    # Searches for linked genes and separates them into separate entries
    if "_link_" in row['customer_line_item_id']:
        linked_string = row['customer_line_item_id'].strip()
        linked_genes = re.match(
            r'([A-Za-z0-9_-]+)_link_([A-Za-z0-9_-]+)',
            linked_string).groups()
        for gene in linked_genes:
            if gene[-2:] != "_1":
                gene += "_1"
            name.append(gene)
            plate.append(row['Plate'])
            well.append(row['Well'])
            seq.append(row['Insert Sequence'])
    else:
        name.append(row['customer_line_item_id'].strip())
        plate.append(row['Plate'])
        well.append(row['Well'])
        seq.append(row['Insert Sequence'])
    print("name:", name)

# Rebuilds the plate map with all of the linked genes split up
new_plate_map = pd.DataFrame({ "Name" : name, "Plate" : plate, "Well" : well, "Sequence" : seq })

# Converts all of the fragment names into their corresponding ID#'s using the dictionary created earlier
for index,row in new_plate_map.iterrows():
    if row['Name'][:6] != "BBF10K":
        if row['Name'][:-2] in dictionary.keys():
            frag_num = row['Name'][-1]
            idnum = dictionary[row['Name'][:-2]]
        else:
            print("{} not in dictionary".format(row['Name'][:-2]))
            not_in_dict.append(row['Name'])
            continue
    else:
        print("{} already in proper format".format(row['Name'][:-2]))
        idnum = row['Name'][:-2]
        frag_num = row['Name'][-1]

    # Uses the ID# to open the correct .json file and then adds the location
    for file in glob.glob("{}/{}/{}.json".format(DATA_PATH,idnum,idnum)):
        print(file)
        # Open the json file
        with open(file,"r") as json_file:
            data = json.load(json_file)
        data["location"]["fragments"][idnum + "_" + frag_num] = row['Plate'] + "_" + row['Well']

        if idnum + "_" + frag_num not in data["sequence"]["fragment_sequences"]:
            no_frag_seq.append(idnum)
            continue

        # Checks if the synthesized sequence matches the requested sequence
        if data["sequence"]["fragment_sequences"][idnum + "_" + frag_num] != row['Sequence']:
            print('Fragment sequence did not match the synthesized sequence')
            seq_disc.append(idnum)
            data["sequence"]["fragment_sequences"][idnum + "_" + frag_num] = row['Sequence']

        # Iterate through all of the fragments
        num_locations = 0
        for frag in data['location']['fragments']:
            # Adds to the counter as it iterates through the fragments
            if data['location']['fragments'][frag] != "":
                print("{} has location".format(frag))
                num_locations += 1
            else:
                print("has no location")

        # Determines if the gene is build ready
        if num_locations == len(data['location']['fragments']):
            print("buildable")
            data['status']["build_ready"] = True
        else:
            print("not buildable")

        with open(file,'w') as json_file:
                json.dump(data,json_file,indent=2)

print("Not in dictionary")
print(not_in_dict)
print("no fragment sequence")
print(no_frag_seq)












# something
