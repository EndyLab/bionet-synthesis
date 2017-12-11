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

for file in glob.glob("./gene_id_dict.json"):
        print(file)
        with open(file,"r") as json_file_dict:
            print(json_file_dict)
            dictionary = json.load(json_file_dict)

#print(dictionary)
for plate in glob.glob("../plate_maps/*.csv"):
    plates = pd.read_csv(plate)

    plates['customer_line_item_id'] = plates['customer_line_item_id'].str.strip()

    ids = plates['customer_line_item_id'].str.rsplit('_', n=1)
    plates['Gene'] = ids.str[0]
    plates['Fragment'] = ids.str[1]
    plates['Location'] = plates['Plate'] + "_" + plates["Well"]

    for index, row in plates.iterrows():
        print(index)
        name = row['Gene']
        print(name)

        if 'link' in name:
            print('link')
            name = name[:11]

        if name not in dictionary:
            continue

        gene = dictionary[name]
        print(gene)

        fragment_num = row['Fragment']
        print(fragment_num)

        fragment = gene + "_" + fragment_num
        print(fragment)

        for file in glob.glob("../data/{}/{}.json".format(gene,gene)):
            num_locations = 0
            print(file)
            with open(file,"r") as json_file:
                data = json.load(json_file)
            for frag in data['location']['fragments']:
                if frag == fragment:
                    print(frag)
                    print("match")
                    print(data['location']['fragments'][frag])
                    data['location']['fragments'][frag] = row['Location']
                    print(data['location']['fragments'][frag])
                else:
                    print("no match")

                if data['location']['fragments'][frag] != "":
                    print("has location")
                    num_locations = num_locations + 1
                else:
                    print("has no location")

                print(num_locations)
                print(len(data['location']['fragments']))

            if num_locations == len(data['location']['fragments']):
                print("buildable")
                data['status']["build_ready"] = "TRUE"
            else:
                print("not buildable")
            with open(file,'w') as json_file:
                    json.dump(data,json_file,indent=2)
