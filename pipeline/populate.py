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

# Imports the csv containing all of our data
data = pd.read_csv("./testing/data_testing/10K_CDS.csv")

# Takes in the
for file in glob.glob("./testing/json/template.json"):
        print(file)
        with open(file,"r") as template_json:
            template = json.load(template_json)
            #print(template)

# Iterates through all of the items in the spreadsheet
for index, row in data.iterrows():

    # Forces it to only run a certain number of times
    #counter = counter + 1
    #if counter == 11:
    #    break

    # Generate the ID#
    gene = row['gene_name']
    seq = row['sequence']
    number = str(index + 1).zfill(6)
    idnum = "BBF10K_" + number

    #State the path to house the new set of directories
    path = "../data/{}".format(idnum)

    # Only create a new directory if there are no existing directories
    if os.path.exists(path):
        print("Directory {} already exists".format(idnum))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)

        # Generate the fasta file with the gene sequence
        fasta = open("./{}/{}.fasta".format(path,idnum),"w+")
        fasta.write(">{}".format(gene))
        fasta.write("\n")
        fasta.write(seq)
        fasta.close()

        # Fill in the template json file with the information from each row
        template["gene_id"] = idnum
        print(idnum)
        template["gene_name"] = gene
        template["sequence"]["original_sequence"] = seq
        template["sequence"]["optimized_sequence"] = seq
        template["author"]["name"] = row["author"]
        template["author"]["email"] = row["author_email"]
        template["author"]["affiliation"] = row["author_affiliation"]
        template["author"]["project"] = row["author_project"]
        template["info"]["type"]["cloning_method"] = row["cloning_method"]
        template["info"]["type"]["part_type"] = row["part_type"]
        template["info"]["type"]["build_type"] = row["build_type"]
        template["info"]["safety"] = row["safety"]
        template["info"]["collection"] = row["collection"]
        template["info"]["other_tags"] = row["other_tags"]
        template["status"]["ordered"] = row["ordered"]
        template["status"]["will_build"] = row["will_build"]
        template["dates"]["ordered"] = row["date_ordered"]

        with open("./{}/{}.json".format(path,idnum),"w+") as json_file:
            json.dump(template,json_file,indent=2)

    
