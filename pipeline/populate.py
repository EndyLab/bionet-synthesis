import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

import shutil
import datetime
from config import *

counter = 0
last_id = 0
idnum_list = []

# Imports the csv containing all of our data
data = pd.read_csv(BASE_PATH + "/pipeline/testing/data_testing/10K_CDS.csv")

# Takes in the template json file to initialize the json files
for file in glob.glob(BASE_PATH + "/pipeline/testing/json/template.json"):
        print(file)
        with open(file,"r") as template_json:
            template = json.load(template_json)

# Iterates through all of the items in the spreadsheet
for index, row in data.iterrows():

    # Checks if it has already been assigned a number
    #if row['idnum'] != "None":
    #    print(row['idnum'][7:13])
    #    last_id = int(row["idnum"][7:13])
    #    idnum_list.append(row["idnum"])
    #    print("BBF10K_{} has already been entered".format(str(last_id).zfill(6)))
    #    continue

    # Generate the ID#
    gene = row['gene_name']
    seq = row['sequence']
    idnum = row['idnum']
    #last_id += 1
    #number = str(last_id).zfill(6)
    #idnum = "BBF10K_" + number
    #idnum_list.append(idnum)

    #State the path to house the new set of directories
    path = "{}/data/{}".format(BASE_PATH,idnum)

    # Only create a new directory if there are no existing directories
    if os.path.exists(path):
        print("Directory {} already exists".format(idnum))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)
        print("Making directory for {}".format(idnum))

        # Generate the fasta file with the gene sequence
        fasta = open("/{}/{}.fasta".format(path,idnum),"w+")
        fasta.write(">{}".format(gene))
        fasta.write("\n")
        fasta.write(seq)
        fasta.close()

        # Fill in the template json file with the information from each row
        template["gene_id"] = idnum
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

        # Takes in the date that it was ordered and calculates the week number
        date = row["date_ordered"]
        year, month, day = date.split(".")
        year = int(year)
        month = int(month)
        day = int(day)
        week = datetime.date(year, month, day).isocalendar()[1]
        order_week = str(year) + "." + str(week)
        template["info"]["order_week"] = order_week

        template["info"]["order_number"] = row["order_number"]
        template["dates"]["submitted"] = row["date_ordered"]
        template["dates"]["ordered"] = row["date_ordered"]

        # Writes the information to a json file in the desired directory
        with open("/{}/{}.json".format(path,idnum),"w+") as json_file:
            json.dump(template,json_file,indent=2)

# Update the .csv file with the new id numbers
#data["idnum"] = idnum_list
#new_data = data[["gene_name","sequence","author","author_email","author_affiliation","author_project","cloning_method","part_type","build_type","safety","collection","other_tags","ordered","will_build","date_ordered","order_number","idnum"]]
#new_data.to_csv('./testing/data_testing/10K_CDS.csv')
