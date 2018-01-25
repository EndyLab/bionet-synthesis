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
import datetime

number = input("BBF10K_ ...? ")
id_num = "BBF10K_{}".format(number.zfill(6))
print(id_num)

# Query the database and iterate through each json file
for file in glob.glob("../data/{}/{}.json".format(id_num,id_num)):
    print(file)

    # Open and store the data within the json file
    with open(file,"r") as json_file:
        data = json.load(json_file)
        print(data)
    gene_name = data["gene_name"]
    author_name = data["author"]["name"]
    original_sequence = data["sequence"]["original_sequence"]
    optimized_sequence = data["sequence"]["optimized_sequence"]

    if data["status"]["build_complete"] == "Good_Sequence":
        status = "Sequence Verified"
    elif data["status"]["build_complete"] == "":
        status = "Build not yet attempted"
    else:
        status = "Build attempt not successful"

    submitted = data["dates"]["submitted"]
    ordered = data["dates"]["ordered"]

    if data["dates"]["build_ready"] == "":
        build_ready = "Not yet"
    else:
        build_ready = data["dates"]["build_ready"]

    if data["dates"]["build_ready"] == "":
        complete = "Not yet"
    else:
        complete = data["dates"]["build_ready"]

    now = datetime.datetime.now()
    now = "{}".format(now)

    path = "../docs/{}".format(id_num)

    if os.path.exists(path):
        print("Directory {} already exists".format(id_num))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)

    # Generate the fasta file with the gene sequence
    index = open("{}/index.md".format(path),"w+")
    index.write("# {}\n".format(id_num))
    index.write("### Gene Name: **{}**\n".format(gene_name))
    index.write("### Author Name: **{}**\n".format(author_name))
    index.write("### Current Status: **{}**\n".format(status))
    index.write("\n")
    index.write("Status breakdown:")
    index.write("\n\n")
    index.write("Step | Date ")
    index.write("\n")
    index.write("--- | --- ")
    index.write("\n")
    index.write("Submitted | {} ".format(submitted))
    index.write("\n")
    index.write("Ordered | {} ".format(ordered))
    index.write("\n")
    index.write("Received | {} ".format(build_ready))
    index.write("\n")
    index.write("Sequence Verified | {} ".format(complete))
    index.write("\n\n")
    index.write("##### Last updated: {}".format(now[:-16]))
    index.write("\n")
    index.write("### Sequences:\n")
    index.write("* Original: \n")
    index.write("   {} \n".format(original_sequence))
    index.write("* Optimized: \n")
    index.write("   {} \n".format(optimized_sequence))
    index.write("\n")
    index.close()
