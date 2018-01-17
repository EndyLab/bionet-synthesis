## Queries the database and pulls sequences to send to Sagacious to be analyzed

import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re

gene_ids = []
gene_names = []
sequences = []
max_num = 0

now = datetime.datetime.now()
date = "{}".format(now)
date = date[0:10]

if glob.glob("../sagacious/sag_submissions/*.csv"):
    print("previous submissions")
    for submission in glob.glob("../sagacious/sag_submissions/sagacious_order*.csv"):
        current_num = int(submission[46:48])
        if current_num > max_num:
            max_num = current_num
            sub_num = str(max_num + 1).zfill(3)
else:
    print("no previous builds")
    sub_num = '001'
counter = 0

path = "../sagacious/sag_submissions/sagacious_order_{}_{}.csv".format(sub_num,date)

# Reads into each of the files in the database to find candidates to send out
for file in glob.glob("../data/*/*.json"):

    #print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    # Only selects for genes that have been successfully built
    if data["status"]["build_complete"] != "Good_Sequence":
        continue

    # Checks if the genes have already been submitted
    if data["info"]["IP"]["submitted"]:
        print("submitted")
        continue

    counter += 1
    if counter == 51:
        break

    gene_ids.append(data["gene_id"])
    gene_names.append(data["gene_name"])
    sequences.append(data["sequence"]["optimized_sequence"])

    data["info"]["IP"]["submitted"] = True
    data["info"]["IP"]["submission_number"] = sub_num
    data["info"]["IP"]["results"] = "Awaiting results"


    with open(file,"w+") as json_file:
        json.dump(data,json_file,indent=2)

submission_data = pd.DataFrame({
    "Gene ID" : gene_ids,
    "Gene Name" : gene_names,
    "Sequence" : sequences
})
submission_data = submission_data[["Gene ID","Gene Name","Sequence"]]
submission_data.set_index("Gene ID")
submission_data.to_csv(path)





#
