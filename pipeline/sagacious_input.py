## Queries the database and pulls sequences to send to Sagacious to be analyzed
import argparse
import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re

parser = argparse.ArgumentParser(description="Pull together sequences to submit to Sagacious.")
parser.add_argument('-w', '--write', required=False, action="store_false", help="Writes to all of the data files.")
args = parser.parse_args()

gene_ids = []
gene_names = []
sequences = []
max_num = 0

#num_entries = input("How many genes to be submitted: ")
num_entries = 49

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
    print("no previous submissions")
    sub_num = '001'
counter = 0

order_path = "../sagacious/sag_submissions/sagacious_order_{}_{}.csv".format(sub_num,date)

targets = ["BBF10K_000625","BBF10K_000562","BBF10K_000563","BBF10K_000564","BBF10K_000565","BBF10K_000566","BBF10K_000567","BBF10K_000568","BBF10K_000569","BBF10K_000570","BBF10K_000571","BBF10K_000572","BBF10K_000573","BBF10K_000574","BBF10K_000575","BBF10K_000576","BBF10K_000577","BBF10K_000578","BBF10K_000579","BBF10K_000580","BBF10K_000581","BBF10K_000582","BBF10K_000635","*"]

# Reads into each of the files in the database to find candidates to send out
for entry in targets:
    for file in glob.glob("../data/{}/{}.json".format(entry,entry)):
        print(file)
        with open(file,"r") as json_file:
            data = json.load(json_file)
        print(data["gene_name"])

        # Checks if the genes have already been submitted
        if data["info"]["IP"]["submitted"]:
            continue
        if counter == num_entries:
            break
        gene_ids.append(data["gene_id"])
        gene_names.append(data["gene_name"])
        sequences.append(data["sequence"]["optimized_sequence"])

        data["info"]["IP"]["submitted"] = True
        data["info"]["IP"]["submission_number"] = sub_num
        data["info"]["IP"]["results"] = "Awaiting results"

        counter += 1

        if args.write == True:
            with open(file,"w+") as json_file:
                json.dump(data,json_file,indent=2)

submission_data = pd.DataFrame({
    "Gene ID" : gene_ids,
    "Gene Name" : gene_names,
    "Sequence" : sequences
})
submission_data = submission_data[["Gene ID","Gene Name","Sequence"]]
submission_data.set_index("Gene ID")
submission_data.to_csv(order_path)
print(order_path)





#
