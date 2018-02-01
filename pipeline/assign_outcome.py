from datetime import datetime
start = datetime.now()
print("Starting run at: ",start)

import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

import shutil

import_time = datetime.now()
print("Time to import: ",import_time - start)

build_num = "build003"

SEQ_FILE = "{}_seq_files".format(build_num)
RESULTS_PATH = "../builds/{}/{}_alignment_results.csv".format(build_num,build_num)

results = pd.read_csv(RESULTS_PATH)
print(results)

for index, row in results.iterrows():
    if row["Target"] != "gene":
        continue
    if row["Manual"] == "0":
        row["Manual"] = row["Outcome"]
        if row["Manual"] == "Perfect":
            row["Manual"] = "Good_sequence"
    elif row["Manual"] == "1":
        row["Manual"] = "Point Mutation"
    elif row["Manual"] == "Perfect":
        row["Manual"] = "Good_sequence"

    with open("../data/{}/{}.json".format(row["Gene"],row["Gene"]),"r") as json_file:
        data = json.load(json_file)
    for attempts in data["status"]["build_attempts"]:
        attempts["build_number"] = "build003"
        attempts["build_outcome"] = row["Manual"]
    data["status"]["build_complete"] = row["Manual"]

    with open("../data/{}/{}.json".format(row["Gene"],row["Gene"]),"w+") as json_file:
        json.dump(data,json_file,indent=2)

stop = datetime.now()
runtime = stop - start
print("Total runtime is: ", runtime)






#
