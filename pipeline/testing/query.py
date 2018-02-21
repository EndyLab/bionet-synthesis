import numpy as np
import pandas as pd

import json
import os
import glob
import re

id_nums = []
names = []
plates = ["pSHPs0826B426849MU","pSHPs0807B412039MU"]

# Pulls in all of the previously attempted genes
previous_genes = []
for file in glob.glob("../../builds/*/build*_20*.csv"):
    print(file)
    build = pd.read_csv(file)
    previous_genes += list(build['Gene'])

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    #if data["author"]["affiliation"] == "Addgene":
    if data["gene_id"] in previous_genes:
        continue
    if data["status"]["build_ready"] != True:
        continue
    for fragment in data["location"]["fragments"]:
        plate = data["location"]["fragments"][fragment].split("_")[0]
        print("fragment",plate)
        if plate == plates[0] or plate == plates[1]:
            print("positive")
            id_num = data["gene_id"]
            id_nums.append(id_num)
            gene_name = data["gene_name"]
            names.append(gene_name)

#names = pd.unique(names)
#id_nums = pd.unique(id_nums)

results = pd.DataFrame({
    "Gene ID" : id_nums,
    "Name" : names
})
results.to_csv("./select_plates.csv")
print(results)






#
