import numpy as np
import pandas as pd

import json
import os
import glob
import re

id_nums = []
names = []
build_type = []
built_ids = []
built_names = []
all_seq_ids = []
plates = ["pSHPs0826B426849MU","pSHPs0807B412039MU"]

# Pulls in all of the previously attempted genes
# previous_genes = []
# for file in glob.glob("../../builds/*/build*_20*.csv"):
#     print(file)
#     build = pd.read_csv(file)
#     previous_genes += list(build['Gene'])

path = "./Addgene_fasta"
os.makedirs(path)

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data["author"]["affiliation"] != "Addgene":
        continue
    # if data["gene_id"] in previous_genes:
        # continue
    all_seq_ids.append(data["gene_id"])
    if data["status"]["build_complete"] != "Good_sequence":
        continue
    built_ids.append(data["gene_id"])
    built_names.append(data["gene_name"])
    with open("{}/{}.fasta".format(path,data["gene_id"]),"w+") as fasta:
        fasta.write(">{}\n{}".format(data["gene_id"],data["sequence"]["optimized_sequence"]))




    # if data["status"]["build_ready"] != True:
    #     continue
    # if data["info"]["type"]["build_type"] != "10K_MoClo-EntryCDS-BbsI":
    #     id_nums.append(data["gene_id"])
    #     build_type.append(data["info"]["type"]["build_type"])

    # for fragment in data["location"]["fragments"]:
    #     plate = data["location"]["fragments"][fragment].split("_")[0]
    #     print("fragment",plate)
    #     if plate == plates[0] or plate == plates[1]:
    #         print("positive")
    #         id_num = data["gene_id"]
    #         id_nums.append(id_num)
    #         gene_name = data["gene_name"]
    #         names.append(gene_name)

#names = pd.unique(names)
#id_nums = pd.unique(id_nums)
print(all_seq_ids)
print()
print(built_ids)
print(built_names)

results = pd.DataFrame({
    "Gene ID" : built_ids,
    "Gene Name" : built_names
})
results.to_csv("./Addgene_genes.csv")
print(results)






#
