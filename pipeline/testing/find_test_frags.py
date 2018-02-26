import numpy as np
import pandas as pd

import json
import os
import glob
import re

names = []
cand = []
plate = []
well = []
frag_nums = []
gene_lengths = []
gene_ids = []

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv("./data_testing/10K_CDS.csv")
dictionary = dict(zip(data['gene_name'], data['idnum']))

previous_genes = []
for file in glob.glob("../../builds/*/build*_20*.csv"):
    print(file)
    build = pd.read_csv(file)
    previous_genes += list(build['Gene'])

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    if data["status"]["build_complete"] == "Good_sequence" and data["gene_id"] in previous_genes:
        gene_name = data["gene_name"]
        names.append(gene_name)

print("Names: ", names)
print(len(names))

total_plates = []

for plate_map in glob.glob("../../plate_maps/*.csv"):
    details = pd.read_csv(plate_map)
    print("plate_map", plate_map)
    for index,row in details.iterrows():
        gene_name = row["customer_line_item_id"][:-2].strip()
        print(gene_name)
        total_plates.append(row["Plate"])
        if gene_name in names:
            print("match")
            cand.append(gene_name)
            plate.append(row["Plate"])
            well.append(row["Well"])

candidates = pd.DataFrame({
        "Candidate" : cand,
        "Plate" : plate,
        "Well" : well
    })

for index, row in candidates.iterrows():
    for file in glob.glob("../../data/{}/{}.json".format(dictionary[row["Candidate"]],dictionary[row["Candidate"]])):
        frag_num = 0
        print(file)
        with open(file,"r") as json_file:
            data = json.load(json_file)
        print(data["location"]["fragments"])
        for fragment in data["location"]["fragments"]:
            frag_num += 1
        frag_nums.append(str(frag_num))
        length = len(data["sequence"]["optimized_sequence"])
        gene_lengths.append(length)
        gene_ids.append(data["gene_id"])

candidates["Fragments"] = frag_nums
candidates["Length"] = gene_lengths
candidates["Gene ID"] = gene_ids

candidates = candidates[["Candidate","Gene ID","Plate","Well","Fragments","Length"]]

print(candidates)

unique_plates = pd.unique(candidates["Plate"])
print(unique_plates)

complete_per_plate = []

for plate in unique_plates:
    complete = len(candidates[candidates.Plate == plate])
    complete_per_plate.append(complete)
    print("{} has {} genes left".format(plate, complete))

plate_breakdown = pd.DataFrame({
    "Plate" : unique_plates,
    "Successful_constructs" : complete_per_plate
})
print(plate_breakdown)
candidates.to_csv("./complete_constructs.csv")


#
