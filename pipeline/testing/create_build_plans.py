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

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv("./data_testing/10K_CDS.csv")
dictionary = dict(zip(data['gene_name'], data['idnum']))

# Pulls in all of the previously attempted genes
previous_genes = []
for file in glob.glob("../../builds/*/build*_20*.csv"):
    print(file)
    build = pd.read_csv(file)
    previous_genes += list(build['Gene'])

# Extracts all of the gene names that are currently buildable and have not been attempted
for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data["status"]["build_ready"] == True and data["gene_id"] not in previous_genes:
        gene_name = data["gene_name"]
        names.append(gene_name)

print(previous_genes)
print("Names: ", names)
print(len(names))

input("continue?")

unknown_names = []

# Pulls all of the plate_maps and stores the locations of all of the flagged fragments
for plate_map in glob.glob("../../plate_maps/*.csv"):
    details = pd.read_csv(plate_map)
    print("plate_map", plate_map)
    for index,row in details.iterrows():
        if re.match(r'.+_[0-9]',row["customer_line_item_id"]):
            gene_name = row["customer_line_item_id"][:-2].strip()
        else:
            unknown_names.append(row["customer_line_item_id"])
        print(gene_name)
        if gene_name in names and gene_name not in previous_genes:
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
        #print(data["location"]["fragments"])
        for fragment in data["location"]["fragments"]:
            frag_num += 1
        frag_nums.append(str(frag_num))

candidates["Fragments"] = frag_nums

candidates = candidates[["Candidate","Plate","Well","Fragments"]]

print(candidates)

unique_plates = pd.unique(candidates["Plate"])
print(unique_plates)

remaining_per_plate = []

for plate in unique_plates:
    remaining = len(candidates[candidates.Plate == plate])
    remaining_per_plate.append(remaining)

plate_breakdown = pd.DataFrame({
    "Plate" : unique_plates,
    "Remaining_constructs" : remaining_per_plate
})
print(plate_breakdown)

#candidates.to_csv("./remaining_constructs.csv")

def find_combinations(list, sum):
    if not list:
        if sum == 0:
            return [[]]
        return []
    return find_combinations(list[1:], sum) + \
        [[list[0]] + tail for tail in
         find_combinations(list[1:], sum - list[0])]

x = 0
result_counter = 0
target = 96
max_plates = 3
#numbers = [5, 53, 19, 2, 17, 2, 23, 90, 75, 94, 15]
#numbers = list(plate_breakdown["Remaining_constructs"])
numbers = remaining_per_plate

combinations = []
reactions = []

while x == 0:
    result = find_combinations(numbers, target)
    if len(result) > 0:
        for comb in result:
            if len(comb) <= max_plates:
                combinations.append(comb)
                reactions.append(target)
                result_counter += 1
        else:
            target -= 1
    else:
        target -= 1
    #if len(combinations) > 15:
    #    break
    if target == 50:
        break

solution = pd.DataFrame({
    "Reactions" : reactions,
    "Combinations" : combinations
})

solution = solution[["Reactions","Combinations"]]

print(solution)

for index, row in solution.iterrows():
    current_comb = list(row["Combinations"])
    plates = plate_breakdown.loc[plate_breakdown['Remaining_constructs'].isin(current_comb)]
    if len(plates) > max_plates:
        print("Unique Total reactions = ", plates['Remaining_constructs'].unique().sum())
    else:
        if plates['Remaining_constructs'].sum() > target:
            print("Unique Total reactions = ", plates['Remaining_constructs'].unique().sum())
        else:
            print("Total reactions = ", plates['Remaining_constructs'].sum())
    print(plates,"\n")






#
