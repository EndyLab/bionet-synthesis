## Generates a dictionary json file linking the new ID#'s to their gene names across the whole database
## Script to extract information from the entire database

import json
import os
import glob
import re

dictionary = {}
duplicates = []

for file in glob.glob("../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    gene = data["gene_name"]
    idnum = data["id"]
    if gene in dictionary:
        duplicates.append(gene)
    row = {gene: idnum}
    print(row)
    dictionary.update(row)

with open('./gene_id_dict.json','w') as json_file:
    json.dump(dictionary,json_file,indent=2)

print(duplicates)
