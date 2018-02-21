import numpy as np
import pandas as pd

import json
import os
import glob
import re

id_nums = []
names = []
authors = []
seqs = []

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    id_nums.append(data["gene_id"])
    names.append(data["gene_name"])
    authors.append(data["author"]["affiliation"])
    seqs.append(data["sequence"]["optimized_sequence"])

results = pd.DataFrame({
    "Gene ID" : id_nums,
    "Name" : names,
    "Author" : authors,
    "Sequence" : seqs
})

results = results[["Gene ID","Name","Author","Sequence"]]
results.to_csv("./sagacious_data.csv")
print(results)






#
