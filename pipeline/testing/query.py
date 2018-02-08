import numpy as np
import pandas as pd

import json
import os
import glob
import re

id_nums = []
names = []
plates = []

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data["author"]["affiliation"] == "Addgene":
        print("positive")
        id_num = data["gene_id"]
        id_nums.append(id_num)
        gene_name = data["gene_name"]
        names.append(gene_name)
        plate = data["location"]["fragments"].values()
        plates.append(plate)

results = pd.DataFrame({
    "Gene ID" : id_nums,
    "Name" : names,
    "Plate" : plates
})
results.to_csv("./addgene_samples.csv")
print(results)






#
