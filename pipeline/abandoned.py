## Queries the database and presents the status

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re

complete_order = input("What is the last complete order: ")
abandoned = []
no_frag = []

# Query every file in the database
for file in glob.glob("../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    idnum = data["gene_id"]
    order_num = data["info"]["order_number"]
    if int(order_num) <= int(complete_order) and data["status"]["build_ready"] != "TRUE":
        if data["location"]["fragments"] == {}:
            no_frag.append(idnum)
            continue
        abandoned.append(idnum)
        data["status"]["abandoned"] = "TRUE"
        with open(file,'w') as json_file:
            json.dump(data,json_file,indent=2)

print(no_frag)
print(len(no_frag))

print(abandoned)
print("Number of abandoned genes: ",len(abandoned))

    #if data["status"]["build_ready"] == "TRUE":
    #    continue
