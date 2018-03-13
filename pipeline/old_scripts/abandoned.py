## Queries the database and presents the status

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re

from config import *


complete_order = input("What is the last complete order: ")
abandoned = []
no_frag = []
not_abandoned = []

# Query every file in the database
for file in glob.glob(BASE_PATH + "/data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    idnum = data["gene_id"]
    order_num = data["info"]["order_number"]
    if int(order_num) <= int(complete_order) and data["status"]["build_ready"] != True:
        if data["location"]["fragments"] == {}:
            no_frag.append(idnum)
            continue
        abandoned.append(idnum)
        data["status"]["abandoned"] = True
    else:
        not_abandoned.append(idnum)
        data["status"]["abandoned"] = False

    with open(file,'w') as json_file:
        json.dump(data,json_file,indent=2)

print(no_frag)
print(len(no_frag))

print(not_abandoned)
print("Number of not abandoned genes: ",len(not_abandoned))

print(abandoned)
print("Number of abandoned genes: ",len(abandoned))
