import numpy as np
import pandas as pd

import json
import os
import glob
import re

# from config import *

all_jcvi = []
jcvi = []
good = []
jc_good = []
jc_point = []
jc_bad = []
jc_blank = []
jc_pro = []
other = []

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    id_num = data["gene_id"]
    if "MMSYN1" in data["gene_name"]:
        all_jcvi.append(id_num)
    if data["status"]["build_ready"] == True:
        if "MMSYN1" in data["gene_name"]:
            jcvi.append(id_num)
            if "Good" in data["status"]["build_complete"]:
                jc_good.append(id_num)
            elif "Point" in data["status"]["build_complete"]:
                jc_point.append(id_num)
            elif "Original" in data["status"]["build_complete"]:
                jc_bad.append(id_num)
            elif data["status"]["build_complete"] == "":
                if data["status"]["build_attempts"][0]["build_outcome"] == "In Process":
                    jc_pro.append(id_num)
                else:
                    jc_blank.append(id_num)
            else:
                other.append(data["status"]["build_complete"])
        if "Good" in data["status"]["build_complete"]:
            good.append(id_num)
print(jcvi)
print(good)
print(jc_good)
print()
print(other)
print("\n\n",jc_blank,"\n\n")


print("all_jcvi",len(all_jcvi),"JCVI: ",len(jcvi),"Good: ",len(good),"JC_Good",len(jc_good),"JC_point",len(jc_point),"Jc_bad",len(jc_bad),"Jc_blank",len(jc_blank),"jc_pro",len(jc_pro),"other",len(other))
