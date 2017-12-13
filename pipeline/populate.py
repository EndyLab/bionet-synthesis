import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re
import math

import shutil

counter = 0


data = pd.read_csv("./testing/data_testing/10K_CDS.csv")

for file in glob.glob("./testing/json/template.json"):
        print(file)
        with open(file,"r") as template_json:
            template = json.load(template_json)
            print(template)

for index, row in data.iterrows():
    counter = counter + 1
    if counter == 11:
        break
    gene = row['Gene']
    seq = row['Sequence']
    number = str(index + 1).zfill(6)
    idnum = "BBF10K_" + number
    path = "../old_files/test_db_upload/{}".format(idnum)
    if os.path.exists(path):
        print("Directory {} already exists".format(idnum))
    else:
        os.makedirs(path)

        fasta = open("./{}/{}.fasta".format(path,idnum),"w+")
        fasta.write(">{}".format(gene))
        fasta.write("\n")
        fasta.write(seq)
        fasta.close()

        #json = template
        template["gene_id"] = idnum
        print(template["gene_id"])
        template["gene_name"] = gene
        print(template["gene_name"])
        template["sequence"]["original_sequence"] = seq
        with open("./{}/{}.json".format(path,idnum),"w+") as json_file:
            json.dump(template,json_file,indent=2)
            
