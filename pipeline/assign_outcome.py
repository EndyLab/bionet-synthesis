from datetime import datetime
start = datetime.now()
print("Starting run at: ",start)

import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

import shutil

import_time = datetime.now()
print("Time to import: ",import_time - start)

SEQ_FILE = "CM_530609_zipfile-1"
RESULTS_PATH = "../sequencing_files/{}/alignment_results.csv".format(SEQ_FILE)

results = pd.read_csv(RESULTS_PATH)
print(results)

for index, row in results.iterrows():
    with open("../data/{}/{}.json".format(row["Gene"],row["Gene"),"r") as json_file:
        data = json.load(json_file)
    for attempts in data["status"]["build_attempts"]:
        print(attempts["build_number"])






#
