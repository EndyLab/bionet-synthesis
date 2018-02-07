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

data = pd.read_csv("build000.csv")

new_dir = "../../sequencing_files/build000_seq_files"

if os.path.exists(new_dir):
    print("{} already exists".format(new_dir))
else:
    os.makedirs(new_dir)
    print("Made directory")

input("Continue")

for index, row in data.iterrows():
    name = row["Gene Name"]

    for file in glob.glob("../../sequencing_files/*/*{}*.ab1".format(name)):
        print(file)
        shutil.copy(file, new_dir)






#
