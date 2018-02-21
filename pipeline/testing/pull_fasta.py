import numpy as np
import pandas as pd

import json
import os
import glob
import re

import shutil

build_num = "build006"

built = []
for file in glob.glob("../../builds/{}/{}_alignment_results.csv".format(build_num,build_num)):
    print(file)
    results = pd.read_csv(file)
    built += list(results['Gene'])

new_dir = "./{}/{}_fasta_files".format(build_num,build_num)
if os.path.exists(new_dir):
    print("{} already exists".format(new_dir))
else:
    os.makedirs(new_dir)


for gene in built:
    for fasta in glob.glob("../../data/{}/{}.fasta".format(gene,gene)):
        shutil.copy(fasta, new_dir)










#
