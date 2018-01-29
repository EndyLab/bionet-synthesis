import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

PATH = "../sequencing_files/"
target = "CM_530609_zipfile-1"

ref = pd.read_csv("../builds/build003_2018-01-22 12:58:16.csv")
genes = list(ref["Gene"])
rev_list = genes[::-1]
dictionary = dict(zip(genes, rev_list))
print(dictionary)

for seqfile in glob.glob("../sequencing_files/{}/*.ab1".format(target)):
    print(seqfile)

    if "Forward" in seqfile:
        initials, order_number, plate_number, well_number, sample_name, sample_number, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_([A-Za-z0-9]+)_([0-9]+)_M13-Forward---20-([A-H][0-9]{2}).ab1',
            seqfile).groups()
        primer = "_M13-Forward---20-"
    else:
        initials, order_number, plate_number, well_number, sample_name, sample_number, well_address = re.match(
            r'.*/([A-Z]+)_([0-9]+)-([0-9])([0-9]+)_([A-Za-z0-9]+)_([0-9]+)_M13-Reverse_([A-H][0-9]{2}).ab1',
            seqfile).groups()
        primer = "_M13-Reverse_"
    id_num = str(sample_name + "_" + sample_number)
    print(id_num)

    new_name = dictionary[id_num]
    print(new_name)

    new_file = str(PATH + target + "/" + initials + "_" + order_number + "-" + plate_number + well_number + "_" + new_name + primer + well_address + ".ab1")
    print(new_file)

    os.rename(seqfile, new_file)






#
