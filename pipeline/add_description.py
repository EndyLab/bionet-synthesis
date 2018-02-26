
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

import datetime
from config import *

data = pd.read_csv(BASE_PATH + "pipeline/testing/data_testing/10K_CDS.csv")

for index,row in data.iterrows():
    idnum = row["idnum"]
    print(idnum)
    description = row["description"]

    for file in glob.glob(BASE_PATH + "/data/{}/{}.json".format(idnum,idnum)):
        print(file)
        with open(file,"r") as json_file:
            data = json.load(json_file)
        print(description)
        print(type(description))
        if type(description) is str:
            data["description"] = description
        else:
            data["description"] = "None"

        with open(file,"w+") as json_file:
            json.dump(data,json_file,indent=2)
