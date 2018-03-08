import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re

counter = 0

for file in glob.glob("../data/*/*.json"):
    print(file)
    #counter += 1
    #if counter == 11:
    #    break
    #print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    #del data["info"]["IP"]
    #print("deleted")

    data["info"]["IP"] = {}
    data["info"]["IP"]["submitted"] = False
    data["info"]["IP"]["submission_number"] = ""
    data["info"]["IP"]["results"] = ""

    with open(file,"w+") as json_file:
        json.dump(data,json_file,indent=2)
