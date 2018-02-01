import numpy as np
import pandas as pd

import json
import os
import glob
import re

build_date = "2018-01-31 13:55:42"

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    if data["status"]["build_attempts"][0]["build_date"] == build_date:
        data["status"]["build_attempts"][0]["build_number"] = "build005"
        data["status"]["build_complete"] = "In_Process"

    with open(file,"w+") as json_file:
        json.dump(data,json_file,indent=2)








#
