import numpy as np
import pandas as pd

import json
import os
import glob
import re

for file in glob.glob("../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    data["status"]["build_complete"] = ""

    with open(file,"w+") as json_file:
        json.dump(data,json_file,indent=2)








#
