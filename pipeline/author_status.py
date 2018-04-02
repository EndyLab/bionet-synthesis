import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re


author_set=set()

for file in glob.glob("./../data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    author = data["author"]["name"]
    author_set.add(author)

print(len(author_set))

