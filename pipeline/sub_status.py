## Queries the database and presents the status

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re

info = []
contributors = []

ordered = 0
abandoned = 0
build_ready = 0
attempted = 0
unknown = 0
vector = 0
complete = 0

total_subs = 6

for sub in range(total_subs):
    sub += 1
    print("sub", sub)

    info = []
    contributors = []

    ordered = 0
    abandoned = 0
    build_ready = 0
    attempted = 0
    unknown = 0
    vector = 0
    complete = 0

    for file in glob.glob("../data/*/*.json"):
        #print(file)
        with open(file,"r") as json_file:
            data = json.load(json_file)

        if data["info"]["order_number"] != sub:
            continue

        author = data["author"]["name"]
        contributors.append(author)
        if data["status"]["ordered"] == True:
            ordered += 1
        if data["status"]["abandoned"] == True:
            abandoned += 1
        if data["status"]["build_ready"] == True:
            build_ready += 1
        if data["status"]["build_complete"] != "":
            attempted += 1
        if data["status"]["build_complete"] == "Good_Sequence":
            complete += 1
        elif data["status"]["build_complete"] == "Original_Vector_Sequence":
            vector += 1
        elif data["status"]["build_complete"] == "Unknown_Sequence":
            unknown += 1

    contributors = pd.Series(contributors)
    unique_cont = len(contributors.unique())
    production = ordered - (build_ready + abandoned)
    failures = (vector + unknown)
    not_attempted = (build_ready - attempted)
    print("Contributers :", unique_cont)
    print("Ordered :", ordered)
    print("Abandoned :", abandoned)
    print("Received :", build_ready)
    print("In Production : ", production)
    print("Build Attempted :", attempted)
    print("Verified Success :", complete)
    print("Failures :", failures)
    print("Not Yet Attempted :", not_attempted)
    print()
