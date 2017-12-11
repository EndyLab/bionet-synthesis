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

ordered = 0
build_ready = 0
attempted = 0
unknown = 0
vector = 0
complete = 0


for file in glob.glob("../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)

    if data["status"]["ordered"] == "TRUE":
        ordered = ordered + 1
    if data["status"]["build_ready"] == "TRUE":
        build_ready = build_ready + 1
    if data["status"]["build_complete"] != "":
        attempted = attempted + 1
    if data["status"]["build_complete"] == "Good_Sequence":
        complete = complete + 1
    elif data["status"]["build_complete"] == "Original_Vector_Sequence":
        vector = vector + 1
    elif data["status"]["build_complete"] == "Unknown_Sequence":
        unknown = unknown + 1

print("Ordered :", ordered)
print("Received :", build_ready)
print("Build Attempted :", attempted)
print("Verified Success :", complete)
print("Failures :", (vector + unknown))
