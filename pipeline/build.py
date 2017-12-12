## FOR GIT REPOSITORY -- Queries the database for build ready constructs and then sets up the plates

from opentrons import robot, containers, instruments
import argparse
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import json
import os
import glob
import re
import math

## Take in required information

targets = []

counter = 0
gene = 0
plate = 0
well = 0
row = 0
unique = 0
plates = ["pSHPs0807B412037MU", "pSHPs0807B412038MU"]
#plates = []

# Query the database and iterate through each json file
for file in glob.glob("../data/BBF10K*/*.json"):
    print(file)

    # Open and store the data within the json file
    with open(file,"r") as json_file:
        data = json.load(json_file)

    # Determine if it is build ready
    if data["status"]["build_ready"] != "TRUE":
        continue
    print("build ready")

    # Pull general information about the gene
    gene = data["id"]
    locations = data["location"]["fragments"]
    frag_num = len(locations)
    print("number of fragments: ", frag_num)

    number = 0

    # Iterate through all of the fragments
    for fragment in locations:
        print(fragment)

        # Pulls out the location for a specific fragment
        frag_loc = data["location"]["fragments"][fragment]

        # Separate the well from the plate and determine how many unique plates there are
        plate_loc, well = frag_loc.split("_")
        plates.append(plate_loc)
        plates_pd = pd.Series(plates)
        unique = len(plates_pd.unique())
        print("number on unique plates:", unique)

        # If there are too many unique plates it roles back the counter and series of plates and moves to the next file
        if unique > 2:

            counter = counter - 1 + number
            plates = plates[:-1]
            print("counter: ", counter)
            print("too many unique plates")
            continue
        else:
            row = [gene, plate_loc, well]
            targets.append(row)
        number = number + 1

    counter = counter + 1
    print("counter:", counter)
    if counter == 9:
        break

targets = np.array(targets)
plan = pd.DataFrame({
    "Plate" : targets[:,1],
     "Gene" : targets[:,0],
     "Well" : targets[:,2]
    })
#plan = plan.set_index('Plate')
#print()
#print("plates used: ", plan["Plate"].unique())

plan
