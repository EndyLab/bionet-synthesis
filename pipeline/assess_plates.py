import pandas as pd
import json
import glob
import numpy as np
import math

from config import *
import interact_db

##
### Type 'e' in the results column of the build map to specify that it is
### empty or has no pickable colonies.
##

def find_map(path):
    '''Checks to make sure that the build exists and if not asks for a different one'''
    if glob.glob(path):
        return glob.glob(path)
    else:
        build_num = input("Not a valid build number. Please try again: ")
        build_name = "build{}".format(str(build_num).zfill(3))
        path = "{}/builds/{}/{}_20*.csv".format(BASE_PATH,build_name,build_name)
        return find_map(path)

# Determine which build to assess
build_num = input("Enter build number: ")
build_name = "build{}".format(str(build_num).zfill(3))
path = "{}/builds/{}/{}_20*.csv".format(BASE_PATH,build_name,build_name)
csv_path = find_map(path)

# Load the build map
build_map = pd.read_csv(csv_path)
build_map["Result"] = build_map["Result"].str.replace("e","Replace")
build_map["Result"] = build_map["Result"].str.replace("Pick","Replace")
build_map["Result"] = build_map["Result"].fillna("Good")

# Query the database for all current entries in a dataframe
db = interact_db.query_db()

# Refine the dataframe to only include candidates that should be repicked
mut_db = db[['data.gene_id','data.status.build_ready','data.status.abandoned','data.status.build_complete','data.location.fragments']]
exclude = ["Good_sequence","","Original_vector_sequence","Original Vector Sequence","Unknown_sequence","In_process","Unknown Sequence"]
for result in exclude:
    mut_db = mut_db[mut_db['data.status.build_complete'] != result]

# Pull out all of the mutated constructs
mut_db = mut_db.reset_index()
mut_list = mut_db["data.gene_id"].tolist()

# Maintain good constructs while replacing bad/empty ones with mutated
# constructs. If it runs out of mutated constructs, it will pick duplicates
picked = []
plates = []
wells = []
mut_counter = 0
for index,row in build_map.iterrows():
    if mut_counter >= len(mut_list):
        print("start over")
        mut_counter = 0
    if row["Result"] == "Good":
        picked.append(row["Gene"])
    elif row["Result"] == "Replace":
        picked.append(mut_list[mut_counter])
        mut_counter += 1

# Pull the loation information from the json files
build_nums = []
wells = []
for gene in picked:
    with open("{}/data/{}/{}.json".format(BASE_PATH,gene,gene),"r") as json_file:
        data = json.load(json_file)
    build_nums.append(data["status"]["build_attempts"][-1]["build_number"])
    wells.append(data["status"]["build_attempts"][-1]["build_well"])

# Append them to the dataframe
build_map["Picked"] = picked
build_map["Build"] = build_nums
build_map["Well"] = wells

# Sort the constructs based on build so that they are easier to pick
build_map = build_map.sort_values(["Build"],ascending=False)

build_map.to_csv(csv_path)
