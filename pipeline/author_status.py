import json
import os
import glob
import re

part_set=set()
#name = input("name to find (lowercase) : ")

for file in glob.glob("./../data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    #part = data["author"]["name"]
    #if data["sequence"]["fragment_sequences"] == "": 
    if str(data["sequence"]["fragment_sequences"]) == "{}":
        print(data["gene_id"] + " " + data["info"]["documentation"]["gene_name"] + " " + str(data["info"]["order_number"]))


#print(part_set)
