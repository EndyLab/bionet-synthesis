import json
import os
import glob
import re

part_set=set()
#name = input("name to find (lowercase) : ")

for file in glob.glob("./../data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    part = data["author"]["name"]
    print(part)
    part_set.add(part)


print(len(part_set))
print(part_set)
