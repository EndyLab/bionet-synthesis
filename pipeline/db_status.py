## Queries the database and presents the status

import datetime

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

sub1 = 0
sub2 = 0
sub3 = 0
sub4 = 0
sub5 = 0

for file in glob.glob("../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    author = data["author"]["name"]
    contributors.append(author)

    if data["info"]["order_number"] == 1:
        sub1 = sub1 + 1
    elif data["info"]["order_number"] == 2:
        sub2 = sub2 + 1
    elif data["info"]["order_number"] == 3:
        sub3 = sub3 + 1
    elif data["info"]["order_number"] == 4:
        sub4 = sub4 + 1
    elif data["info"]["order_number"] == 5:
        sub5 = sub5 + 1

    if data["status"]["ordered"] == bool("true"):
        ordered = ordered + 1
    if data["status"]["abandoned"] == "TRUE":
        abandoned = abandoned + 1
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

contributors = pd.Series(contributors)
unique_cont = len(contributors.unique())
production = ordered - (build_ready + abandoned)
failures = (vector + unknown)
not_attempted = (build_ready - attempted)

print("Contributers :", unique_cont)
print("Ordered :", ordered)

print("2017.42 :", sub1)
print("2017.44 :", sub2)
print("2017.46 :", sub3)
print("2017.48 :", sub4)
print("2017.50 :", sub5)

print("Abandoned :", abandoned)
print("Received :", build_ready)
print("In Production : ", production)
print("Build Attempted :", attempted)
print("Verified Success :", complete)
print("Failures :", failures)
print("Not Yet Attempted :", not_attempted)
print()

now = datetime.datetime.now()

submission1 = "2017.42 [{}] Total Ordered #2b5aa5.7".format(sub1)
submission2 = "2017.44 [{}] Total Ordered #2b5aa5.7".format(sub2)
submission3 = "2017.46 [{}] Total Ordered #2b5aa5.7".format(sub3)
submission4 = "2017.48 [{}] Total Ordered #2b5aa5.7".format(sub4)
submission5 = "2017.50 [{}] Total Ordered #2b5aa5.7".format(sub5)
aband = "Total Ordered [{}] Abandoned #ef8c81".format(abandoned)
rec = "Total Ordered [{}] Received #418fba.9".format(build_ready)
prod = "Total Ordered [{}] In Production #f7c862".format(production)
att = "Received [{}] Build Attempted #6e9db7.7".format(attempted)
ver = "Build Attempted [{}] Verified Success #43bc68.7".format(complete)
fail = "Build Attempted [{}] Failures #a8252b.7".format(failures)
#pro = "Build Attempted [{}] In Process  #e8c620.4".format()
not_att = "Received [{}] Not Yet Attempted #e2bec0".format(not_attempted)

sankey = submission1+"\n"+submission2+"\n"+submission3+"\n"+submission4+"\n"+submission5+"\n"+aband+"\n"+rec+"\n"+prod+"\n"+att+"\n"+ver+"\n"+fail+"\n"+not_att
print(sankey)

name = "../raw_files/sankey_diagrams/sankey_{}.txt".format(str(now))
print(name)

f = open(name,'w')
f.write(sankey)
f.close()

#note
