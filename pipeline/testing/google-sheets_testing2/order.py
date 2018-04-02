from Bio import SeqIO
from io import StringIO
import requests
import json
import pandas as pd
import re
import os
import glob
import datetime
import shutil 
import sys
import fragment
import optimize as optimize_sequence
from zipfile import ZipFile
from io import BytesIO
import requests
import zipfile
import io

import optimize
#from config import *

## ====================================================
## Query the database for all of the small sequences
## and sequences to attach them to
## ====================================================
small_seq_ids = []
small_seqs = []
large_seq_ids = []
large_seqs = []
for file in glob.glob("./data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data["status"]["ordered"]:
        continue
    if data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BtgZI":
        small_seq_ids.append(data["gene_id"])
        small_seqs.append(data["sequence"]["fragment_sequences"]["{}_1".format(data["gene_id"])])
    elif data["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] == "BbsI":
        print("Num frags: ",len(data["sequence"]["fragment_sequences"]))
        if len(data["sequence"]["fragment_sequences"]) > 1:
            print("too many frags")
            continue
        large_seq_ids.append(data["gene_id"])
        large_seqs.append(data["sequence"]["fragment_sequences"]["{}_1".format(data["gene_id"])])

# Generate dataframes that are sorted in opposite directions based on length
# which pairs the smallest large fragment with the largest small fragment
small_df = pd.DataFrame({
    "Gene ID" : small_seq_ids,
    "Sequence" : small_seqs,
    "Length" : [len(seq) for seq in small_seqs]
})
small_df = small_df.sort_values("Length",ascending=False)
large_df = pd.DataFrame({
    "Gene ID" : large_seq_ids,
    "Sequence" : large_seqs,
    "Length" : [len(seq) for seq in large_seqs]
})
large_df = large_df.sort_values("Length")

small_counter = 0
print("Total small sequences: ",len(small_df))

## ====================================================
## Join Fragments
## ====================================================
joined_seqs = []
joined_ids = []
fragment_names = []

# Pair sequences together until it runs out of either type of sequence
for index,row in large_df.iterrows():
    print("small counter: ",small_counter)
    if len(small_df) == small_counter:
        print("ran out of small")
        break
    small_row = small_df.iloc[small_counter]
    joined_seq = row["Sequence"] + small_row["Sequence"]
    joined_ids.append(row["Gene ID"])
    joined_seqs.append(joined_seq)
    fragment_names.append(row["Gene ID"] + "_link_" + small_row["Gene ID"])
    joined_ids.append(small_row["Gene ID"])
    joined_seqs.append(joined_seq)
    fragment_names.append(row["Gene ID"] + "_link_" + small_row["Gene ID"])
    small_counter += 1
joined_df = pd.DataFrame({
    "Gene ID" : joined_ids,
    "Sequence" : joined_seqs,
    "Fragment Name" : fragment_names
})

# Change the files in the database to reflect the joined sequences
for index,row in joined_df.iterrows():
    with open("{}/data/{}/{}.json".format(BASE_PATH,row["Gene ID"],row["Gene ID"]),"r") as json_file:
        data = json.load(json_file)
    template["status"]["will_build"] = True
    data["sequence"]["fragment_sequences"] = {}
    data["sequence"]["fragment_sequences"][row["Fragment Name"]] = row["Sequence"]
    with open("{}/data/{}/{}.json".format(BASE_PATH,row["Gene ID"],row["Gene ID"]),"w+") as json_file:
        json.dump(data,json_file,indent=2)

## Find all of the sequences that have yet to be ordered
will_order = []
will_order_seqs = []
for file in glob.glob("{}/data/*/*.json".format(BASE_PATH)):
    with open(file,"r") as json_file:
        data = json.load(json_file)

    # Excludes sequences that have already been ordered and small sequences
    # that haven't been paired yet
    if data["status"]["ordered"] or not template["status"]["will_build"]:
        continue
    # Only pulls the sequence to order from the large fragment
    if data["info"]["type"]["build_type"] == "BbsI":
        for fragment in data["sequence"]["fragment_sequences"]:
            print("fragment",fragment)
            will_order.append(fragment)
            will_order_seqs.append(data["sequence"]["fragment_sequences"][fragment])
    data["status"]["ordered"] = True
    data["info"]["order_number"] = next_sub_num
    with open(file,"w+") as json_file:
        json.dump(data,json_file,indent=2)

# Output DNA in Twist order format
twist_dna = pd.DataFrame({
        'gene name': will_order,
        'FASTA_seq': will_order_seqs,
        }, columns=['gene name','FASTA_seq'])

# Update the previous submissions file so that the submissions don't get logged again
current_data.to_csv(BASE_PATH + '/raw_files/previous_submissions/previous_submissions.csv',index=False)

previous_submissions = (sorted(glob.glob(BASE_PATH + "/submissions/*.csv")))
string = (previous_submissions[-1])
sub_num = int(string[-6:-4])
print("Last ID number: ",sub_num)
twist_dna.to_csv('{}/submissions/submission{}.csv'.format(BASE_PATH,str(sub_num)),index=False)





