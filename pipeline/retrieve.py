from Bio import SeqIO

import requests

import json
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas as pd
import re
import os
import glob
import datetime
import shutil

import fragment
import optimize

from config import *

def strip_df(df):
    for column in df:
        new_col = []
        for item in df[column]:
            new_col.append(str(item).strip())
        df[column] = new_col
    return df

# Takes in the template json file to initialize the json files
for file in glob.glob(BASE_PATH + "/pipeline/testing/json/template.json"):
        with open(file,"r") as template_json:
            template = json.load(template_json)

# Find the next id number
previous_genes = (sorted(glob.glob(BASE_PATH + "/data/*")))
string = (previous_genes[-1])
id_num = int(string[-6:])
print("Last ID number: ",id_num)

# Use creds to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name('Free Genes-153172d39301.json', scope)
client = gspread.authorize(creds)

# Find a workbook by name and open the first sheet
sheet = client.open("Free Genes Project Rolling Submissions (Responses)").sheet1

# Extract and print all of the values
current_data = pd.DataFrame(sheet.get_all_records())
current_data = current_data.fillna('')
current_data = strip_df(current_data)
# print("======================================\ndata\n=================================\n",current_data)
previous = pd.read_csv(BASE_PATH + '/raw_files/previous_submissions/previous_submissions.csv')
previous = previous.fillna('')
previous = strip_df(previous)
# print("======================================\nprevious\n=================================\n",previous)

# Extract unique current_data
united_data = pd.concat([current_data, previous])
# print("united w/ na\n",united_data)
united_data.fillna('')
# print("united\n",united_data)
united_data_grouped = united_data.groupby(list(united_data.columns))
uniq_data_idx = [x[0] for x in united_data_grouped.indices.values() if len(x) == 1]
uniq_data = united_data.iloc[uniq_data_idx]
print("===============\nunique\n===============\n",uniq_data)
uniq_data.to_csv("./test_unique.csv")
input("pause")
# Iterate over each line in the Google Sheets
for index, row in uniq_data.iterrows():

    # Import the genbank sequence as a text file and rewrite it as a gb file
    genbank_url = row['Genbank file']
    link,file_id = re.match(
        r'(https:\/\/drive.google.com\/)open\?id=([A-Za-z0-9-_]+)',
        genbank_url).groups()
    file_url = link + "uc?id=" + file_id
    # input("use URL {}".format(file_url))
    gb_data = requests.get(file_url).text
    gb_name = "temp.gb"
    print("working on {}".format(row["Collection name"]))
    with open(gb_name,"w+") as genbank:
        genbank.write(gb_data)

    for gb_record in SeqIO.parse(open(gb_name,"r+"), "genbank"):
        template["gene_name"] = gb_record.description
        print("Current Part: ",template["gene_name"])

        ## ====================================================
        ## Pull in information from submission and genbank
        ## ====================================================
        id_num += 1
        id_num_str = "BBF10K_" + str(id_num).zfill(6)

        # Pull out sequence specific information
        template["gene_id"] = id_num_str
        template["info"]["type"]["part_type"] = row["Part type"]
        template["project_description"] = row["Detailed description"]
        template["sequence"]["original_sequence"] = str(gb_record.seq)

        # Extract information from Sheet and populate the template json file
        template["author"]["name"] = row["Name (first and last)"]
        template["author"]["email"] = row["Email Address"]
        template["author"]["affiliation"] = row["Affiliation"]
        template["author"]["orcid"] = row["ORCID "]
        template["author"]["project"] = row["Project"]
        template["info"]["safety"] = row["Safety information"]
        template["info"]["other_tags"] = row["Other tags"]
        template["info"]["database_links"] = row["Database links"]

        # Hard code certain parameters
        template["status"]["will_build"] = True
        template["info"]["type"]["cloning_method"] = "10K_MoClo"
        template["status"]["ordered"] = False

        # Add in timestamps
        timestamp = row["Timestamp"]
        date,time = timestamp.split(" ")
        month,day,year = date.split("/")
        week = datetime.date(int(year), int(month), int(day)).isocalendar()[1]
        order_week = str(year) + "." + str(week)
        template["info"]["order_week"] = order_week

        ## ====================================================
        ## Optimize sequence
        ## ====================================================
        if template["info"]["type"]["part_type"] == "cds":
            print("original\n",str(gb_record.seq))
            print("optimized\n",optimize.optimize_gene(str(gb_record.seq),row["Codon usage "]))
            template["sequence"]["optimized_sequence"] = optimize.optimize_gene(str(gb_record.seq),row["Codon usage "])
        else:
            print("not a cds")
            template["sequence"]["optimized_sequence"] = str(gb_record.seq)

        ## ====================================================
        ## Fragment sequence
        ## ====================================================
        template["sequence"]["fragment_sequences"] = {}
        fragments,enzyme = fragment.fragment_gene(str(gb_record.seq),row["Part type"])
        template["info"]["type"]["build_type"] = enzyme
        for index,frag in enumerate(fragments):
            fragment_name = id_num_str + "_" + str(index + 1)
            template["sequence"]["fragment_sequences"][fragment_name] = frag

        ## ====================================================
        ## Generate directory and desired files for each part
        ## ====================================================
        path = "{}/data/{}".format(BASE_PATH,id_num_str)
        os.makedirs(path)
        print("Making directory for {}".format(id_num_str))
        with open("{}/{}.json".format(path,id_num_str),"w+") as json_file:
            json.dump(template,json_file,indent=2)
        with open("{}/{}.gb".format(path,id_num_str),"w+") as genbank_single:
            SeqIO.write(gb_record, genbank_single, "genbank")

    # Remove the temporary file
    os.remove(gb_name)

## ====================================================
## Query the database for all of the small sequences and sequences to attach them to
## ====================================================
small_seq_ids = []
small_seqs = []
large_seq_ids = []
large_seqs = []
for file in glob.glob(BASE_PATH + "/data/*/*.json"):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data["status"]["ordered"]:
        continue
    if data["info"]["type"]["build_type"] == "BtgZI":
        small_seq_ids.append(data["gene_id"])
        small_seqs.append(data["sequence"]["fragment_sequences"]["{}_1".format(data["gene_id"])])
    elif data["info"]["type"]["build_type"] == "BbsI":
        print("Num frags: ",len(data["sequence"]["fragment_sequences"]))
        if len(data["sequence"]["fragment_sequences"]) > 1:
            print("too many frags")
            continue
        large_seq_ids.append(data["gene_id"])
        large_seqs.append(data["sequence"]["fragment_sequences"]["{}_1".format(data["gene_id"])])

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

joined_seqs = []
joined_ids = []
fragment_names = []

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

for index,row in joined_df.iterrows():
    with open("{}/data/{}/{}.json".format(BASE_PATH,row["Gene ID"],row["Gene ID"]),"r") as json_file:
        data = json.load(json_file)
    data["sequence"]["fragment_sequences"] = {}
    data["sequence"]["fragment_sequences"][row["Fragment Name"]] = row["Sequence"]
    with open("{}/data/{}/{}.json".format(BASE_PATH,row["Gene ID"],row["Gene ID"]),"w+") as json_file:
        json.dump(data,json_file,indent=2)

will_order = []
will_order_seqs = []

for file in glob.glob("{}/data/*/*.json".format(BASE_PATH)):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data["status"]["ordered"]:
        continue
    if data["info"]["type"]["build_type"] == "BbsI":
        for fragment in data["sequence"]["fragment_sequences"]:
            print("fragment",fragment)
            will_order.append(fragment)
            will_order_seqs.append(data["sequence"]["fragment_sequences"][fragment])
    data["status"]["ordered"] = True
    with open(file,"w+") as json_file:
        json.dump(data,json_file,indent=2)

# Output DNA in Twist order format
twist_dna = pd.DataFrame({
        'gene name': will_order,
        'FASTA_seq': will_order_seqs,
        }, columns=['gene name','FASTA_seq'])

print(twist_dna)

current_data.to_csv(BASE_PATH + '/raw_files/previous_submissions/previous_submissions.csv',index=False)

previous_submissions = (sorted(glob.glob(BASE_PATH + "/submissions/*.csv")))
string = (previous_submissions[-1])
sub_num = int(string[-6:-4])
print("Last ID number: ",sub_num)
twist_dna.to_csv('{}/submissions/submission{}.csv'.format(BASE_PATH,str(sub_num)),index=False)





#
