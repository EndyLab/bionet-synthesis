from Bio import SeqIO

import requests

import json
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas
import re
import os
import glob
import datetime
import shutil

# from config import *

# Takes in the template json file to initialize the json files
for file in glob.glob("../json/template.json"):
        print(file)
        with open(file,"r") as template_json:
            template = json.load(template_json)

# Find the next id number and make a list of previously entered genes
id_num = 0
previous_genes = []
for file in glob.glob("../../../data/*/*.json"):
    nothing,current_id_num = re.match(
        r'.*/([A-Za-z0-9]+)_([0-9]+).json',
        file).groups()
    current_id_num = int(current_id_num)
    if current_id_num > id_num:
        id_num = current_id_num
    with open(file,"r") as json_file:
        data = json.load(json_file)
    previous_genes.append(data["gene_name"])
print(id_num)
print(len(previous_genes))

# use creds to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name('Free Genes-153172d39301.json', scope)
client = gspread.authorize(creds)

# Find a workbook by name and open the first sheet
sheet = client.open("Free Genes Project Rolling Submissions (Responses)").sheet1

# Extract and print all of the values
data = pandas.DataFrame(sheet.get_all_records())

# Iterate over each line in the Google Sheets
for index, row in data.iterrows():

    # Import the genbank sequence as a text file and rewrite it as a gb file
    genbank_url = row['Genbank file']
    link,file_id = re.match(
        r'(https:\/\/drive.google.com\/)open\?id=([A-Za-z0-9]+)',
        genbank_url).groups()
    file_url = link + "uc?id=" + file_id
    gb_data = requests.get(file_url).text
    gb_name = "temp.gb"
    with open(gb_name,"w+") as genbank:
        genbank.write(gb_data)

    # Parse the
    for gb_record in SeqIO.parse(open(gb_name,"r+"), "genbank"):

        # Checks to see if the gene has already been entered
        if gb_record.name in previous_genes:
            print("Already entered {}".format(gb_record.name))
            continue

        # Assign unique id numbers to each entry
        id_num += 1
        id_num_str = "BBF10K_" + str(id_num).zfill(6)

        # Pull out sequence specific information
        template["gene_id"] = id_num_str
        template["gene_name"] = gb_record.definition
        print("name",template["gene_name"])
        template["info"]["type"]["part_type"] = row["Part type"]
        template["project_description"] = row["Detailed description"]
        template["part_description"] = gb_record.description
        print(gb_record.description)
        template["sequence"]["original_sequence"] = str(gb_record.seq)

        # REVIEW: change to account for different part types
        template["info"]["type"]["build_type"] = "10K_MoClo-EntryCDS-BbsI"


        # Extract information from Sheet and populate the template json file
        template["author"]["name"] = row["Name (first and last)"]
        template["author"]["email"] = row["Email Address"]
        template["author"]["affiliation"] = row["Affiliation"]
        template["author"]["orcid"] = row["ORCID"]
        template["author"]["project"] = row["Project"]
        template["info"]["safety"] = row["Safety information"]
        template["info"]["other_tags"] = row["Other tags"]
        template["info"]["database_links"] = row["Database links"]

        # Hard code certain parameters
        template["status"]["will_build"] = True
        template["info"]["type"]["cloning_method"] = "10K_MoClo"

        # Add in timestamps
        timestamp = row["Timestamp"]
        date,time = timestamp.split(" ")
        month,day,year = date.split("/")
        week = datetime.date(int(year), int(month), int(day)).isocalendar()[1]
        order_week = str(year) + "." + str(week)
        template["info"]["order_week"] = order_week

        #State the path to house the new set of directories
        path = "../../../data/{}".format(id_num_str)
        os.makedirs(path)
        print("Making directory for {}".format(id_num_str))
        with open("{}/{}.json".format(path,id_num_str),"w+") as json_file:
            json.dump(template,json_file,indent=2)
        with open("{}/{}.gb".format(path,id_num_str),"w+") as genbank_single:
            SeqIO.write(gb_record, genbank_single, "genbank")

    # Remove the temporary file
    os.remove(gb_name)







#
