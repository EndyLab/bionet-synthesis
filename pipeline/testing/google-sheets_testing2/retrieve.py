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

## ==========
## Parameters
## ==========

single_submission_id = '1rR4RQ1GLN3eMkOHVRc3jXZ6ZR2FfbUGU9JUXkm8OFHk'
bulk_submission_id = '1UtZkWkogPifDyD9sw46YM0PrSGXCh_2KJyfCdYqkIJs'
ID_next = 0
for file in glob.glob("template.json"): #(BASE_PATH + "/pipeline/testing/json/template.json"):
        with open(file,"r") as template_json:
            template = json.load(template_json)

## =========
## Functions
## =========

def NextCollection():
    data = glob.glob("./data/*/*.json")
    collection_list = [1]
    for json_file in data:
        with open(json_file,"r") as json_data:
            collection_json = json.load(json_data)
            collection_list.append(collection_json["info"]["gene_metadata"]["collection_id"])
    return(int(max([e for e in collection_list if isinstance(e, int)])) + 1)

def NextID(counter):
    data = (sorted(glob.glob("./data/*")))
    string = (data[-1])
    number = int(string[-6:]) + 1 + counter
    string_number = str(number)
    id_number = (string_number.zfill(6))
    full_id = "BBF10K_" + id_number
    return full_id

def csvtext_to_pandas(csvtext, byte_string=False):
    if byte_string == True:
        csv_file= StringIO(str(csvtext, 'utf-8'))
        data = pd.read_csv(csv_file)
        return data
    else: 
        csv_file= StringIO(csvtext)
        data = pd.read_csv(csv_file)
        return data

def extract_google_form(google_id):
    response = requests.get('https://docs.google.com/spreadsheet/ccc?key=' + google_id + '&output=csv')
    assert response.status_code == 200, 'Wrong status code'
    byte_string=(response.content)
    data = csvtext_to_pandas(byte_string, True)
    #csv_file= StringIO(str(byte_string, 'utf-8'))
    #data = pd.read_csv(csv_file)
    return data

def fill_strip(data):
    data = data.fillna('')
    data = strip_df(data)
    return data

#previous = pd.read_csv('previous_submissions.csv')
def uniq_data(current_data, previous):
    current_data = current_data.fillna('')
    current_data = strip_df(current_data)
    previous = previous.fillna('')
    previous = strip_df(previous) 
    united_data = pd.concat([current_data, previous])
    united_data.fillna('')
    united_data_grouped = united_data.groupby(list(united_data.columns))
    uniq_data_idx = [x[0] for x in united_data_grouped.indices.values() if len(x) == 1]
    uniq_data = united_data.iloc[uniq_data_idx]
    return uniq_data

def url_fixer(url):
    link,file_id = re.match(r'(https:\/\/drive.google.com\/)open\?id=([A-Za-z0-9-_]+)',url).groups()
    file_url = link + "uc?id=" + file_id
    return file_url

def get_google_textfile(url, zip_file=False):
    file_url = url_fixer(url)
    if zip_file == True:
        data = requests.get(file_url)
        return data
    else: 
        data = requests.get(file_url)
        data.encoding = 'utf-8'
        data = data.text
        return data

def strip_df(df):
    for column in df:
        new_col = []
        for item in df[column]:
            new_col.append(str(item).strip())
        df[column] = new_col
    return df

## =======
## Classes
## =======
class FreeGene:
    """FreeGene class"""
    def __init__(self, gene_id, collection_id, timestamp, author_name, author_email, author_affiliation, author_orcid, gene_name, who, what, where, why, future_app, database_links, part_type, target_organism, optimize, safety, genbank_file, template_json):
        self.gene_id = gene_id
        self.collection_id = collection_id
        self.submission_timestamp = timestamp
        self.author_name = author_name
        self.author_email = author_email
        self.author_affiliation = author_affiliation
        self.author_orcid = author_orcid
        self.gene_name = gene_name
        self.who = who
        self.what = what
        self.where = where
        self.why = why
        self.future_app = future_app
        self.database_links = database_links
        self.part_type = part_type.lower()
        self.target_organism = target_organism
        self.optimize = optimize
        self.safety = safety
        self.genbank_file = genbank_file
        self.original_sequence = str(SeqIO.read(StringIO(genbank_file), "genbank").seq)
        if part_type == "CDS":
            self.optimized = optimize_sequence.optimize_gene(self.original_sequence, self.target_organism)
        else:
            self.optimized = self.original_sequence
        self.fragments = fragment.fragment_gene(self.optimized,self.part_type)
        if len(self.optimized) < 300:
            self.cloning_enzyme = "BtgZI"
        else:
            self.cloning_enzyme = "BbsI"
        self.retrieval_enzyme = "BsaI"
    def write(self): 
        path = "{}/data/{}".format(".",gene_id)
        os.makedirs(path)
        template["gene_id"] = self.gene_id
        template["author"]["name"] = self.author_name
        template["author"]["email"] = self.author_email
        template["author"]["affiliation"] = self.author_affiliation
        template["author"]["orcid"] = self.author_orcid
        template["info"]["documentation"]["gene_name"] = self.gene_name
        template["info"]["documentation"]["who"] = self.who
        template["info"]["documentation"]["what"] = self.what
        template["info"]["documentation"]["where"] = self.where
        template["info"]["documentation"]["why"] = self.why
        template["info"]["documentation"]["future_app"] = self.future_app
        template["info"]["documentation"]["database_links"] = self.database_links
        template["info"]["gene_metadata"]["cloning"]["part_type"] = self.part_type
        template["info"]["gene_metadata"]["cloning"]["cloning_enzyme"] = self.cloning_enzyme
        template["info"]["gene_metadata"]["cloning"]["retrieval_enzyme"] = self.retrieval_enzyme
        #template["info"]["gene_metadata"]["target_organism"]["organism_name"] = self.target_organism
        template["info"]["gene_metadata"]["safety"] = self.safety
        template["info"]["gene_metadata"]["collection_id"] = self.collection_id
        template["dates"]["submitted"] = self.submission_timestamp
        # Write taxid, target_organism
        template["sequence"]["original_sequence"] = self.original_sequence
        template["sequence"]["optimized_sequence"] = self.optimized
        template["sequence"]["fragment_sequences"] = {}
        for index,frag in enumerate(self.fragments):
            fragment_name = self.gene_id + "_" + str(index + 1)
            template["sequence"]["fragment_sequences"][fragment_name] = frag
            print(frag)
        # Figure out how to write fragments 
        with open("{}/{}.json".format(path,gene_id),"w+") as json_file:
            json.dump(template,json_file,indent=2)
        with open("{}/{}.gb".format(path,gene_id),"w+") as genbank_single:
            genbank_single.write(self.genbank_file)


## ======
## Script
## ======

print ("=== BULK DATA ===")
# Bulk data
current_data_bulk = fill_strip(extract_google_form(bulk_submission_id))
previous_bulk = fill_strip(pd.read_csv('previous_bulk.csv'))
bulk_data = uniq_data(current_data_bulk, previous_bulk)
print ("========== Current_data ===========")
print (current_data_bulk)
print ("========== Previous_bulk ==========")
print (previous_bulk)
print ("========== Uniq_data ==============")
print (bulk_data)
for index, row in bulk_data.iterrows():
    csv_url = row['CSV file']
    csv_data = csvtext_to_pandas(get_google_textfile(csv_url))
    for index_csv, row_csv in csv_data.iterrows():
        # Setup import into object
        gene_id = NextID(ID_next) # New ID
        zip_file = zipfile.ZipFile(io.BytesIO(requests.get(url_fixer(row['Genbank files in zip file'])).content)) # Download zip file
        for name in zip_file.namelist():
            if row_csv['Genbank file'] in name:
                with zip_file.open(name) as myfile:
                    genbank_file="\n".join(str(myfile.read(), 'utf-8').splitlines()) # Format genbank file all nice 
        # Import into object 
        freegene = FreeGene(gene_id, NextCollection(), row['Timestamp'], row['Name (first and last)'], row['Email address'], row['Affiliation'], row['ORCID'], row_csv['Gene name'], row_csv['Who would be interested or who would they be useful to?'], row_csv['What exactly is this gene? What does it do?'], row_csv['Where does this gene come from? How was it discovered?'], row_csv['Why are you synthesizing this gene?'], row_csv['What are some future applications of this gene?'], row_csv['Database links'], row_csv['Part type'], row_csv['Target organism'], row_csv['Optimize?'], row_csv['Safety information'], genbank_file, template)
        freegene.write()
current_data_bulk.to_csv('previous_bulk.csv', index=False)

# Single submissions
current_data_single = fill_strip(extract_google_form(single_submission_id))
previous_single = fill_strip(pd.read_csv('previous_single.csv'))
single_data = uniq_data(current_data_single, previous_single)
print ("========== Current_data ===========")
print (current_data_single)
print ("========== Previous_bulk ==========")
print (previous_single)
print ("========== Uniq_data ==============")
print (single_data)
for index, row in single_data.iterrows():
    gene_id = NextID(ID_next) # New ID
    genbank_file = get_google_textfile(row['Genbank file'])
        # Import into object 
    freegene = FreeGene(gene_id, NextCollection(), row['Timestamp'], row['Name (first and last)'], row['Email Address'], row['Affiliation'], row['ORCID'], row['Gene name'], row['Who would be interested or who would they be useful to?'], row['What exactly is this gene? What does it do?'], row['Where does this gene come from? How was it discovered?'], row['Why are you synthesizing this gene?'], row['What are some future applications of this gene?'], row['Database links'], row['Part type'], row['Target organism'], row['Optimize?'], row['Safety information'], genbank_file, template)
    freegene.write()
current_data_single.to_csv('previous_single.csv', index=False)
sys.exit()




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
        # template["status"]["will_build"] = True
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
        unknown_enzyme = []
        template["sequence"]["fragment_sequences"] = {}
        fragments,enzyme = fragment.fragment_gene(str(gb_record.seq),row["Part type"])
        template["info"]["type"]["build_type"] = enzyme

        # Makes sure that the small fragment won't get ordered until paired
        if enzyme == "BtgZI":
            template["status"]["will_build"] = False
        elif enzyme == "BbsI":
            template["status"]["will_build"] = True
        else:
            unknown_enzyme.append(id_num_str)

        # Fragments are named based on the gene id and the fragment number
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

    # Remove the temporary genbank file
    os.remove(gb_name)

## ====================================================
## Query the database for all of the small sequences
## and sequences to attach them to
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





