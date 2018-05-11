import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime

from config import *
from db_config import *

confirmed = [part.part_id for part in session.query(Part).filter(Part.status == 'sequence_confirmed').order_by(Part.id)]
# print(confirmed)
author_info = []
ignore = ['Stanford BIOE80 class 2017','Stanley Qi','Keoni Gandall','Scott Pownall']
for part in confirmed:
    for file in glob.glob('../data/{}/{}.json'.format(part,part)):
#         print(file)
        with open(file,"r") as json_file:
            data = json.load(json_file)
        if data['author']['name'] in ignore:
            continue
        if data['version'] == '2.0.0':
            author_info.append([data['author']['name'],data['author']['email'],data['author']['affiliation'],data['gene_id'],data['info']['documentation']['gene_name'],data['info']['documentation']['what']])
        elif data['version'] == '3.0.0':
            author_info.append([data['author']['name'],data['author']['email'],data['author']['affiliation'],data['gene_id'],data['info']['documentation']['gene_name'],data['info']['documentation']['description']])

author_info = np.array(author_info)
authors = pd.DataFrame({
    'Name' : author_info[:,0],
    'Email' : author_info[:,1],
    'Affiliation' : author_info[:,2],
    'Gene ID' : author_info[:,3],
    'Gene Name' : author_info[:,4],
    'Description' : author_info[:,5]
})

authors = authors[['Name','Email','Affiliation','Gene ID','Gene Name','Description']]
authors.to_csv('./author_info.csv')

name_email = authors[['Name','Email']].drop_duplicates()
names = pd.Series(authors['Name'])
unique_names = names.unique()
print("Current authors:\n",names.value_counts())
current = str(datetime.now()).split(" ")[0]

genes = []
for file in sorted(glob.glob('../data/*/*.json')):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data['author']['name'] not in unique_names:
        continue
    genes.append([data['author']['name'],data['gene_id']])

genes = np.array(genes)
df_genes = pd.DataFrame({
    'Author' : genes[:,0],
    'Gene ID' : genes[:,1]
})

for name in unique_names:
    temp = df_genes[df_genes['Author'] == name]
    new_rows = []
    for i,row in temp.iterrows():
        gene = row['Gene ID']
        print("Gene",gene)
        part = session.query(Part).filter(Part.part_id == gene).one()
        attempts = [well.seq_outcome for well in part.wells if well.plates.plate_type == 'seq_plate']
        if len(attempts) == 0:
            attempts += ['N/A','N/A']
        elif len(attempts) == 1:
            attempts += ['N/A']
        print(part.part_name,part.part_id,part.status,attempts[0],attempts[1])
        new_rows.append([part.part_name,part.part_id,part.status,attempts[0],attempts[1]])
    new_rows = np.array(new_rows)
    new_df = pd.DataFrame({
        "Gene Name" : new_rows[:,0],
        "Gene ID" : new_rows[:,1],
        "Current Status" : new_rows[:,2],
        "Attempt 1" : new_rows[:,3],
        "Attempt 2" : new_rows[:,4],
    })
    new_df = new_df[['Gene Name','Gene ID','Current Status','Attempt 1','Attempt 2']]
    file_name = "{}_status_{}.csv".format(name,current)
    path = '{}/authors/{}'.format(BASE_PATH,name)
    if os.path.exists(path):
        print("Directory for {} already exists".format(name))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)
        print("Making directory for {}".format(name))
    file_path = path + "/" + file_name
    new_df.to_csv(file_path,index=False)
