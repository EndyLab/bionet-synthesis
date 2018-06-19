from datetime import datetime
import pandas as pd
from config import *
import ot_functions as ot
from db_config import *
session,engine = connect_db()

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import json
import glob
import os

author_dict = []
for file in sorted(glob.glob('../data/*/*.json')):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    author_dict.append([data['gene_id'],data['author']['name']])
authors = [a for g,a in author_dict]
author_dict = dict(author_dict)

authors = pd.Series(authors).unique()
print('Current authors:\n',authors)
author = ot.request_info('Enter author name: ',select_from=authors)

query = "SELECT parts.part_id,parts.status,builds.build_name,part_wells.address,part_wells.seq_outcome,fragments.fragment_name,frag_plates.plate_name,frag_wells.address FROM parts \
        INNER JOIN wells AS part_wells ON parts.id = part_wells.part_id\
        INNER JOIN plates AS part_plates ON part_wells.plate_id = part_plates.id\
        INNER JOIN builds ON part_plates.build_id = builds.id\
        INNER JOIN part_frag ON parts.id = part_frag.part_id\
        INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
        INNER JOIN wells AS frag_wells ON fragments.id = frag_wells.fragment_id\
        INNER JOIN plates AS frag_plates ON frag_wells.plate_id = frag_plates.id\
        WHERE part_wells.plate_type = 'seq_plate'"

data = pd.read_sql_query(query,con=engine)
data['author'] = data.part_id.apply(lambda x: author_dict[x])
data = data.sort_values('part_id')
data = data.sort_values('build_name')

data = data[data.author == author]

now = str(datetime.now()).split(" ")[0]
author_path = '{}/authors/{}'.format(BASE_PATH,author)
ot.make_directory(author_path)
data.to_csv('{}/{}_well_locations_{}'.format(author_path,author,now),index=False)
print(data)








#
