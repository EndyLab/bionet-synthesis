from config import *
from db_config import *

import json
import glob
import numpy as np
import pandas as pd

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

# prev_builds = ['build000','build003','build004','build005','build006','build007']

# for build in session.query(Build).order_by(Build.id):
#     # if build.build_name in prev_builds:
#     #     build.status = 'complete'
#     for plate in build.plates:
#         print(build.build_name,build.status,plate.plate_name)

# session.commit()

# for frag in session.query(Fragment).order_by(Fragment.id):
#     print(frag.fragment_name)
#     if "_link_" in frag.fragment_name:
#         print("Linked")
#         continue
#     elif int(frag.fragment_name[-1]) > 1:
#         print("Multi")
#         continue
#     if frag.seq[8:12] == "GGAG":
#         frag.cloning_method = 'entry_1'
#     elif frag.seq[8:12] == "AATG":
#         frag.cloning_method = 'entry_2'
#     else:
#         print(frag.seq)
#         print(frag.seq[8:12])
#         input("Doesn't fit either method")
not_synced = []
not_sql = []
for file in sorted(glob.glob('{}/data/*/*.json'.format(BASE_PATH))):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    print(data['gene_id'],data['version'])
    bionet_id = ''
    try:
        bionet_id = data['virtual_id']
    except:
        print('not synced')
        not_synced.append(data['gene_id'])
    if bionet_id != '':
        try:
            part = session.query(Part).filter(Part.part_id == data['gene_id']).one()
        except:
            print("not in sql")
            not_sql.append(data['gene_id'])
        print('before',part.part_id,part.bionet_id)
        part.bionet_id = bionet_id
        print('after',part.part_id,part.bionet_id)

print(not_synced)
print(not_sql)
input('check not synced')

commit = int(input("Commit changes (1-yes, 2-no): "))
if commit == 1:
    session.commit()

table = pd.read_sql_query('select parts.part_id,parts.bionet_id from parts', con=engine)
print(table)
# row = []
# counter = 0
# for file in sorted(glob.glob('{}/data/*/*.json'.format(BASE_PATH))):
#     with open(file,"r") as json_file:
#         data = json.load(json_file)
#     print([data['gene_id'],data['virtual_id']])
#     if counter == 10:
#         break
#     counter += 1
#     row.append([data['gene_id'],data['virtual_id']])
#
# # array = np.array(row)
# id_dict = dict(row)
# print(id_dict)

# table = pd.read_sql_query('select parts.part_id,parts.bionet_id from parts', con=engine)


#
