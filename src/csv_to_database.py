import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import re
import pandas as pd
from config import *
from db_config import *
import ot_functions as ot

session,engine = connect_db()

query_parts = 'SELECT * FROM parts'

data = pd.read_sql_query(query_parts,con=engine)

print('Currently there are {} parts in the database'.format(len(data)))

if len(data) == 0:
    prefix = ot.request_info('Enter the desired prefix for the unique id numbers (i.e. "initials"_)')
else:
    prefix = data.part_id.tolist()[0].split('_')[0] + '_'

start_id = len(data) + 1

files = []
for i,file in enumerate(glob.glob('{}/src/*.csv'.format(BASE_PATH))):
    print('{}: {}'.format(i,file))
    files.append(file)

ans = ot.request_info('Which file would you like to upload: ',type='int',select_from=[num for num in range(len(files))])
target_file = files[ans]
data = pd.read_csv(target_file)

ot.print_center('...Adding samples to the database...')

for i,part in data.iterrows():
    current_id = prefix + str(start_id + i).zfill(3)
    new_part = Part(
        part_id=current_id,
        part_name=part['part_name'],
        part_type=part['part_type'],
        original_seq=part['original_sequence'].upper()
    )
    session.add(new_part)



session.commit()
