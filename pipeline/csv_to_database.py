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

print(len(data))

prefix = 'FG_'
start_id = len(data) + 1

data = pd.read_csv('{}/pipeline/sample_data.csv'.format(BASE_PATH))

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
