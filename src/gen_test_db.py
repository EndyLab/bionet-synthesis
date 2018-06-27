import pandas as pd
from config import *
from db_config import *
session,engine = connect_db()

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

from datetime import datetime
import json
import glob
import shutil
import os

data = pd.read_csv('./parts_frags.csv')
counter = 0
for part,df in data.groupby('part_id'):
    counter += 1
    print(counter,part)
    new_part = Part(part_id=df.part_id.values[0],
        part_name=df.part_name.values[0],
        part_type=df.part_type.values[0],
        cloning_enzyme=df.cloning_enzyme.values[0],
        seq=df.seq.values[0],
        status='ordered')

    for i,row in df.iterrows():
        print(row.fragment_name)
        new_frag = Fragment(
            fragment_name=row.fragment_name,
            seq=row['seq.1'].upper(),
            has_part='True')
        session.add(new_frag)
        new_part.fragments.append(new_frag)

    session.add(new_part)

session.commit()
















#
