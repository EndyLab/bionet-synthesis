import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import re
import pandas as pd
from config import *
import db_config
session,engine = db_config.connect_db()

j_counter = 0
linked = []
linked_seq = []
no_frag = []
old_parts = [part.part_id for part in session.query(Part)]
for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    print(data['gene_id'])
    if data['gene_id'] in old_parts:
        print('Found old part')
        continue
    new = Part(part_id=data['gene_id'],
        part_name=data['info']['documentation']['gene_name'],
        part_type=data['info']['gene_metadata']['cloning']['part_type'],
        cloning_enzyme=data['info']['gene_metadata']['cloning']['cloning_enzyme'],
        seq=data['sequence']['optimized_sequence'].upper(),
        status='ordered')
    frags = []
    if data['sequence']['fragment_sequences'] == {}:
        no_frag.append(data['gene_id'])
    elif data['info']['gene_metadata']['cloning']['cloning_enzyme'] == 'BtgZI':
        linked.append(new)
        frag = [frag for frag in data['sequence']['fragment_sequences']][0]
        linked_seq.append(data['sequence']['fragment_sequences'][frag])
    else:
        for fragment in data['sequence']['fragment_sequences']:
            new_frag = Fragment(
                fragment_name=fragment,
                seq=data['sequence']['fragment_sequences'][fragment].upper(),
                parts=[new],cloning_method='entry_1')
            session.add(new_frag)
    session.add(new)
## Add the extra connections in for the small linked genes
for part,seq in zip(linked,linked_seq):
    frag = session.query(Fragment).filter(Fragment.seq == seq).one()
    part.fragments.append(frag)
    session.add(part)

print(no_frag)


commit = int(input("Commit changes (1-yes, 2-no): "))
if commit == 1:
    session.commit()





    #
