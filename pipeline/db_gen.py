import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import pandas as pd
from config import *
from db_config import *

def build_db():
    print("\n============ Building the database ============\n")
    no_frag = []
    no_frag_loc = []
    linked = []
    linked_seq = []
    names = []
    id_nums = []

    for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
        with open(file,"r") as json_file:
            data = json.load(json_file)
        names.append(data['info']['documentation']['gene_name'])
        id_nums.append(data['gene_id'])
    dictionary = dict(zip(names,id_nums))

    ## Generate part and fragment objects from JSON database
    j_counter = 0
    for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
        with open(file,"r") as json_file:
            data = json.load(json_file)
        new = Part(part_id=data['gene_id'],
            part_name=data['info']['documentation']['gene_name'],
            part_type=data['info']['gene_metadata']['cloning']['part_type'],
            seq=data['sequence']['optimized_sequence'],
            status='ordered')
        frags = []
        if data['sequence']['fragment_sequences'] == {}:
            no_frag.append(data['gene_id'])
        elif data['info']['gene_metadata']['cloning']['cloning_enzyme'] == 'BtgZI':
            linked.append(new)
            linked_seq.append(data['sequence']['fragment_sequences'][data['gene_id']+"_1"])
        else:
            for fragment in data['sequence']['fragment_sequences']:
                if data['location']['fragments'] == {}:
                    no_frag_loc.append(data['gene_id'])
                    session.add(Fragment(
                        fragment_name=fragment,
                        seq=data['sequence']['fragment_sequences'][fragment],
                        parts=[new]))
                else:
                    session.add(Fragment(
                        fragment_name=fragment,
                        location=data['location']['fragments'][fragment],
                        seq=data['sequence']['fragment_sequences'][fragment],
                        parts=[new]))
        session.add(new)

    for part,seq in zip(linked,linked_seq):
        frag = session.query(Fragment).filter(Fragment.seq == seq).first()
        part.fragments.append(frag)
        session.add(part)
    print('\nCreated database\n')
    return session
