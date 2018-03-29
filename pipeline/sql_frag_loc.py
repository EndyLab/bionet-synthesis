import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import pandas as pd
from config import *
from db_config import *

# ## Generate the SQL Database from json files and plate maps
# session = build_db()

def frag_assign(session):
    print("\n======= Beginning fragment location assignment =======\n")
    ## Create a dictionary to including all of the gene names and their ID numbers
    ## To check for whether the fragments can link to a buildable part
    names = []
    id_nums = []
    for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
        with open(file,"r") as json_file:
            data = json.load(json_file)
        names.append(data['info']['documentation']['gene_name'])
        id_nums.append(data['gene_id'])
    dictionary = dict(zip(names,id_nums))

    ## Build a list of the previously entered synthesis plates
    plates_made = [plate.plate_name for plate in session.query(Plate)]
    print("Plates that have already been parsed:\n",plates_made)

    ## Take in plate maps from twist and generate fragment objects
    not_found = []
    for file in glob.glob('{}/plate_maps/*.csv'.format(BASE_PATH)):
        plate_map = pd.read_csv(file)
        plate_map['customer_line_item_id'] = plate_map['customer_line_item_id'].str.strip()
        plates_made = []
        for index,row in plate_map.iterrows():
            if row['customer_line_item_id'][:-2] not in names and row['customer_line_item_id'][:-2] not in id_nums:
                not_found.append(row['customer_line_item_id'])
                continue
            if row['Plate'] not in plates_made:
                current_plate = Plate('','syn_plate',plate_name=row['Plate'])
                plates_made.append(row['Plate'])
            current_frag = ''
            current_frag = session.query(Fragment).filter(Fragment.seq == row['Insert Sequence']).first()
            if current_frag:
                current_plate.add_item(current_frag,row['Well'],syn_yield=row['Yield (ng)'])
                session.add(current_plate)

    ## Determie which samples are build ready
    for part in session.query(Part):
        no_build = False
        for frag in part.fragments:
            if not frag.wells:
                no_build = True
        if not no_build:
            part.change_status('received')
        session.add(part)

    # for part in session.query(Part).filter(Part.part_id.like('BBF10K_0000%')).order_by(Part.part_id):
    #     print(part.part_id,part.status)

    return session







#
