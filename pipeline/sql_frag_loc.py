import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import pandas as pd
from config import *
from db_config import *

def frag_assign():
    ## Build a list of the previously entered synthesis plates
    plates_made = [plate.plate_name for plate in session.query(Plate).filter(Plate.plate_type == 'syn_plate')]
    print("Plates that have already been parsed:\n",plates_made)

    ## Take in the previous list of fragments that weren't placed
    if glob.glob('{}/raw_files/not_found.json'.format(BASE_PATH)):
        with open('{}/raw_files/not_found.json'.format(BASE_PATH),'r') as json_file:
            not_found = json.load(json_file)
    else:
        not_found = {}

    ## Take in plate maps from twist and generate fragment objects
    for file in glob.glob('{}/plate_maps/*.csv'.format(BASE_PATH)):
        plate_map = pd.read_csv(file)
        ## Have to account for two different csv formats
        try:
            plate_map['Plate']
            version = 1
            plate_key = 'Plate'
            well_key = 'Well'
            name_key = 'customer_line_item_id'
            sequence_key = 'Insert Sequence'
        except:
            try:
                plate_map['Plate ID']
                version = 2
                plate_key = 'Plate ID'
                well_key = 'Well Location'
                name_key = 'Name'
                sequence_key = 'Insert sequence'
            except:
                print("Doesn't fit the standard formats")
        unique = plate_map[plate_key].unique()
        for plate in unique:
            # Create a subset of rows pertaining to a specific plate
            current_plate_map = plate_map[plate_map[plate_key] == plate]
            current_plate = Plate('','syn_plate',plate_name=plate)
            session.add(current_plate)
            current_missing = []
            for index,row in current_plate_map.iterrows():
                current_frag = ''
                current_frag = session.query(Fragment).filter(Fragment.seq == row[sequence_key]).first()
                # Checks if a corresponding fragment is found and then adds it to the plate
                if current_frag:
                    current_plate.add_item(current_frag,address=row[well_key].strip(),syn_yield=row['Yield (ng)'])
                else:
                    current_missing.append(row[name_key].strip())
                    print(row[name_key].strip())
                    print(row[sequence_key])
                    print(len(row[sequence_key]))
            current_missing = {plate : current_missing}
            print(current_missing)
            not_found.update(current_missing)

    ## Write the new fragments that were not placed to a file
    with open('{}/raw_files/not_found.json'.format(BASE_PATH),'w+') as json_file:
        json.dump(not_found,json_file,indent=2)
    print(not_found)    ## Determie which samples are build ready
    for part in session.query(Part):
        no_build = False
        for frag in part.fragments:
            if not frag.wells:
                no_build = True
        if not no_build:
            part.change_status('received')
        session.add(part)

    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()
    return

if __name__ == "__main__":
    frag_assign()




#
