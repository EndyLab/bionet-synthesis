import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import sys
import glob
import json
import numpy as np
import pandas as pd
from config import *
from db_config import *
import ot_functions as ot

def frag_assign():
    ## Build a list of the previously entered synthesis plates
    plates_made = [plate.plate_name for plate in session.query(Plate).filter(Plate.plate_type == 'syn_plate')]
    print("Plates that have already been parsed:\n",plates_made)

    ## Stores the last used plate ID so that it can increment
    previous_plates = [plate for plate in session.query(Plate).filter(Plate.plate_type == 'syn_plate')\
        .order_by(Plate.plate_id)]
    if len(previous_plates) == 0:
        print('No previous plate ID, starting at 001')
        current_id = 1
    else:
        last_plate = previous_plates[-1]
        print('last ID ',last_plate.plate_id)
        current_id = int(last_plate.plate_id[-3:])+1

    ## Take in plate maps from twist and generate fragment objects
    for file in glob.glob('{}/src/data/plate_maps/*.csv'.format(BASE_PATH)):
        plate_map = pd.read_csv(file)
        ## Have to account for two different csv formats
        try:
            plate_map['Plate']
            version = 1
            print('version',version)
            plate_key = 'Plate'
            well_key = 'Well'
            name_key = 'customer_line_item_id'
            sequence_key = 'Insert Sequence'
            yield_key = 'Yield (ng)'
        except:
            try:
                plate_map['Plate ID']
                version = 2
                print('version',version)
                plate_key = 'Plate ID'
                well_key = 'Well Location'
                name_key = 'Name'
                sequence_key = 'Insert sequence'
                yield_key = 'Yield (ng)'
            except:
                print("Doesn't fit the standard formats, please check column names in {}".format(file))
                sys.exit(0)
        unique = plate_map[plate_key].unique()
        for plate in unique:
            # Skips all of the plates that have already been parsed and added to the database
            if plate in plates_made:
                print('Plate {} has already been added'.format(plate))
                continue
            ot.print_center("...Parsing plate {}...".format(plate))
            # Create a subset of rows pertaining to a specific plate
            current_plate_map = plate_map[plate_map[plate_key] == plate]
            current_plate = Plate('','syn_plate',plate,plate_id='syn_plate'+str(current_id).zfill(3))
            current_id += 1

            session.add(current_plate)
            for index,row in current_plate_map.iterrows():
                current_frag = ''
                current_frag = session.query(Fragment).filter(Fragment.seq == row[sequence_key]).first()
                # Checks if a corresponding fragment is found and then adds it to the plate
                if current_frag:
                    current_plate.add_item(current_frag,address=row[well_key].strip(),syn_yield=row[yield_key])
                else:
                    new_frag = Fragment(
                        fragment_name=row[name_key],
                        seq=row[sequence_key].upper(),
                        has_part='False')
                    current_plate.add_item(new_frag,address=row[well_key].strip(),syn_yield=row[yield_key])

    ot.print_center('...Updating part status...')
    for part in session.query(Part).join(Fragment,Part.fragments).join(Well,Fragment.wells)\
            .join(Plate,Well.plates).filter(Plate.plate_name.notin_(plates_made)):
        print(part.part_id)
        no_build = False
        for frag in part.fragments:
            if not frag.wells:
                no_build = True
        if not no_build:
            part.change_status('received')
        session.add(part)

    commit = int(ot.request_info("Commit changes (1-yes, 2-no): ",type='int'))
    if commit == 1:
        session.commit()
        current_ids = pd.read_sql_query("SELECT plates.plate_name,plates.plate_id FROM plates WHERE plates.plate_type = 'syn_plate'", con=engine)
        print(current_ids)
        print("Write the ID on the bottom left corner of the front edge of each plate")


if __name__ == "__main__":
    session,engine = connect_db()
    frag_assign()




#
