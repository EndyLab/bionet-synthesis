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
    def well_addresses():
        '''Generates a list of well address A1-H12'''
        letter = ["A","B","C","D","E","F","G","H"]
        number = ["1","2","3","4","5","6","7","8","9","10","11","12"]
        target_well = []
        temp_well = 0
        for n in number:
            for l in letter:
                temp_well = l + n
                target_well.append(temp_well)
        return target_well

    wells = well_addresses()

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
    plates_made = [plate.plate_name for plate in session.query(Plate).filter(Plate.plate_type == 'syn_plate')]
    print("Plates that have already been parsed:\n",plates_made)

    ## Take in plate maps from twist and generate fragment objects
    not_found = []
    making = []
    current_plate = Plate([],'','')
    for file in glob.glob('{}/plate_maps/*.csv'.format(BASE_PATH)):
        plate_map = pd.read_csv(file)
        for index,row in plate_map.iterrows():
            ## Have to account for two different csv formats
            try:
                row['Plate']
                plate_id = row['Plate'].strip()
                fragment_name = row['customer_line_item_id'].strip()
                insert = row['Insert Sequence']
                well_address = row['Well']
                format = 1
            except:
                try:
                    row['Plate ID']
                    plate_id = row['Plate ID'].strip()
                    fragment_name = row['Name'].strip()
                    insert = row['Insert sequence']
                    well_address = row['Well Location']
                    format = 2
                except:
                    print("Doesn't fit either acceptable formats")
                    continue

            if plate_id in plates_made:
                continue
            elif plate_id not in making:
                # if format == 1:
                #     subset = plate_map[plate_map["Plate"] == plate_id]
                #     ori_wells = subset['Well'].tolist()
                # elif format == 2:
                #     subset = plate_map[plate_map["Plate ID"] == plate_id]
                #     ori_wells = subset['Well Locations'].tolist()
                # print("original\n",ori_wells)
                # used_wells = [well.address for well in current_plate.wells]
                # used_wells = [well.address for well in session.query(Well).join(Plate,Well.plates).filter(Plate.plate_name == plate_id)]
                # print("Used\n",used_wells)
                # missing = [well for well in ori_wells if well not in used_wells]
                # print("Missing\n",missing)
                # input('inspect missing')
                current_plate = Plate('','syn_plate',plate_name=plate_id)
                session.add(current_plate)
                print('Adding plate: ',plate_id)
                making.append(plate_id)

            if fragment_name[:-2] not in names and fragment_name[:-2] not in id_nums:
                not_found.append(fragment_name)
                continue
            current_frag = ''
            current_frag = session.query(Fragment).filter(Fragment.seq == insert).first()
            if current_frag:
                current_plate.add_item(current_frag,well_address,syn_yield=row['Yield (ng)'])

    ## Determie which samples are build ready
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
