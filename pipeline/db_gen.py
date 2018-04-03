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
    old_parts = [part.part_id for part in session.query(Part)]
    for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
        with open(file,"r") as json_file:
            data = json.load(json_file)
        if data['gene_id'] in old_parts:
            print('Found old part')
            continue
        new = Part(part_id=data['gene_id'],
            part_name=data['info']['documentation']['gene_name'],
            part_type=data['info']['gene_metadata']['cloning']['part_type'],
            seq=data['sequence']['optimized_sequence'].upper(),
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
                    new_frag = Fragment(
                        fragment_name=fragment,
                        seq=data['sequence']['fragment_sequences'][fragment].upper(),
                        parts=[new])
                    session.add(new_frag)
                else:
                    new_frag = Fragment(
                        fragment_name=fragment,
                        location=data['location']['fragments'][fragment],
                        seq=data['sequence']['fragment_sequences'][fragment].upper(),
                        parts=[new])
                    session.add(new_frag)
                print(new_frag.fragment_name,new_frag.parts[0].part_id)
        session.add(new)
    ## Add the extra connections in for the small linked genes
    for part,seq in zip(linked,linked_seq):
        frag = session.query(Fragment).filter(Fragment.seq == seq).first()
        part.fragments.append(frag)
        session.add(part)

    ## Convert build maps into build objects
    old_builds = [build.build_name for build in session.query(Build)]
    for file in glob.glob('{}/builds/build*/build*_20*.csv'.format(BASE_PATH)):
        # if 'build007' in file:
        #     continue
        build_name = re.match(r'.+\/builds\/(build[0-9]{3})\/.+',file).groups()[0]
        if build_name in old_builds:
            print('Found old build')
            continue
        print("Adding build: ",build_name)
        plate_map = pd.read_csv(file)
        parts = plate_map['Gene'].tolist()
        part_objs = [session.query(Part).filter(Part.part_id == part).one() for part in parts if part != 'Empty' and 'Unk' not in part]
        current_build = Build(part_objs,build_name=build_name)
        current_build.status = 'building'
        session.add(current_build)
        for plate in current_build.plates:
            for well in plate.wells:
                well.parts.change_status('building')
                session.add(well)

    for target_build in session.query(Build):
        if target_build in old_builds:
            print('Found old target build')
            continue
        seq_plate = Plate([],'seq_plate',target_build.build_name)
        target_build.plates.append(seq_plate)
        target_build.status = 'sequencing'
        session.add(seq_plate)
        for well in target_build.plates[0].wells:
            seq_plate.add_item(well.parts,well.address)
            well.parts.change_status('sequencing')
            session.add(seq_plate)
    session.commit()
    print('\nCreated database\n')
    return


if __name__ == "__main__":
    build_db()




    #
