import argparse
import sys

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import re
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime

from config import *
from db_config import *

to_build = []
priority = ['pSHPs0325B569005MU','pSHPs0325B569008MU','pSHPs0325B569010MU']
acceptable_status = ['received']

to_build = [part for part in session.query(Part).join(Fragment,Part.fragments).\
            join(Well,Fragment.wells).join(Plate,Well.plates).filter(Plate.plate_name.in_(priority))\
            .filter(Part.status.in_(acceptable_status)).filter(Fragment.cloning_enzyme == 'BbsI').order_by(Plate.id)]
# to_build = [part for part in session.query(Part).join(Fragment,Part.fragments).\
#             join(Well,Fragment.wells).join(Plate,Well.plates).filter(Plate.plate_name.in_(priority))\
#             .filter(Part.status.in_(acceptable_status)).order_by(Plate.id)]
for part in to_build:
    print(part.part_id,part.status)
    for frag in part.fragments:
        print(frag.fragment_name,frag.cloning_enzyme)
# print(to_build)
print(len(to_build))



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

for instance in session.query(Plate).filter(Plate.plate_name == 'pSHPs0325B569005MU'):
    print(instance.plate_name)
    fragment_names = [well.fragments.fragment_name for well in instance.wells]
    used = [well.address for well in instance.wells]
    # if '38' in instance.plate_name:
    print(len(instance.wells))
    # for well in instance.wells:
    #     print(well.address,well.syn_yield,well.volume)
    print(instance.plate_name)
    # input("plate done")
#
# plate_map = pd.read_csv('../plate_maps/O-008.1_A-001_2018.03.26.csv')
# original = plate_map['Name'].tolist()
# ori_wells = plate_map['Well Location'].tolist()
#
# missing = [name for name in original if name not in fragment_names]
# print(missing)
# missing_wells = [well for well in ori_wells if well not in used]
# print(missing_wells)
#
# plate = session.query(Plate).filter(Plate.plate_name == 'pSHPs0325B569008MU').one()
# print(plate.plate_name,plate.wells[0].fragments.fragment_name)

# part = session.query(Part).filter(Part.part_name == 'WD_RS05740').one()
# print(part.part_id,part.part_name)

# names = []
# id_nums = []
#
# for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
#     with open(file,"r") as json_file:
#         data = json.load(json_file)
#     names.append(data['info']['documentation']['gene_name'])
#     id_nums.append(data['gene_id'])
#     if data['info']['documentation']['gene_name'] == 'WRI_RS02680':
#         print(data['gene_id'])
#         input()
# dictionary = dict(zip(names,id_nums))
# print(dictionary)
# print(names)
