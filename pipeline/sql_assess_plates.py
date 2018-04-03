from opentrons import robot, containers, instruments
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
import getch

from config import *
from db_config import *
from db_gen import *
from sql_frag_loc import *
from sql_resuspension import *
from sql_build import *
from sql_plating import *

session = build_db()
# session = frag_assign(session)
# session = resuspension(session)
# session = resuspension(session)
# session = run_build(session)
# session = plate(session)

assemblies = []
print("Choose which Build you would like to assess:")
for index,assembly in enumerate(session.query(Build).filter(Build.status == 'building')):
    print("{}. {}-{}".format(index,assembly.id,assembly.build_name))
    assemblies.append(assembly)

build_num = int(input("Enter build here: "))
target_build = assemblies[build_num]

## =============================================
## QUERY DATABASE FOR REPLACEMENTS
## =============================================

print("Enter well addresses that didn't work. e.g. A1 H6")
failed = [x for x in input("List: ").split(" ")]
print(failed)

# TODO: Need to set the criteria that I am going to use to pull the replacements
replacements = [part for part in session.query(Part).limit(len(failed))]
for part in replacements:
    print(part.part_id)


seq_plate = Plate([],'seq_plate',target_build.build_name)
target_build.plates.append(seq_plate)
target_build.status = 'sequencing'
session.add(seq_plate)
for well in target_build.plates[0].wells:
    print(well.parts.part_id,well.address)
    if well.address in failed:
        print("found failure")
    seq_plate.add_item(well.parts,well.address)
    well.parts.change_status('sequencing')
    session.add(seq_plate)

for well in seq_plate.wells:
    print(well.parts.part_id,well.address)

for build in session.query(Build):
    for plate in build.plates:
        print(build.build_name,plate.plate_type)










#
