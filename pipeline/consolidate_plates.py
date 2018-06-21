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

import ot_functions as ot
from config import *
from db_config import *
session,engine = connect_db()

# Create a csv file that has a build_name and failed column
# The failed column should be a string of well addresses with a space
# separating them
build_failures = pd.read_csv('./testing/plate_consolidation.csv')

# Separate the well addresses into a list
build_failures['wells'] = build_failures.failed.apply(lambda x: str(x).split(' '))

# State the build to be consolidated into the others
builds = build_failures.build_name.tolist()
for i, build in enumerate(builds):
    print('{}: {}'.format(i,build))
ans = ot.request_info('Choose which build to consolidate into the others: ',type='int')
extra_build = builds[ans]
    
# Pull out the build to consolidate
consolidate = build_failures[build_failures.build_name == extra_build]
consol_obj = session.query(Build).filter(Build.build_name == extra_build).one()

# Add in the result from the transformation on the consolidation plate
extras = []
for well in consol_obj.plates[0].wells:
    if well.address in consolidate.iloc[0].wells:
        print('fail')
        well.trans_outcome = 'trans_failure'
    else:
        well.trans_outcome = 'trans_success'
        extras.append(well.parts) # Store the parts that can be transferred

extras = extras[::-1] # Invert the string to the pop items from the beginning

other = build_failures[build_failures.build_name != extra_build]

for i,row in other.iterrows():

    # Find the corresponding build and link it to a seq plate
    target_build = session.query(Build).filter(Build.build_name == row['build_name']).one()
    print(target_build.build_name)
    seq_plate = Plate([],'seq_plate','{}-seq'.format(target_build.build_name))
    target_build.plates.append(seq_plate)
    target_build.status = 'sequencing'
    session.add(seq_plate)

    # Update the status of each well in the assembly plate
    for well in session.query(Well).join(Plate,Well.plates).join(Build,Plate.builds).filter(Build.build_name == row['build_name']):
        if well.address in row['wells']:
            well.trans_outcome = 'trans_failure'
            if len(extras) > 0:
                # If the well didn't work, replace it with one of the available extra wells
                new_part = extras.pop()
                print(new_part.part_id)
                seq_plate.add_item(new_part,address=well.address)
        else:
            well.trans_outcome = 'trans_success'
            new_part = well.parts
            seq_plate.add_item(new_part,address=well.address)

session.commit()


builds = build_failures.build_name.tolist()

for build in builds:

    # Creates a df by going from build to seq plate to well to parts to find what parts
    # are in each plate
    query_build = "SELECT parts.part_id,builds.build_name,wells.address from builds\
                    INNER JOIN plates ON plates.build_id = builds.id\
                    INNER JOIN wells ON plates.id = wells.plate_id\
                    INNER JOIN parts ON wells.part_id = parts.id\
                    WHERE builds.build_name = '{}'\
                        AND plates.plate_type = 'seq_plate'".format(build)

    df = pd.read_sql_query(query_build,con=engine)

    # Create a list of parts to the query in the next step
    parts = df.part_id.tolist()
    df = df.set_index('part_id')

    # Use the parts in the last query to find the original builds that they came from
    query_original = "SELECT parts.part_id,builds.build_name AS original_build,wells.address AS original_well from parts\
                    INNER JOIN wells ON parts.id = wells.part_id\
                    INNER JOIN plates ON plates.id = wells.plate_id\
                    INNER JOIN builds ON plates.build_id = builds.id\
                    WHERE wells.plate_type = 'assembly_plate'\
                        AND parts.part_id IN ({})\
                        AND builds.build_name IN ({})".format(ot.list_to_string(parts),ot.list_to_string(builds))

    df_ori = pd.read_sql_query(query_original,con=engine)
    df_ori = df_ori.set_index('part_id')

    # Merge the two dataframes to get both the original well and the destination in one df
    full = pd.concat([df,df_ori],axis=1)

    # Sort the wells based to go A1-H12
    well_addresses = ot.well_addresses()
    well_rank = dict(zip(well_addresses,range(len(well_addresses))))
    full['well_rank'] = full.address.apply(lambda x: well_rank[x])
    full = full.sort_values('well_rank')

    # Generate the seq maps
    full = full[['original_build','original_well','build_name','address']]
    path = '{}/builds/{}'.format(BASE_PATH,build)
    ot.make_directory(path)
    full.to_csv('{}/{}_seq_plate_map.csv'.format(path,build))





    #
