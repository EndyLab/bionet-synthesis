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

def assess_plate():
    # Determine which build/plate to assess
    assemblies = []
    print("Choose which build you would like to assess:")
    for index,assembly in enumerate(session.query(Build).filter(Build.status == 'building').order_by(Build.build_name)):
        print("{}. {}".format(index,assembly.build_name))
        assemblies.append(assembly)
    build_num = int(input("Enter build here: "))
    target_build = assemblies[build_num]
    print("Will assess: ",target_build.build_name)

    ## =============================================
    ## QUERY DATABASE FOR REPLACEMENTS
    ## =============================================

    # Take in the wells that didn't work
    print("Enter well addresses that didn't work. e.g. A1 H6")
    failed = [x for x in input("List: ").split(" ")]
    print("Will exclude the following:",failed)

    ot.print_center('...Generating sequencing plate...')

    # Generate an empty sequencing plate and update the build status
    seq_plate = Plate([],'seq_plate','{}-seq'.format(target_build.build_name))
    target_build.plates.append(seq_plate)
    target_build.status = 'sequencing'
    session.add(seq_plate)

    # Change the status of the wells and add successful ones to seq_plate
    for well in target_build.plates[0].wells:
        if well.address in failed:
            well.trans_outcome = 'trans_failure'
        else:
            well.trans_outcome = 'trans_success'
            seq_plate.add_item(well.parts,address=well.address)

    # Query the plate to find how many wells it has already used
    used_wells = [well for well in seq_plate.wells]
    remaining = 96 - len(used_wells)
    print("Remaining: ",remaining)
    if remaining > 0:
        ans = input("backfill empty wells? (y/n): ")
        if ans == 'n':
            remaining = 0
        # Pull parts that need to the retried
        acceptable_status = ['sequence_failure','cloning_mutation'] # List with all of the status that you want to pull from
        replacement_parts = [part for part in session.query(Part).join(Well,Part.wells)\
                    .join(Plate,Well.plates).join(Build,Plate.builds).filter(Part.status\
                    .in_(acceptable_status)).order_by(Build.id.desc(),Well.id)]

        # Find the wells that house the parts with the desired status
        rep_wells = []
        for part in replacement_parts:
            skip = False
            for well in part.wells:
                if skip:
                    continue
                if well.seq_outcome == part.status:
                    rep_wells.append(well)
                    skip = True
                    continue

        # Limit the number of wells to the number required and add it to the seq plate
        rep_wells = rep_wells[:remaining]
        part_to_well = dict(zip(replacement_parts[:remaining],rep_wells))
        for well in rep_wells:
            seq_plate.add_item(well.parts)
            print(well.parts.part_id,well.plates.builds.build_name,well.address)

        # Generate a dictionary to link the parts to the well with the desired status
        part_to_well = dict(zip(replacement_parts[:remaining],rep_wells))
    else:
        part_to_well = {}

    # Generate a dictionary to sort on well addresses A1-H1 instead of A1-A12
    well_index = ot.well_addresses()
    well_index = enumerate(well_index)
    well_index = [(y,x) for x,y in well_index]
    well_index = dict(well_index)

    # Iterate through the plate and add the necessary information for the map
    parts = []
    build_names = []
    source_wells = []
    target_wells = []
    target_index = []
    # Add all of the needed info to build a map
    for well in seq_plate.wells:
        parts.append(well.parts.part_id)

        # Finds parts coming from a different build and specifies where its coming from
        if well.parts in part_to_well.keys():
            build_names.append('__'+part_to_well[well.parts].plates.builds.build_name+'__')
            source_wells.append(part_to_well[well.parts].address)
        else:
            build_names.append(well.plates.builds.build_name)
            source_wells.append(well.address)
        target_wells.append(well.address)
        target_index.append(well_index[well.address])

    # Build the plate map and sort based on the well index
    seq_plate_map = pd.DataFrame({
        'Part' : parts,
        'Build' : build_names,
        'Source Well' : source_wells,
        'Target Well' : target_wells,
        'Target Index' : target_index
    })
    seq_plate_map = seq_plate_map.set_index('Target Index')
    seq_plate_map = seq_plate_map.sort_index()
    print(seq_plate_map)

    # Print the plates that are required
    unique_builds = pd.Series(seq_plate_map['Build']).unique()
    print("You will need the following builds:")
    print(unique_builds)

    # Generate a sequence plate map
    path = '{}/builds/{}/{}_seq_plate.csv'.format(BASE_PATH,target_build.build_name,target_build.build_name)
    seq_plate_map.to_csv(path,index=False)

    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()
    return

if __name__ == "__main__":
    assess_plate()








#
