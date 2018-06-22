from opentrons import robot, containers, instruments
import argparse
import sys

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime
import getch
import shutil
import atexit

from config import *
import ot_functions as ot

import sql_resuspension

from db_config import *

def run_build(session,engine):

    ot.print_center("============ Beginning build ============")

    # Take in the 'run' argument from the command line
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    args = parser.parse_args()

    # Verify that the correct robot is being used
    if args.run:
        ot.check_robot()

    # Choose which build plan to build
    build_options = []
    for num,build in enumerate([build for build in session.query(Build).filter(Build.status == 'planning')]):
        print('{} - {}'.format(num,build.build_name))
        build_options.append(build)
    ans = ot.request_info("Enter desired plan number: ",type='int')
    target_build = build_options[ans]

    # Use that build name to create a dataframe with the information from the plan
    query = "SELECT parts.part_id,builds.build_name,part_wells.address as destination,fragments.fragment_name,frag_plates.plate_name,frag_plates.plate_id,frag_wells.address as source,frag_wells.volume FROM parts \
            INNER JOIN wells AS part_wells ON parts.id = part_wells.part_id\
            INNER JOIN plates AS part_plates ON part_wells.plate_id = part_plates.id\
            INNER JOIN builds ON part_plates.build_id = builds.id\
            INNER JOIN part_frag ON parts.id = part_frag.part_id\
            INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
            INNER JOIN wells AS frag_wells ON fragments.id = frag_wells.fragment_id\
            INNER JOIN plates AS frag_plates ON frag_wells.plate_id = frag_plates.id\
            WHERE builds.build_name = '{}'".format(target_build.build_name)

    build_plan = pd.read_sql_query(query,con=engine)
    print(build_plan)
    frags = build_plan.groupby('part_id').agg(len)
    if len(frags) == len([frag for frag in frags.fragment_name.tolist() if frag == 2]):
        print('Build only contains 2 fragment assemblies')
        num_frags = 2
        rxn_vol = 0.6
    else:
        print('Using MM for single fragment')
        rxn_vol = 0.8
        num_frags = 1

    unique_plates = build_plan.plate_id.unique().tolist()

    # Give each row a rank based on the order of the plates to sort on later
    plate_dict = dict([[y,x] for x,y in enumerate(unique_plates)])
    build_plan['plate_rank'] = build_plan.plate_id.apply(lambda x: plate_dict[x])
    build_plan = build_plan.sort_values('plate_rank')

    # Currently available spots on the OT-one deck
    SOURCE_SLOTS = ['D2','D3','B2']

    ## Generate a list of unique plates that are needed
    plate_index = [(y,x) for x,y in enumerate(unique_plates)]
    plate_index = dict(plate_index)

    ## Group the plates so that they can be swapped in batches
    ot.print_center("...Grouping plates...")
    group_plates = [unique_plates[n:n+len(SOURCE_SLOTS)] for n in range(0, len(unique_plates), len(SOURCE_SLOTS))]
    for num,group in enumerate(group_plates):
        print("Group{}: {}".format(num+1,group))

    ot.print_center("...Checking if plates need to be resuspended...")

    query_resuspend = "SELECT plates.plate_id,plates.resuspended FROM plates\
                        WHERE plates.resuspended = 'not_resuspended'\
                            AND plates.plate_id IN ({})".format(ot.list_to_string(unique_plates))

    resuspended = pd.read_sql_query(query_resuspend,con=engine)
    if len(resuspended) == 0:
        print('All plates are resuspended')
    else:
        for i,plate in resuspended.iterrows():
            ans = ot.request_info('Plate {} is not resuspended, would you like to resuspend it? y/n: '.format(plate['plate_id']),type='string')
            if ans == 'y':
                sql_resuspension.resuspension(session,engine,session.query(Plate).filter(Plate.plate_id == plate['plate_id']).one())

    input("Press enter to continue")

    ## =============================================
    ## SETUP THE OT-1 DECK
    ## =============================================

    # Specify the locations of each object on the deck
    locations = {
                "tiprack-200" : "A3",
                "tiprack-10" : "E1",
                "tiprack-10s1" : "E3",
                "tiprack-10s2" : "E2",
                "trash" : "D1",
                "PCR-strip-tall" : "C3",
                "DEST_PLATE" : "C2",
                "Tube_rack" : "B1"
            }
    # Sets the first group of plates
    used_plates = []
    plate_counter = 0
    current_group = group_plates[plate_counter]
    source_plates = ot.change_plates(locations,current_group,SOURCE_SLOTS)

    ## =============================================
    ## SETUP THE MASTER MIX
    ## =============================================

    vol = int(ot.request_info('Enter desired reaction volume (i.e. 5,10,20): ',type='int'))

    # Set the proportion of master mix to fragment to 4:1
    master_volume = rxn_vol * vol
    frag_vol = 0.2 * vol

    ot.print_center('...Calculating master mix volumes...')

    # Set a multiplier to account for pipetting error and evaporation
    extra_master = 1.3

    unique_frag = build_plan[['part_id','fragment_name','destination']].drop_duplicates()

    frag_df = unique_frag.groupby('destination').agg(len).part_id
    frag_df = frag_df.reset_index()
    frag_df = frag_df.rename(columns={'part_id':'frag_num'})
    frag_dict = dict(zip(frag_df.destination.tolist(),frag_df.frag_num.tolist()))

    build_plan['frag_num'] = build_plan.destination.apply(lambda x: frag_dict[x])

    unique_df = build_plan[['part_id','destination','frag_num']].drop_duplicates()

    total_rxns = unique_df.frag_num.sum()

    need_extra = unique_df[unique_df.frag_num > 1]

    num_wells = len(build_plan.part_id.unique().tolist())
    num_rows = num_wells // 8
    master_reactions = math.ceil((total_rxns) * extra_master)
    print("Total rxns: {}".format(total_rxns,master_reactions))

    # Generate the dataframe to present the master mix composition
    master_mix = ot.make_gg_rxns(master_reactions,master_volume)
    print("Use the table below to create the master mix")
    print()
    print(master_mix)
    print()
    print("Place the master mix in the 'A1' position of the tube rack")
    print("Also place a tube of water in the 'B1' position ")
    input("Press enter to continue")

    ## =============================================
    ## INITIALIZE THE OT-1
    ## =============================================
    # Determine whether to simulate or run the protocol

    if args.run:
        port = os.environ["ROBOT_DEV"]
        print("Connecting robot to port {}".format(port))
        robot.connect(port)
    else:
        print("Simulating protcol run")
        robot.connect()

    # Declare components on the deck
    p200_tipracks = [
        containers.load('tiprack-200ul', locations["tiprack-200"]),
    ]
    p10_tipracks = [
        containers.load('tiprack-10ul', locations["tiprack-10"]),
    ]
    p10s_tipracks = [
        containers.load('tiprack-10ul', locations["tiprack-10s1"]),
        containers.load('tiprack-10ul', locations["tiprack-10s2"])
    ]
    trash = containers.load('point', locations["trash"], 'holywastedplasticbatman')
    centrifuge_tube = containers.load('tube-rack-2ml',locations["Tube_rack"])
    master = containers.load('PCR-strip-tall', locations["PCR-strip-tall"])
    dest_plate = containers.load('96-PCR-tall', locations["DEST_PLATE"])

    p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)

    # Update database status
    def exit_handler():
        ans = ot.request_info('Restore the plan before quitting? (y/n): ',type='string')
        if ans != 'n':
            target_build.status = 'planning'
        else:
            target_build.status = 'abandoned'
            ot.print_center('...Unstaging all parts in build plan...')
            for part in session.query(Part).filter(Part.part_id.in_(build_plan.part_id.unique().tolist())):
                part.change_status('received')
        session.commit()


    target_build.status = 'building'
    session.commit()

    atexit.register(exit_handler)


    ## =============================================
    ## OT-1 PROTOCOL
    ## =============================================

    # Start timer
    start = datetime.now()
    print("Starting run at: ",start)

    # Home the robot to start
    robot.home()

   # Aliquot the master mix into the PCR tube strip
    vol_per_tube = round((num_rows * master_volume * extra_master),2)
    print("Aliquoting MM into PCR tubes")
    print("{}ul into each tube".format(vol_per_tube))
    p200.pick_up_tip()
    for well in range(8):
        print("Transferring {}ul to well {}".format(vol_per_tube,well))
        p200.transfer(vol_per_tube, centrifuge_tube['A1'].bottom(),master.wells(well).bottom(), mix_before=(3,50),new_tip='never')
    p200.drop_tip()

    # Aliquot the master mix into all of the desired wells
    p10.pick_up_tip()
    for row in range(num_rows):
        print("Transferring {}ul of master mix to row {}".format(master_volume,int(row)+1))
        p10.transfer(master_volume, master['A1'].bottom(), dest_plate.rows(row).bottom(), mix_before=(1,8), blow_out=True, new_tip='never')
    p10.drop_tip()

    # Aliquot master mix into the last row if not a complete row
    if num_wells % 8 > 0:
        p10s.pick_up_tip()
        print("need single channel for {}".format(num_wells % 8))
        for missing in range(num_wells % 8):
            current_well = (8 * num_rows) + (missing)
            print("Transferring {}ul of extra MM to {}".format(master_volume,current_well))
            p10s.transfer(master_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
        p10s.drop_tip()

    # Aliquot extra master mix into wells with multiple fragments
    if len(need_extra) != 0:
        p10s.pick_up_tip()
        for i,transfer in need_extra.iterrows():
            extra_volume = (int(transfer['frag_num']) - 1) * master_volume
            current_well = transfer['destination']
            print("Transferring {}ul of extra MM to {}".format(extra_volume,current_well))
            p10s.transfer(extra_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
        p10s.drop_tip()
    else:
        print('No extra MM must be aliquoted')

   ## Add the fragments from the source plates to the destination plate
   ## ============================================

    # Sets the volume of water to dilute with, if needed
    dil_vol = 5

    build_plan = build_plan.sort_values('plate_rank')
    for i,row in build_plan.iterrows():
        start_well = row['source']
        dest_well = row['destination']
        gene = row['part_id']
        plate = row['plate_id']
        volume = row['volume']

        if plate not in current_group:
            plate_counter += 1
            current_group = group_plates[plate_counter]
            source_plates = ot.change_plates(locations,current_group,SOURCE_SLOTS)

        p10s.pick_up_tip()

        # Only dilutes wells that have low starting volume
        if volume < 30:
            print("Diluting sample in plate {} well {} with {}uL of water".format(plate,start_well,dil_vol))
            p10s.transfer(dil_vol,centrifuge_tube['B1'].bottom(),source_plates[plate].wells(start_well).bottom(),new_tip='never')

        print("Transferring {} of {} from plate {} well {} to well {}".format(frag_vol,gene,plate,start_well,dest_well))
        p10s.mix(3, 8, source_plates[plate].wells(start_well).bottom())

        # Checks the calibration to make sure that it can aspirate correctly
        p10s.aspirate(frag_vol,source_plates[plate].wells(start_well).bottom())
        # if plate not in used_plates:
        #     ot.change_height(p10s,source_plates[plate],source_plates[plate].wells(start_well))
        p10s.dispense(frag_vol,dest_plate.wells(dest_well).bottom())
        used_plates.append(plate)

        p10s.drop_tip()

    robot.home()

    # commit = int(ot.request_info("Commit changes (1-yes, 2-no): ",type='int'))
    # if commit == 1:
    input("About to commit")
    for part in to_build:
        part.change_status('building')
    session.commit()
    return

if __name__ == "__main__":
    session,engine = connect_db()
    run_build(session,engine)










#
