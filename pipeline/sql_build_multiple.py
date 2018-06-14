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

    build_options = []
    for num,file in enumerate(sorted(glob.glob('{}/builds/*/*_plan.csv'.format(BASE_PATH)))):
        print('{} - {}'.format(num,file.split("/")[-1]))
        build_options.append(file)
    ans = ot.request_info("Enter desired plan number: ",type='int')
    target_file = build_options[int(ans)]
    build_plan = pd.read_csv(target_file)
    build_parts = build_plan.Gene.tolist()

    to_build = [part for part in session.query(Part).filter(Part.part_id.in_(build_parts))]

    parts_df = ot.query_for_fragments(build_parts,engine)
    frags = parts_df.groupby('part_id').agg(len)
    if len(frags) == len([frag for frag in frags.fragment_name.tolist() if frag == 2]):
        print('Build only contains 2 fragment assemblies')
        num_frags = 2
        rxn_vol = 0.6
    else:
        rxn_vol = 0.8
        num_frags = 1

    build_name = target_file.split('/')[-1].split('_')[0]
    # print(build_name)

    target_build = Build(to_build,build_name=build_name)
    session.add(target_build)

    building = []
    addresses = []
    for well in target_build.plates[0].wells:
        building.append(well.parts.part_id)
        addresses.append(well.address)
    export_map = pd.DataFrame({
        'Gene' : building,
        'Destination' : addresses
    })
    export_map = export_map[['Gene','Destination']]

    build_path = '{}/builds/{}'.format(BASE_PATH,build_name)

    if args.run:
        ot.make_directory(build_path)
        export_map.to_csv('{}/{}_trans_map.csv'.format(build_path,build_name),index=False)
        os.remove(glob.glob('{}/{}_plan.csv'.format(build_path,build_name))[0])

    ## =============================================
    ## CREATE PLATE GROUPS
    ## =============================================

    ot.print_center("...Finding the required plates...")

    parts_df = ot.query_for_plates(build_parts,engine)
    unique_plates = parts_df.plate_id.unique().tolist()

    ## Add required columns to the dataframe

    # Rename the columnn names to reflect the values
    parts_df = parts_df.rename(columns={'id':'destination','plate_name':'plate_rank'})

    # Add the destination wells to each row which represents a single fragment
    dest_dict = dict(zip(building,addresses))
    parts_df.destination = parts_df.part_id.apply(lambda x: dest_dict[x])

    # Give each row a rank based on the order of the plates to sort on later
    plate_dict = dict([[y,x] for x,y in enumerate(unique_plates)])
    parts_df.plate_rank = parts_df.plate_id.apply(lambda x: plate_dict[x])

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
    for plate in unique_plates:
        current_plate = session.query(Plate).filter(Plate.plate_id == plate).one()
        if current_plate.resuspended == 'not_resuspended':
            ans = input("Plate {} is not resuspended, would you like to resuspend it? y/n: ".format(plate))
            if ans == 'y':
                sql_resuspension.resuspension(current_plate)
        else:
            print(plate,current_plate.resuspended)
    input("Press enter to continue")

    # Create a means to index to the plates
    plate_dict = dict(enumerate(unique_plates))

    ## =============================================
    ## SETUP THE OT-1 DECK
    ## =============================================
    SOURCE_SLOTS = ['D2','D3','B2']

    # Specify the locations of each object on the deck
    locations = {
                "tiprack-200" : "A3",
                "tiprack-10" : "E2",
                "tiprack-10s1" : "E3",
                "tiprack-10s2" : "E1",
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

    need_extra = []
    for plate in target_build.plates:
        need_extra.append([{'plate':well.plates.plate_name,'well':well.address,\
        'frags':len(well.parts.fragments)} for well in plate.wells\
        if len(well.parts.fragments) > num_frags])

    total_extra = 0
    for plates in need_extra:
        for transfer in plates:
            total_extra += transfer['frags']
    num_reactions = len(to_build)
    num_rows = num_reactions // 8
    master_reactions = math.ceil((total_extra + num_reactions) * extra_master)
    print("Extra needed: {}, Total rxns: {}".format(total_extra,master_reactions))

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
    if num_reactions % 8 > 0:
        p10s.pick_up_tip()
        print("need single channel for {}".format(num_reactions % 8))
        for missing in range(num_reactions % 8):
            current_well = (8 * num_rows) + (missing)
            print("Transferring {}ul of extra MM to {}".format(master_volume,current_well))
            p10s.transfer(master_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
        p10s.drop_tip()

    # Aliquot extra master mix into wells with multiple fragments
    if len(need_extra[0]) != 0:
        p10s.pick_up_tip()
        for transfer in need_extra[0]:
            extra_volume = (int(transfer['frags']) - 1) * master_volume
            current_well = transfer['well']
            print("Transferring {}ul of extra MM to {}".format(extra_volume,current_well))
            p10s.transfer(extra_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
        p10s.drop_tip()
    else:
        print('No extra MM must be aliquoted')

    ## Add the fragments from the source plates to the destination plate
    ## ============================================

    # Sets the volume of water to dilute with, if needed
    dil_vol = 5

    parts_df = parts_df.sort_values('plate_rank')
    for i,row in parts_df.iterrows():
        start_well = row['address']
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

    commit = int(ot.request_info("Commit changes (1-yes, 2-no): ",type='int'))
    if commit == 1:
        for part in to_build:
            part.status = 'building'
        session.commit()
    return

if __name__ == "__main__":
    session,engine = connect_db()
    run_build(session,engine)










#
