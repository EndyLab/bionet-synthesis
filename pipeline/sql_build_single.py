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
session,engine = connect_db()

def run_build():

    ## ============================================
    ## ESTABLISH INITIAL FUNCTIONS
    ## ============================================

    def change_plates(locations,current_plates):
        '''Allows the user to swap plates in the middle of a protocol'''
        source_plates = {}
        temp = locations
        plate_locations = list(zip(pd.unique(current_plates),SOURCE_SLOTS[:len(current_plates)]))
        temp.update(plate_locations)
        ot.print_layout(temp)
        for plate, slot in plate_locations:
            source_plates[plate] = containers.load('96-flat', slot)
        return source_plates

    def make_gg_rxns(num_rxns,rxn_vol):
        '''Calculates the amount of each reagent to add to reach the desired master mix'''
        cutsmart = 1 * num_rxns
        atp = 1 * num_rxns
        vector = 0.25 * num_rxns
        ligase = 0.5 * num_rxns
        enzyme = 0.25 * num_rxns
        water = (rxn_vol - ((cutsmart + atp + vector + ligase + enzyme)/num_rxns)) * num_rxns
        master_mix = pd.DataFrame(
            {'Component':['H2O','Cutsmart','ATP','Vector','T4 Ligase','Restriction Enzyme','Total'],
            'Amount':[water,cutsmart,atp,vector,ligase,enzyme,rxn_vol*num_rxns]},
            columns=["Component","Amount"]
        )
        return master_mix

    ot.print_center("============ Beginning build ============")

    # Take in the 'run' argument from the command line
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    args = parser.parse_args()

    # Verify that the correct robot is being used
    if args.run:
        robot_name = str(os.environ["ROBOT_DEV"][-5:])
        robot_number = str(input("Run on this robot: {} ? 1-Yes, 2-No ".format(robot_name)))
        if robot_number == "1":
            print("Proceeding with run")
        else:
            sys.exit("Run roboswitch 'ROBOT_NAME' to change to the correct robot")

    ## =============================================
    ## CREATE A BUILD OF DESIRED GENES
    ## =============================================
    # max_rxns = 96 # Sets the maximal number of clones you want to run
    max_rxns = int(input("Enter the desired number of parts: "))
    enzyme = input("Which enzyme: (1 - BbsI or 2 - BtgZI)")
    if enzyme == '':
        enzyme = 'BbsI'
    elif int(enzyme) == 1:
        enzyme == 'BbsI'
    else:
        enzyme = 'BtgZI'
    print("Using: ",enzyme)
    ot.print_center("...Finding genes to build...")

    to_build = []
    acceptable_status = ['received'] # List with all of the status that you want to pull from
    rework = ['cloning_error','cloning_failure','trans_failure']

    ## Pulls in parts that have the desirable status and then sorts them by the plates their fragments are in,
    ## sorts them by their plate numbers and returns the earliest set of parts
    to_build += [part for part in session.query(Part).join(Fragment,Part.fragments).\
                join(Well,Fragment.wells).join(Plate,Well.plates)\
                .filter(Part.status.in_(acceptable_status)).filter(Part.cloning_enzyme == enzyme).order_by(Plate.id)]
    if len(to_build) < 96:
        to_build += [part for part in session.query(Part).join(Fragment,Part.fragments).\
                    join(Well,Fragment.wells).join(Plate,Well.plates)\
                    .filter(Part.status.in_(rework)).filter(Part.cloning_enzyme == enzyme).order_by(Plate.id)]

    to_build = to_build[:max_rxns]
    print("to_build",len(to_build))
    part_names = [[part.part_id,part.status] for part in to_build]

    prev_builds = [build.build_name for build in session.query(Build)]
    build_numbers = [int(name[-3:]) for name in prev_builds if name != '']
    last_build = max(build_numbers)

    build_num = str(last_build+1).zfill(3)
    build_name = 'build' + build_num

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

    path = '{}/builds/build{}'.format(BASE_PATH,build_num)
    if os.path.exists(path):
        print("Directory build{} already exists".format(build_num))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)
        print("Making directory for build{}".format(build_num))
    export_map.to_csv('{}/build{}_trans_map.csv'.format(path,build_num),index=False)
    input("check map")

    ## =============================================
    ## CREATE PLATE GROUPS
    ## =============================================
    # Currently available spots on the OT-one deck
    SOURCE_SLOTS = ['D2','D3','B2']

    ## Generate a list of unique plates that are needed
    ot.print_center("...Finding the required plates...")
    unique_plates = []
    for part in to_build:
        for frag in part.fragments:
            unique_plates.append(frag.wells[0].plates.plate_id)
    unique_plates = pd.Series(unique_plates).unique()
    plate_index = [(y,x) for x,y in enumerate(unique_plates)]
    plate_index = dict(plate_index)

    ## Group the plates so that they can be swapped in batches
    ot.print_center("...Grouping plates...")
    group_plates = [unique_plates[n:n+len(SOURCE_SLOTS)] for n in range(0, len(unique_plates), len(SOURCE_SLOTS))]
    for num,group in enumerate(group_plates):
        print("Group{}: {}".format(num+1,group))

    print("Checking if plates need to be resuspended")
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
    ot.print_layout(locations)

    ## =============================================
    ## SETUP THE MASTER MIX
    ## =============================================

    # Set the volume of master mix to use per reaction
    master_volume = 8
    # Set a multiplier to account for pipetting error and evaporation
    extra_master = 1.3

    need_extra = []
    for plate in target_build.plates:
        need_extra.append([{'plate':well.plates.plate_name,'well':well.address,'frags':len(well.parts.fragments)}\
                          for well in plate.wells if len(well.parts.fragments) > 1])
    total_extra = 0
    for plates in need_extra:
        for transfer in plates:
            total_extra += transfer['frags']
    num_reactions = len(to_build)
    num_rows = num_reactions // 8
    master_reactions = math.ceil((total_extra + num_reactions) * extra_master)
    print("Extra needed: {}, Total rxns: {}".format(total_extra,master_reactions))

    # Generate the dataframe to present the master mix composition
    master_mix = make_gg_rxns(master_reactions,master_volume)
    print("Use the table below to create the master mix")
    print()
    print(master_mix)
    print()
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
    vol_per_tube = math.ceil(num_rows * master_volume * extra_master)
    print("Aliquoting MM into PCR tubes")
    print("{}ul into each tube".format(vol_per_tube))
    p200.pick_up_tip()
    for well in range(8):
        print("Transferring {}ul to well {}".format(vol_per_tube,well))
        p200.transfer(vol_per_tube, centrifuge_tube['A1'].bottom(),master.wells(well).bottom(), mix_before=(3,50),new_tip='never')
    p200.drop_tip()

    # REVIEW:
    input("Run other build")

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
    # TEMP: Only allows for a single destination plate
    p10s.pick_up_tip()
    for transfer in need_extra[0]:
        extra_volume = (int(transfer['frags']) - 1) * master_volume
        current_well = transfer['well']
        print("Transferring {}ul of extra MM to {}".format(extra_volume,current_well))
        p10s.transfer(extra_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
    p10s.drop_tip()

    ## Add the fragments from the source plates to the destination plate
    ## ============================================
    # Sets the volume to pipet of each fragment to 2uL
    frag_vol = 2

    # Sets the volume of water to dilute with, if needed
    dil_vol = 5

    # Sets the first group of plates
    used_plates = []
    plate_counter = 0
    current_group = group_plates[plate_counter]
    source_plates = change_plates(locations,current_group)

    indexes = []
    rows = []
    for plate in target_build.plates:
        for well in plate.wells:
            for frag in well.parts.fragments:
                indexes.append(plate_index[frag.wells[0].plates.plate_id])
                rows += [[frag.wells[0].address,well.address,well.parts.part_id,frag.wells[0].plates.plate_id,frag.wells[0].volume]]

    build_plan = [row for i,row in sorted(zip(indexes,rows))]
    for row in build_plan:
        start_well = row[0]
        dest_well = row[1]
        gene = row[2]
        plate = row[3]
        volume = row[4]

        if plate not in current_group:
            plate_counter += 1
            current_group = group_plates[plate_counter]
            source_plates = change_plates(locations,current_group)

        p10s.pick_up_tip()

        # Only dilutes wells that have low starting volume
        if volume < 30:
            print("Diluting sample in plate {} well {} with {}uL of water".format(plate,start_well,dil_vol))
            p10s.transfer(dil_vol,centrifuge_tube['B1'].bottom(),source_plates[plate].wells(start_well).bottom(),new_tip='never')

        print("Transferring {} from plate {} well {} to well {}".format(gene,plate,start_well,dest_well))
        p10s.mix(3, 8, source_plates[plate].wells(start_well).bottom())

        # Checks the calibration to make sure that it can aspirate correctly
        p10s.aspirate(frag_vol,source_plates[plate].wells(start_well).bottom())
        if plate not in used_plates:
            ot.change_height(p10s,source_plates[plate],source_plates[plate].wells(start_well))
        p10s.dispense(frag_vol,dest_plate.wells(dest_well).bottom())
        used_plates.append(plate)

        p10s.drop_tip()

    input("stop again")

    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()
    return

if __name__ == "__main__":
    run_build()










#
