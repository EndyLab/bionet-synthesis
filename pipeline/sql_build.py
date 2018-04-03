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

def run_build():

    print("\n============ Beginning build ============\n")

    ## ============================================
    ## Take in required information
    ## ============================================

    # Set starting paths
    # BASE_PATH is defined in the config file
    PIPELINE_PATH = BASE_PATH + "/pipeline"
    BUILDS_PATH = BASE_PATH + "/builds"
    DATA_PATH = BASE_PATH + "/data"

    # Load files
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    args = parser.parse_args()

    ## ============================================
    ## ESTABLISH INITIAL FUNCTIONS
    ## ============================================
    def change_height(container,target):
        '''Allows for real-time calibration of the p10 single height'''
        x = 0
        print("Change height - s-g:up h-l:down x:exit")
        while x == 0:
            c = getch.getch()
            if c == "s":
                p10.robot._driver.move(z=20,mode="relative")
            elif c == "d":
                p10.robot._driver.move(z=5,mode="relative")
            elif c == "f":
                p10.robot._driver.move(z=0.5,mode="relative")
            elif c == "g":
                p10.robot._driver.move(z=0.1,mode="relative")
            elif c == "h":
                p10.robot._driver.move(z=-0.1,mode="relative")
            elif c == "j":
                p10.robot._driver.move(z=-0.5,mode="relative")
            elif c == "k":
                p10.robot._driver.move(z=-5,mode="relative")
            elif c == "l":
                p10.robot._driver.move(z=-20,mode="relative")
            elif c == "x":
                x = 1
        p10s.calibrate_position((container,target.from_center(x=0, y=0, z=-1,reference=container)))

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

    def change_plates(current_plates):
        '''Allows the user to swap plates in the middle of a protocol'''
        source_plates = {}
        plate_locations = list(zip(pd.unique(current_plates),SOURCE_SLOTS[:len(current_plates)]))
        print("plate_locations","\n", plate_locations)
        input("Switch out plates for those listed. Press enter when ready.")
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

    # Verify that the correct robot is being used
    if args.run:
        robot_name = str(os.environ["ROBOT_DEV"][-5:])
        robot_number = str(input("Run on this robot: {} ? 1-Yes, 2-No ".format(robot_name)))
        if robot_number == "1":
            print("Proceeding with run")
        else:
            sys.exit("Run roboswitch 'ROBOT_NAME' to change to the correct robot")

    target_well = well_addresses()

    ## =============================================
    ## CREATE A BUILD OF DESIRED GENES
    ## =============================================
    max_rxns = 96 # Sets the maximal number of clones you want to run
    to_build = []
    acceptable_status = ['received'] # List with all of the status that you want to pull from
    ## Pulls in parts that have the desirable status and then sorts them by the plates their fragments are in,
    ## sorts them by their plate numbers and returns the earliest set of parts
    to_build = [part for part in session.query(Part).join(Fragment,Part.fragments).\
                join(Well,Fragment.wells).join(Plate,Well.plates)\
                .filter(Part.status.in_(acceptable_status)).filter(Fragment.cloning_enzyme == 'BbsI').order_by(Plate.id)]
    to_build = to_build[:max_rxns]
    print("to_build",len(to_build))
    target_build = Build(to_build)
    session.add(target_build)

    ## =============================================
    ## CREATE PLATE GROUPS
    ## =============================================
    # Currently available spots on the OT-one deck
    SOURCE_SLOTS = ['D2','D3','B2']

    ## Generate a list of unique plates that are needed
    unique_plates = []
    for part in to_build:
        for frag in part.fragments:
            unique_plates.append(frag.wells[0].plates.plate_name)
    unique_plates = pd.Series(unique_plates).unique()

    ## Group the plates so that they can be swapped in batches
    group_plates = [unique_plates[n:n+len(SOURCE_SLOTS)] for n in range(0, len(unique_plates), len(SOURCE_SLOTS))]
    for num,group in enumerate(group_plates):
        print("Group{}: {}".format(num+1,group))

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

    # Make the dataframe to represent the OT-1 deck
    deck = ['A1','B2','C3','D2','E1']
    slots = pd.Series(deck)
    columns = sorted(slots.str[0].unique())
    rows = sorted(slots.str[1].unique(), reverse=True)
    layout_table = pd.DataFrame(index=rows, columns=columns)
    layout_table.fillna("---", inplace=True)

    # Fill in the data frame with the locations
    for obj in locations:
            layout_table.loc[locations[obj][1], locations[obj][0]] = obj

    # Displays the required plate map and waits to proceed
    print("\n Please arrange the items in the following configuration: \n")
    print(layout_table,"\n")
    input("Press enter to continue")

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

    # Declare all of the pipettes
    p10 = instruments.Pipette(
        axis='a',
        max_volume=10,
        min_volume=0.5,
        tip_racks=p10_tipracks,
        trash_container=trash,
        channels=8,
        name='p10-8',
        aspirate_speed=400,
        dispense_speed=800
    )

    p10s = instruments.Pipette(
        axis='a',
        max_volume=10,
        min_volume=0.5,
        tip_racks=p10s_tipracks,
        trash_container=trash,
        channels=1,
        name='p10-8s',
        aspirate_speed=400,
        dispense_speed=800
    )

    p200 = instruments.Pipette(
        axis='b',
        max_volume=200,
        min_volume=20,
        tip_racks=p200_tipracks,
        trash_container=trash,
        channels=1,
        name='p200-1',
        aspirate_speed=400,
        dispense_speed=800
    )

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
    source_plates = change_plates(current_group)

    for plate in target_build.plates:
        for well in plate.wells:
            for frag in well.parts.fragments:
                start_well = frag.wells[0].address
                dest_well = well.address
                gene = well.parts.part_id
                plate = frag.wells[0].plates.plate_name

                # TEMP: INCLUDE VOLUME ONCE I ACTUALLY RUN THE RESUSPENSION
                # volume = int(frag.wells[0].volume)
                if plate not in current_group:
                    plate_counter += 1
                    current_group = group_plates[plate_counter]
                    source_plates = change_plates(current_group)

                p10s.pick_up_tip()

                # Only dilutes wells that have low starting volume
                # if volume < 15:
                #     print("Diluting sample in plate {} well {} with {}uL of water".format(plate,start_well,dil_vol))
                #     p10s.transfer(dil_vol,centrifuge_tube['B1'].bottom(),source_plates[plate].wells(start_well).bottom(),new_tip='never')

                print("Transferring {} from plate {} well {} to well {}".format(gene,plate,start_well,dest_well))
                p10s.mix(3, 8, source_plates[plate].wells(start_well).bottom())

                # Checks the calibration to make sure that it can aspirate correctly
                p10s.aspirate(frag_vol,source_plates[plate].wells(start_well).bottom())
                if plate not in used_plates:
                    change_height(source_plates[plate],source_plates[plate].wells(start_well))
                p10s.dispense(frag_vol,dest_plate.wells(dest_well).bottom())
                used_plates.append(plate)
    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()
    return

if __name__ == "__main__":
    run_build()










#
