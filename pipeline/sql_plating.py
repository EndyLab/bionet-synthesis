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
session,engine = connect_db()

def plate():
    print("\n============ Beginning plating ============\n")

    ## Establish initial functions
    ## ============================================

    def change_height(container,target):
        x = 0
        counter = 0
        print("Change height - s-g:up h-l:down x:exit")
        while x == 0:
            c = getch.getch()
            if c == "s":
                print("Up 20mm")
                p10.robot._driver.move(z=20,mode="relative")
            elif c == "d":
                print("Up 5mm")
                p10.robot._driver.move(z=5,mode="relative")
            elif c == "f":
                print("Up 0.5mm")
                p10.robot._driver.move(z=0.5,mode="relative")
            elif c == "g":
                print("Up 0.1mm")
                p10.robot._driver.move(z=0.1,mode="relative")
            elif c == "h":
                print("Down 0.1mm")
                p10.robot._driver.move(z=-0.1,mode="relative")
            elif c == "j":
                print("Down 0.5mm")
                p10.robot._driver.move(z=-0.5,mode="relative")
            elif c == "k":
                print("Down 5mm")
                p10.robot._driver.move(z=-5,mode="relative")
            elif c == "l":
                print("Down 20mm")
                p10.robot._driver.move(z=-20,mode="relative")
            elif c == "x":
                print("Exit")
                x = 1
            counter += 1
        if counter > 1:
            print("Will recalibrate")
            redo = True
        else:
            print("Calibrated")
            redo = False
        p10.calibrate_position((container,target.from_center(x=0, y=0, z=-1,reference=container)))

        return redo

    ## Take in required information
    ## ============================================

    PIPELINE_PATH = BASE_PATH + "/pipeline"
    BUILDS_PATH = BASE_PATH + "/builds"
    DATA_PATH = BASE_PATH + "/data"

    BACKBONE_PATH = BASE_PATH + "/sequencing_files/popen_v1-1_backbone.fasta"
    DICTIONARY_PATH = PIPELINE_PATH + "/testing/data_testing/10K_CDS.csv"

    # Take in command line arguments
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    parser.add_argument('-m', '--manual', required=False, action="store_true", help="Maunal entry of parameters.")
    args = parser.parse_args()

    # Verify that the correct robot is being used
    if args.run:
        robot_name = str(os.environ["ROBOT_DEV"][-5:])
        robot_number = str(input("Run on this robot: {} ? 1-Yes, 2-No ".format(robot_name)))
        if robot_number == "1":
            print("Proceeding with run")
        else:
            sys.exit("Run . robot.sh while in the /opentrons/robots directory to change the robot")

    assemblies = []
    print("Choose which plate you would like to transform/plate:")
    for index,assembly in enumerate(session.query(Plate).join(Build,Plate.builds).filter(Plate.plated == 'not_plated').order_by(Build.build_name)):
        print("{}. {}".format(index,assembly.builds.build_name))
        assemblies.append(assembly)

    plate_num = int(input("Enter plate here: "))
    target_plate = assemblies[plate_num]
    build_map = target_plate.wells
    if len(target_plate.wells) > 48:
        print("Too many samples to plate at once")
        portion = int(input("Choose which half to plate, 1 or 2: "))
        if portion == 1:
            build_map = target_plate.wells[:48]
            print(build_map)
        else:
            build_map = target_plate.wells[48:]
            print(build_map)
        num_reactions = len(build_map)
    else:
        portion = 1
        num_reactions = len(build_map)

    num_rows = math.ceil(num_reactions / 8)
    print("num rows: ",num_rows)
    trans_per_plate = 3
    num_plates = num_rows // trans_per_plate

    agar_plate_names = []
    for plate_num in range(num_plates):
        if portion == 2:
            plate_num += 2
        current_plate = target_plate.builds.build_name + "_p" + str(plate_num + 1)
        agar_plate_names.append(current_plate)
    print(agar_plate_names)
    print("You will need {} agar plates".format(len(agar_plate_names)))

    ## =============================================
    ## SETUP THE OT-1 DECK
    ## =============================================
    # Allocate slots for the required agar plates
    AGAR_SLOTS = ['D2','D3']
    layout = list(zip(agar_plate_names,AGAR_SLOTS[:len(agar_plate_names)]))

    # Specify the locations of each object on the deck
    locations = {
                "tiprack-200" : "A3",
                "tiprack-10_2" : "E3",
                "tiprack-10_3" : "E2",
                "tiprack-10_1" : "E1",
                "trash" : "D1",
                "PCR-strip-tall" : "C3",
                "Transformation" : "C2",
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

    ## Initialize the OT-1
    ## ============================================

    # Determine whether to simulate or run the protocol
    if args.run:
        port = os.environ["ROBOT_DEV"]
        print("Connecting robot to port {}".format(port))
        robot.connect(port)
    else:
        print("Simulating protcol run")
        robot.connect()

    # Start timer
    start = datetime.now()
    print("Starting run at: ",start)

    # Start up and declare components on the deck
    robot.home()

    p200_tipracks = [
        containers.load('tiprack-200ul', locations['tiprack-200']),
    ]

    p10_tipracks = [
        containers.load('tiprack-10ul', locations['tiprack-10_2']),
        containers.load('tiprack-10ul', locations['tiprack-10_1']),
        containers.load('tiprack-10ul', locations['tiprack-10_3'])
    ]

    transformation_plate = containers.load('96-PCR-tall', locations['Transformation'])
    trash = containers.load('point', locations['trash'], 'holywastedplasticbatman')
    centrifuge_tube = containers.load('tube-rack-2ml',locations['Tube_rack'])
    master = containers.load('PCR-strip-tall', locations['PCR-strip-tall'])

    agar_plates = {}
    for plate, slot in layout:
        agar_plates[plate] = containers.load('96-deep-well', slot)
        print("agar_plates", agar_plates[plate])
    print("agar_plates",agar_plates,"\n")

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

    num_dilutions = 4
    plate_vol = 7.5
    dilution_vol = 9
    waste_vol = 2.5
    plating_row = 0

    redo = True
    media_per_tube = 150
    plate = agar_plate_names[0]
    plate_counter = 0

    # Aliquot the LB into the PCR tube strip for the dilutions
    #p200.pick_up_tip()
    #for well in range(8):
    #    print("Transferring {}ul to tube {}".format(media_per_tube,well))
    #    p200.transfer(media_per_tube, centrifuge_tube['A1'].bottom(),master.wells(well).bottom(),new_tip='never')
    #p200.drop_tip()

    #input("Start the other run")
    # Iterate through each row of the transformation plate
    for trans_row in range(num_rows):
        if trans_row == 0:
            continue
        # Iterates through each dilution
        for dilution in range(num_dilutions):
            # Resets to switch to the next plate
            if plating_row == 12:
                p200.pick_up_tip()
                for well in range(8):
                    print("Transferring {}ul to tube {}".format(media_per_tube,well))
                    p200.transfer(media_per_tube, centrifuge_tube['B1'].bottom(),master.wells(well).bottom(),new_tip='never')
                p200.drop_tip()
                plate = agar_plate_names[1]
                print("changing to :", plate)
                plating_row = 0
                redo = True
            print("trans_row",trans_row, "dilution",dilution,"plate",plate,"plate row",plating_row)
            p10.pick_up_tip()
            print("Diluting cells in row {} with {}ul".format(trans_row, dilution_vol))
            p10.transfer(dilution_vol, master['A1'].bottom(), transformation_plate.rows(trans_row).bottom(),new_tip='never',mix_before=(1,9))
            print("Plating {}ul from transformation row {} onto {} in row {}".format(plate_vol,trans_row,plate,plating_row))
            p10.aspirate(plate_vol,transformation_plate.rows(trans_row).bottom())
            p10.dispense(plate_vol, agar_plates[plate].rows(plating_row).bottom())

            # Will continue to try and recalibrate until no change is made to the height
            if redo:
                redo = change_height(agar_plates[plate],agar_plates[plate].rows(plating_row)[0])
                p10.blow_out()

            print("Discard {}ul from transformation row {} into waste tube".format(waste_vol,trans_row))
            p10.aspirate(waste_vol,transformation_plate.rows(trans_row).bottom())
            p10.drop_tip()
            plating_row += 1
    stop = datetime.now()
    print(stop)
    runtime = stop - start
    print("Total runtime is: ", runtime)
    print("Rehoming")
    robot.home()

    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()
    return

if __name__ == "__main__":
    plate()







#
