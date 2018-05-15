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

from config import *
from db_config import *
session,engine = connect_db()

# ## Generate the SQL Database
# session = build_db()
def resuspension(target):
    print("\n============ Beginning resuspension ============\n")

    ## Set up initial paths
    PIPELINE_PATH = BASE_PATH + "/pipeline"
    DATA_PATH = BASE_PATH + "/data"

    # Initial Setup
    fmoles = 40

    # Load files
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
            sys.exit("Run . robot.sh while in the /opentrons/robots directory to change the robot")

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
                "PLATE HERE" : "B2",
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

    def calc_vol(amount,length,fmoles=40):
        return math.ceil(((((amount * 1000)/(660*length))*1000) / fmoles) * 2),fmoles

    total = 0
    for well in target.wells:
        length = len(well.fragments.seq)
        amount = well.syn_yield
        volume,conc = calc_vol(amount,length)
        well.volume = volume
        well.concentration = conc
        total += volume
        session.add(well)

    print("total volume of water needed: {}uL".format(total))
    num_tubes = math.ceil(total / 1000)
    print("Prep {} tubes with 1.2mL".format(num_tubes))

    input("Press Enter when you have added them to the tube rack")

    ## Initialize the OT-1

    # Determine whether to simulate or run the protocol
    if args.run:
        #port = robot.get_serial_ports_list()[0]
        port = os.environ["ROBOT_DEV"]
        print("Connecting robot to port {}".format(port))
        robot.connect(port)
    else:
        print("Simulating protcol run")
        robot.connect()

    start = datetime.now()
    print("Starting run at: ",start)

    # Start up and declare components on the deck
    robot.home()

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

    resuspend_plate = containers.load('96-flat', locations["PLATE HERE"])

    p10 = instruments.Pipette(
        axis='a',
        max_volume=10,
        min_volume=0.5,
        tip_racks=p10_tipracks,
        trash_container=trash,
        channels=8,
        name='p10-8'
    )

    p10s = instruments.Pipette(
        axis='a',
        max_volume=10,
        min_volume=0.5,
        tip_racks=p10s_tipracks,
        trash_container=trash,
        channels=1,
        name='p10-8s',
        aspirate_speed=800,
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
        aspirate_speed=800,
        dispense_speed=800
    )

    ## =============================================
    ## OT-1 PROTOCOL
    ## =============================================

    # Start timer
    start = datetime.now()
    print("Starting run at: ",start)


    tubes = dict({1:"A1",2:"B1",3:"C1",4:"D1",5:"A2",6:"B2",7:"C2",8:"D2",9:"A3",10:"B3",11:"C3",12:"D3",13:"A4",14:"B4",15:"C4"})
    tube_count = 1
    current_vol = 1200
    last_pipette = "neither"
    exclude_wells = []

    for target_well in target.wells:
        vol = target_well.volume
        well = target_well.address
        if well in exclude_wells:
            continue
        # Determine which pipette to use
        if vol < 20:
            # Makes sure that the robot doesn't pick up multiple sets of tips
            if last_pipette == "p200":
                print("Changing to p10s")
                p200.drop_tip()
                p10s.pick_up_tip()
            elif last_pipette == "neither":
                p10s.pick_up_tip()

            # Changes tubes of water when one gets low
            if current_vol - vol < 200:
                tube_count += 1
                current_vol = 1200
            current_vol -= vol
            print("Adding {}ul to well {} with the p10".format(vol, well))
            p10s.transfer(vol, centrifuge_tube[tubes[tube_count]].bottom(), resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
            print("Currently {}ul in tube {}".format(current_vol,tubes[tube_count]))
            last_pipette = "p10s"
        else:
            if last_pipette == "p10s":
                print("Changing to p200")
                p10s.drop_tip()
                p200.pick_up_tip()
            elif last_pipette == "neither":
                p200.pick_up_tip()

            # Changes tubes of water when one gets low
            if current_vol - vol < 100:
                tube_count += 1
                current_vol = 1200
            current_vol -= vol
            print("Adding {}ul to well {} with the p200".format(vol, well))
            p200.transfer(vol, centrifuge_tube[tubes[tube_count]].bottom(), resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
            print("currently {}ul in tube {}".format(current_vol,tubes[tube_count]))
            last_pipette = "p200"

    # Last drop tip
    if last_pipette == "p10s":
        p10s.drop_tip()
    elif last_pipette == "p200":
        p200.drop_tip()

    stop = datetime.now()
    print(stop)
    runtime = stop - start
    print("Total runtime is: ", runtime)
    robot.home()

    target.resuspend()
    session.add(target)

    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()
    return

if __name__ == "__main__":
    # Present all available plates to resuspend
    print("Which plate would you like to resuspend:")
    plate_names = []
    for index,plate in enumerate(session.query(Plate).filter(Plate.plate_type == 'syn_plate').filter(Plate.resuspended != 'resuspended')):
        print("{}. {}".format(index,plate.plate_name))
        plate_names.append(plate.plate_name)

    # Asks the user for a number corresponding to the plate they want to resuspend
    number = int(input("Which file: "))
    target = session.query(Plate).filter(Plate.plate_name == plate_names[number]).first()
    print("Will resuspend plate ",target.plate_name)
    resuspension(target)

#
