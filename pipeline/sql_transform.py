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

def transform():
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

    # Specify the locations of each object on the deck
    locations = {
                "tiprack-200" : "A3",
                "tiprack-10_2" : "E2",
                "tiprack-10_3" : "E3",
                "tiprack-10_1" : "E1",
                "trash" : "D1",
                "Transformation" : "C2",
                "Build_plate" : "C2",
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
    build_plate = containers.load('96-PCR-tall', locations['Build_plate'])

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

    for row in range(12):
        p10.pick_up_tip()
        p10.transfer(2, build_plate.rows(row).bottom(), transformation_plate.rows(trans_row).bottom(),new_tip='never',mix_before=(1,9))
        p10.drop_tip()








#
