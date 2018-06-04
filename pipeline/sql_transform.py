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
import ot_functions as ot
from db_config import *

def transform(session,engine):
    # Take in command line arguments
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    # parser.add_argument('-m', '--manual', required=False, action="store_true", help="Maunal entry of parameters.")
    args = parser.parse_args()

    assemblies = []
    print("Choose which plate you would like to transform/plate:")
    for index,assembly in enumerate(session.query(Plate).join(Build,Plate.builds).filter(Build.status == 'building').order_by(Build.build_name)):
        print("{}. {}".format(index,assembly.builds.build_name))
        assemblies.append(assembly)

    plate_num = int(ot.request_info("Enter plate here: ",type='int'))
    target_plate = assemblies[plate_num]

    num_reactions = len(target_plate.wells)
    num_rows = num_reactions // 8

    # Verify that the correct robot is being used
    if args.run:
        ot.check_robot()

    # Specify the locations of each object on the deck
    locations = {
                "tiprack-200" : "A3",
                "tiprack-10_2" : "E2",
                "tiprack-10_3" : "E3",
                "tiprack-10_1" : "E1",
                "trash" : "D1",
                "Transformation" : "C2",
                "Build_plate" : "C3",
                "Tube_rack" : "B1"
            }
    ot.print_layout(locations)

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

    ot.change_speed(robot)

    # Start up and declare components on the deck
    robot.home()

    p200_tipracks = [
        containers.load('tiprack-200ul', locations['tiprack-200']),
    ]

    p10_tipracks = [
        containers.load('tiprack-10ul', locations['tiprack-10_2']),
        containers.load('tiprack-10ul', locations['tiprack-10_1']),
    ]

    p10s_tipracks = [
        containers.load('tiprack-10ul', locations["tiprack-10_3"])
    ]

    transformation_plate = containers.load('96-PCR-tall', locations['Transformation'])
    trash = containers.load('point', locations['trash'], 'holywastedplasticbatman')
    centrifuge_tube = containers.load('tube-rack-2ml',locations['Tube_rack'])
    build_plate = containers.load('96-PCR-tall', locations['Build_plate'])

    p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)


    dna_vol = 2

    for row in range(num_rows):
        p10.pick_up_tip()
        p10.transfer(dna_vol, build_plate.rows(row).bottom(), transformation_plate.rows(row).bottom(),new_tip='never',mix_before=(2,9),blow_out=True)
        print('Transferring DNA from row {}'.format(row))
        p10.drop_tip()
    if num_reactions % 8 > 0:
        p10s.pick_up_tip()
        print("need single channel for {}".format(num_reactions % 8))
        for missing in range(num_reactions % 8):
            current_well = (8 * num_rows) + (missing)
            print("Transferring {}ul of DNA to {}".format(dna_vol,current_well))
            p10s.transfer(dna_vol, build_plate.wells(current_well).bottom(), transformation_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
            p10s.drop_tip()



if __name__ == '__main__':
    session,engine = connect_db()
    transform(session,engine)





#
