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
import ot_functions as ot
from db_config import *

def resuspension(session,engine,target):
    ot.print_center("============ Beginning resuspension ============")

    # Initial Setup
    fmoles = 40

    # Load files
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    args = parser.parse_args()

    # Verify that the correct robot is being used
    if args.run:
        ot.check_robot()

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
                "Trough" : "B1"
            }
    ot.print_layout(locations)

    ot.print_center('...Calculating the volumes to resuspend with...')

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

    # print("total volume of water needed: {}uL".format(total))
    # num_tubes = math.ceil(total / 1000)
    # print("Prep {} tubes with 1.2mL".format(num_tubes))
    # input("Press Enter when you have added them to the tube rack")

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
    centrifuge_tube = containers.load('trough-12row',locations["Trough"])
    # centrifuge_tube = containers.load('tube-rack-2ml',locations["Tube_rack"])

    resuspend_plate = containers.load('96-flat', locations["PLATE HERE"])

    p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)

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
            # if current_vol - vol < 200:
            #     tube_count += 1
            #     current_vol = 1200
            # current_vol -= vol
            print("Adding {}ul to well {} with the p10".format(vol, well))
            #p10s.transfer(vol, centrifuge_tube[tubes[tube_count]].bottom(), resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
            p10s.transfer(vol, centrifuge_tube['A1'], resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
            # print("Currently {}ul in tube {}".format(current_vol,tubes[tube_count]))
            last_pipette = "p10s"
        else:
            if last_pipette == "p10s":
                print("Changing to p200")
                p10s.drop_tip()
                p200.pick_up_tip()
            elif last_pipette == "neither":
                p200.pick_up_tip()

            # Changes tubes of water when one gets low
            # if current_vol - vol < 100:
            #     tube_count += 1
            #     current_vol = 1200
            # current_vol -= vol
            print("Adding {}ul to well {} with the p200".format(vol, well))
            #p200.transfer(vol, centrifuge_tube[tubes[tube_count]].bottom(), resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
            p200.transfer(vol, centrifuge_tube['A1'], resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
            # print("currently {}ul in tube {}".format(current_vol,tubes[tube_count]))
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

    commit = int(ot.request_info("Commit changes (1-yes, 2-no): ",type='int'))
    if commit == 1:
        session.commit()

    ot.print_center('...Completed resuspension...')

    return

if __name__ == "__main__":
    session,engine = connect_db()
    # Present all available plates to resuspend
    print("Which plate would you like to resuspend:")
    plate_names = []
    for index,plate in enumerate(session.query(Plate).filter(Plate.plate_type == 'syn_plate').filter(Plate.resuspended != 'resuspended')):
        print("{}. {} - {}".format(index,plate.plate_id,plate.plate_name))
        plate_names.append(plate.plate_name)

    # Asks the user for a number corresponding to the plate they want to resuspend
    number = int(ot.request_info("Which plate: ",type='int'))
    target = session.query(Plate).filter(Plate.plate_name == plate_names[number]).first()
    print("Will resuspend plate ",target.plate_name)
    resuspension(session,engine,target)

#
