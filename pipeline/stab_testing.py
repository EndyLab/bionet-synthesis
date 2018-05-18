from opentrons import robot, containers, instruments
import argparse
import sys

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




parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
parser.add_argument('-m', '--manual', required=False, action="store_true", help="Maunal entry of parameters.")
args = parser.parse_args()



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
    containers.load('tiprack-10ul', 'E3'),
]

transformation_plate = containers.load('96-PCR-tall', 'C2')
trash = containers.load('point', 'D1', 'holywastedplasticbatman')
master = containers.load('PCR-strip-tall', 'C3')
agar_plate = containers.load('96-deep-well', 'D2')

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

p10.move_to((agar_plate.rows(0),[0,0,1]),strategy='arc')
p10.move_to((agar_plate.rows(0),[0,0,0]),strategy='arc')

pipette.dispense(volume)
if calibrate:
    calibrate,z = change_height(agar_plates[plate],agar_plates[plate].rows(plating_row)[0])
pipette.move_to((row,[0,0,z]),strategy='direct')
input("Successful plating?")
return calibrate,z










#
