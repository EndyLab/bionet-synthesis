## FOR GIT REPOSITORY -- Take in plate map and tells OT volumes to resuspend with
## Read in a plate_map.csv, add locations, and assesses if it is buildable by iterating through the list

from opentrons import robot, containers, instruments
import argparse
import sys

import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns

#import json
import os
import glob
import re
#import math

## Take in required information

# Load files
parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
args = parser.parse_args()

## Get all of the plate maps and display them so you can choose one
counter = 0
plate_map_number = []
for file in glob.glob("../Plate Maps/platemaps/*.csv"):
    counter = counter + 1
    print("{}. {}".format(counter,file[24:]))
    plate_map_number.append(file)

## Asks the user for a number corresponding to the plate they want to resuspend
number = input("Which plate do you want to resuspend: ")
number = int(number) - 1
print("Resuspending plate",plate_map_number[number][24:])

## Import the desired plate map
plates = pd.read_csv(plate_map_number[number])

## Calculate the amount of water to resuspend each well with
amount = plates['Yield (ng)']
length = plates['synthesized sequence length']
fmoles = 40
volume = ((((amount * 1000)/(660*length))*1000) / fmoles) * 2
plan = pd.DataFrame({'Plate': plates['Plate'],
                     'Well': plates['Well'],
                     'Volume': volume})
plan = plan[['Plate','Well','Volume']]
print(plan)

## Setting up the OT-1 deck

robot.home()

p200_tipracks = [
    containers.load('tiprack-200ul', 'A3'),
]

p10_tipracks = [
    containers.load('tiprack-10ul', 'E2'),
]

p10s_tipracks = [
    containers.load('tiprack-10ul', 'E3'),
    containers.load('tiprack-10ul', 'E1'),
]

trash = containers.load('point', 'B1', 'holywastedplasticbatman')
master = containers.load('PCR-strip-tall', 'C3')

dest_plates = [
    containers.load('96-PCR-tall', 'C2'),
    containers.load('96-PCR-tall', 'B2')
]

source_plates = {}
for slot, plate in layout.items():
    source_plates[plate] = containers.load('96-flat', slot)





## Setup the robot
if args.run:
    port = robot.get_serial_ports_list()[0]
    print("Connecting robot to port {}".format(port))
    robot.connect(port)
else:
    print("Simulating protcol run")
    robot.connect()





print()
