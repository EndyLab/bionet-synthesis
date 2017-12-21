## FOR GIT REPOSITORY -- Take in build map and tells OT how to plate the cells

from opentrons import robot, containers, instruments
import argparse
import sys

import numpy as np
import pandas as pd

import os
import glob
import re

## Take in required information

# Load files
parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
args = parser.parse_args()

# Get all of the plate maps and display them so you can choose one
counter = 0
plate_map_number = []

if glob.glob("../builds/build*.csv"):
    print("previous builds")
    for file in glob.glob("../builds/build*.csv"):
        if "bad" in file:
            continue
    counter += 1
    print("{}. {}".format(counter,file[10:]))
    plate_map_number.append(file)
else:
    print("No previous builds")

# Asks the user for a number corresponding to the plate they want to resuspend
number = input("Which file: ")
number = int(number) - 1

# Import the desired plate map
build_map = pd.read_csv(plate_map_number[number])
print(plate_map_number[number][10:18])
print(build_map)

num_reactions = len(build_map)

num_rows = num_reactions // 8

agar_plate_names = []
for row in range(num_rows):
    current_plate = plate_map_number[number][10:18] + "_p" + str(row + 1)
    agar_plate_names.append(current_plate)

print("You will need {} agar plates".format(len(agar_plate_names)))

## Setting up the OT-1 deck

AGAR_SLOTS = ['D2','B2','D1','D3']

layout = list(zip(agar_plate_names,AGAR_SLOTS[:len(agar_plate_names)]))

# Specify locations, note that locations are indexed by the spot in the array
## MAKE INTO A DICTIONARY
locations = np.array([["tiprack-200", "A3"],
                    ["tiprack-10_1", "E2"],
                    ["tiprack-10s", "E3"],
                    ["tiprack-10_2", "E1"],
                    ["trash", "D1"],
                    ["PCR-strip-tall", "C3"],
                    ["Tube Rack","B1"],
                    ["Transformation", "C2"]
                    ])

locations = np.append(locations,layout, axis=0)

# Make the dataframe to represent the OT-1 deck
deck = ['A1','B2','C3','D2','E1']
slots = pd.Series(deck)
columns = sorted(slots.str[0].unique())
rows = sorted(slots.str[1].unique(), reverse=True)
layout_table = pd.DataFrame(index=rows, columns=columns)
layout_table.fillna("", inplace=True)

# Fill in the data frame with the locations
for row,col in locations:
    layout_table.loc[col[1], col[0]] = row

# Displays the required plate map and waits to proceed
print()
print("Please arrange the plates in the following configuration:")
print()
print(layout_table)
print()
input("Press enter to continue")


## Initialize the OT-1

# Determine whether to simulate or run the protocol
if args.run:
    port = robot.get_serial_ports_list()[0]
    print("Connecting robot to port {}".format(port))
    robot.connect(port)
else:
    print("Simulating protcol run")
    robot.connect()

# Start up and declare components on the deck
robot.home()

p200_tipracks = [
    containers.load('tiprack-200ul', locations[0,1]),
]

p10_tipracks = [
    containers.load('tiprack-10ul', locations[1,1]),
    containers.load('tiprack-10ul', locations[3,1])

]

p10s_tipracks = [
    containers.load('tiprack-10ul', locations[2,1])
]

transformation_plate = containers.load('96-PCR-tall', locations[7,1])
trash = containers.load('point', locations[4,1], 'holywastedplasticbatman')
centrifuge_tube = containers.load('tube-rack-2ml',locations[6,1])
master = containers.load('PCR-strip-tall', locations[5,1])

agar_plates = {}
for plate, slot in layout:
    agar_plates[plate] = containers.load('96-deep-well', slot)
    print("agar_plates", agar_plates[plate])

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
    name='p10-8s'
)

p200 = instruments.Pipette(
    axis='b',
    max_volume=200,
    min_volume=20,
    tip_racks=p200_tipracks,
    trash_container=trash,
    channels=1,
    name='p200-1'
)

## Run the protocol

num_dilutions = 12
plate_vol = 7.5
dilution_vol = 7.5
num_rows = 2

for plate in agar_plates:
    p10.pick_up_tip()
    for row in range(num_dilutions):
        print("Plating {}ul onto {} in row {}".format(plate_vol,plate,row))
        p10.transfer(plate_vol, transformation_plate.rows(row), agar_plates[plate],new_tip='never',mix_before=(2,9))
        p10.drop_tip()
        p10.pick_up_tip()
        print("Diluting with {}ul".format(dilution_vol))
        p10.transfer(dilution_vol, master['A1'],transformation_plate.rows(row), new_tip='never',mix_before=(2,9))
