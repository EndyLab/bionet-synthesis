## FOR GIT REPOSITORY -- Take in plate map and tells OT volumes to resuspend with
## Read in a plate_map.csv, add locations, and assesses if it is buildable by iterating through the list

from opentrons import robot, containers, instruments
import argparse
import sys

import numpy as np
import pandas as pd

import os
import glob
import re

import datetime
from datetime import datetime

# Initial Setup
fmoles = 20

## Take in required information

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

# Get all of the plate maps and display them so you can choose one
counter = 0
plate_map_number = []
for file in glob.glob("../plate_maps/*.csv"):
    counter = counter + 1

    # Prints the file name without the preceeding path
    print("{}. {}".format(counter,file[14:]))
    plate_map_number.append(file)

# Asks the user for a number corresponding to the plate they want to resuspend
number = input("Which file: ")
number = int(number) - 1

# Import the desired plate map
plates = pd.read_csv(plate_map_number[number])

# Pulls out all of the unique plate names
plate_names = pd.Series(plates['Plate'].unique())
print(plate_names)

# Selects the desired plate to be resuspended
plate_num = input("Which plate to resuspend: ")
plate_num = int(plate_num)
plate = plate_names[plate_num]

# Truncates the dataframe to only include the desired plate
plates = plates.set_index('Plate')
plates = plates.loc[plate]

# Calculate the amount of water to resuspend each well with
amount = plates['Yield (ng)']
length = plates['synthesized sequence length']
volume = ((((amount * 1000)/(660*length))*1000) / fmoles) * 2
plan = pd.DataFrame({'Well': plates['Well'],
                     'Volume': volume})
plan = plan[['Well','Volume']]

## Setting up the OT-1 deck

# Specify locations, note that locations are indexed by the spot in the array
locations = np.array([["tiprack-200", "A3"],
                    ["tiprack-10", "E2"],
                    ["tiprack-10s1", "E3"],
                    ["tiprack-10s2", "E1"],
                    ["trash", "D1"],
                    ["PCR-strip-tall", "C3"],
                    ["PLATE HERE", "C2"],
                    ["Tube Rack","B1"]])

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
    containers.load('tiprack-200ul', locations[0,1]),
]

p10_tipracks = [
    containers.load('tiprack-10ul', locations[1,1]),
]

p10s_tipracks = [
    containers.load('tiprack-10ul', locations[2,1]),
    containers.load('tiprack-10ul', locations[3,1])
]

trash = containers.load('point', locations[4,1], 'holywastedplasticbatman')
master = containers.load('PCR-strip-tall', locations[5,1])
centrifuge_tube = containers.load('tube-rack-2ml',locations[7,1])

resuspend_plate = containers.load('96-flat', locations[6,1])

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

# Iterate through each well
for i, construct in plan.iterrows():
    vol = construct['Volume']
    well = construct['Well']

    # Determine which pipette to use
    if vol < 20:
        print("Adding {}ul to well {} with the p10".format(vol, well))
        p10s.pick_up_tip()
        p10s.transfer(vol, centrifuge_tube['A1'].bottom(), resuspend_plate.wells(well), blow_out=True, new_tip='never')
        #print("Mixing with {} 3 times".format(vol * 0.5))
        #p10s.mix(3, (vol * 0.5), resuspend_plate.wells(well).bottom())
        p10s.drop_tip()
    else:
        print("Adding {}ul to well {} with the p200".format(vol, well))
        p200.pick_up_tip()
        p200.transfer(vol, centrifuge_tube['A1'].bottom(), resuspend_plate.wells(well), blow_out=True, new_tip='never')
        #print("Mixing with {} 3 times".format(vol * 0.5))
        #p200.mix(3, (vol * 0.5), resuspend_plate.wells(well).bottom())
        p200.drop_tip()

stop = datetime.now()
print(stop)
runtime = stop - start
print("Total runtime is: ", runtime)
