## FOR GIT REPOSITORY -- Take in build map and tells OT how to plate the cells

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

## Take in required information

BASE_PATH = "/Users/conarymeyer/Desktop/GitHub/bionet-synthesis"
PIPELINE_PATH = BASE_PATH + "/pipeline"
BUILDS_PATH = BASE_PATH + "/builds"
DATA_PATH = BASE_PATH + "/data"

BACKBONE_PATH = BASE_PATH + "/sequencing_files/popen_v1-1_backbone.fasta"
DICTIONARY_PATH = PIPELINE_PATH + "/testing/data_testing/10K_CDS.csv"

# Load files
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

# Get all of the plate maps and display them so you can choose one
counter = 0
plate_map_number = []
build_map = 0

# Looks for the possible build maps to plate
for build_map in glob.glob("{}/*/build*.csv".format(BUILDS_PATH)):
    counter += 1
    print("{}. {}".format(counter,build_map[10:]))
    plate_map_number.append(build_map)

# Quits the program if there are no build maps to plate
if build_map == 0:
    print("No plate maps")
    sys.exit()

# Asks the user for a number corresponding to the plate they want to resuspend
number = input("Which file: ")
number = int(number) - 1

# Import the desired plate map
build_map = pd.read_csv(plate_map_number[number])
print(plate_map_number[number][10:18])
print(build_map)

num_reactions = len(build_map)
print("num reactions",num_reactions)

num_rows = num_reactions // 8
print(num_rows)
trans_per_plate = 3
num_plates = num_rows // trans_per_plate
print(num_plates)

agar_plate_names = []
for row in range(num_plates):
    current_plate = plate_map_number[number][10:18] + "_p" + str(row + 1)
    agar_plate_names.append(current_plate)
print(agar_plate_names)
print("You will need {} agar plates".format(len(agar_plate_names)))

## Setting up the OT-1 deck

AGAR_SLOTS = ['D2','B2','D3']

layout = list(zip(agar_plate_names,AGAR_SLOTS[:len(agar_plate_names)]))

# Specify locations, note that locations are indexed by the spot in the array
## MAKE INTO A DICTIONARY
locations = np.array([["tiprack-200", "A3"],
                    ["tiprack-10_2", "E2"],
                    ["tiprack-10_1", "E1"],
                    ["tiprack-10_3", "E3"],
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
    containers.load('tiprack-10ul', locations[2,1]),
    containers.load('tiprack-10ul', locations[3,1])
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

if args.manual:
    num_dilutions = input("Number of dilutions: ")
    plate_vol = input("Volume to plate: ")
    number  = input("Number of dilutions: ")
else:
    num_dilutions = 5
    plate_vol = 7.5
    dilution_vol = 9
    waste_vol = 2.5
    plating_row = 0

p10.pick_up_tip()

# Iterate through each agar plate
for plate in agar_plates:
    # Iterate through each row of cells in the transformation plate
    for trans_row in range(num_rows):
        #Reset the dilution factor and starting volume in the tubes
        dil_factor = [1]
        tube_vol = 15

        for plating_counter in range(num_dilutions):
            if plating_counter == 0:
                print("Discard {}ul from transformation row {} into waste tube".format(plate_vol,trans_row))
                p10.transfer(plate_vol, transformation_plate.rows(trans_row).bottom(), master['A12'].bottom(),new_tip='never',mix_before=(2,9))
                tube_vol -= plate_vol
                print("Volume after initial discard: ",tube_vol)
                plating_row -= 1
            else:
                print("Plating {}ul from transformation row {} onto {} in row {}".format(plate_vol,trans_row,plate,plating_row))
                p10.transfer(plate_vol, transformation_plate.rows(trans_row).bottom(), agar_plates[plate].rows(plating_row).bottom(),new_tip='never',mix_before=(2,9))
                print("Counter {} starts with {}uL in the transformation tube".format(plating_counter,tube_vol))
                tube_vol -= plate_vol
                print("After plating the volume is", tube_vol)
            if plating_counter != 0 and plating_counter != num_dilutions-1:
                print("Discard {}ul from transformation row {} into waste tube".format(waste_vol,trans_row))
                p10.aspirate(waste_vol,transformation_plate.rows(trans_row).bottom())
                #p10.transfer(waste_vol, transformation_plate.rows(trans_row).bottom(), master['A12'].bottom(),new_tip='never',mix_before=(2,9))
                tube_vol -= waste_vol
                print("Volume after discard: ", tube_vol)
            p10.drop_tip()
            print("Drop tip")
            if plating_counter == num_dilutions-1:
                plating_row += 1
                print("Last dilution in series, moving to next row")
                if trans_row != num_rows-1:
                    p10.pick_up_tip()
                    print("Pick up tip")
                continue
            print("Pick up tip")
            p10.pick_up_tip()
            print("Diluting cells in row {} with {}ul".format(trans_row, dilution_vol))
            p10.transfer(dilution_vol, master['A1'].bottom(), transformation_plate.rows(trans_row).bottom(),new_tip='never',mix_before=(2,9))
            previous_vol = tube_vol
            print("previous volume: ",previous_vol)
            tube_vol += dilution_vol
            print("current volume: ", tube_vol)
            dil_factor.append((tube_vol / previous_vol) * dil_factor[-1])
            print("dilution factors: ",dil_factor)

            print(tube_vol)
            plating_row += 1

        print("transformation row:",trans_row," --- dilution factors: ",dil_factor)
    plating_row = 0
    print("plating_row reset to 0")



stop = datetime.now()
print(stop)
runtime = stop - start
print("Total runtime is: ", runtime)
