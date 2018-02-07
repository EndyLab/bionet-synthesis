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
import math

import datetime
from datetime import datetime

#BASE_PATH = "/Users/conarymeyer/Desktop/GitHub/bionet-synthesis"
BASE_PATH = "../"
PIPELINE_PATH = BASE_PATH + "/pipeline"
BUILDS_PATH = BASE_PATH + "/builds"
DATA_PATH = BASE_PATH + "/data"

BACKBONE_PATH = BASE_PATH + "/sequencing_files/popen_v1-1_backbone.fasta"
DICTIONARY_PATH = PIPELINE_PATH + "/testing/data_testing/10K_CDS.csv"

# Initial Setup
fmoles = 40

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
for file in glob.glob("{}/plate_maps/*.csv".format(BASE_PATH)):
    counter = counter + 1
    plate_map_name = re.match(
        r'.+(Order[0-9]+_Attempt[0-9]+)_.+',
        file).groups()
    # Prints the file name without the preceeding path
    print("{}. {}".format(counter,plate_map_name[0]))
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
print("plan", plan)

total_vol = plan['Volume'].sum()
print("total volume of water needed: {}uL".format(total_vol))
num_tubes = math.ceil(total_vol / 1200)
print("Prep {} tubes with 1.2mL".format(num_tubes))

input("Add tubes of water")

## Setting up the OT-1 deck

# Specify locations, note that locations are indexed by the spot in the array
locations = np.array([["tiprack-200", "A3"],
                    ["tiprack-10", "E2"],
                    ["tiprack-10s1", "E3"],
                    ["tiprack-10s2", "E1"],
                    ["trash", "D1"],
                    ["PCR-strip-tall", "C3"],
                    ["PLATE HERE", "B2"],
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
    containers.load('tiprack-10ul', locations[3,1]),
    containers.load('tiprack-10ul', locations[2,1])
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

tubes = dict({1:"A1",2:"B1",3:"C1",4:"D1",5:"A2",6:"B2",7:"C2",8:"D2",9:"A3",10:"B3",11:"C3",12:"D3",13:"A4",14:"B4",15:"C4"})
tube_count = 1
current_vol = 1200
last_pipette = "neither"

exclude_wells = ["A2","B2","C2","D2","E2"]

## Run the protocol

# Iterate through each well
for i, construct in plan.iterrows():
    vol = construct['Volume']
    well = construct['Well']
    if well in exclude_wells:
        continue
    # Determine which pipette to use
    if vol < 20:
        print("Adding {}ul to well {} with the p10".format(vol, well))

        if last_pipette == "p200":
            p200.drop_tip()
            print("p200 drop")
            p10s.pick_up_tip()
            print("p10s pick")
        elif last_pipette == "neither":
            p10s.pick_up_tip()
            print("p10s pick")

        if current_vol - vol < 200:
            tube_count += 1
            current_vol = 1200
        current_vol -= vol
        print("currently {}ul in tube {}".format(current_vol,tubes[tube_count]))

        p10s.transfer(vol, centrifuge_tube[tubes[tube_count]].bottom(), resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
        last_pipette = "p10s"
    else:
        if last_pipette == "p10s":
            p10s.drop_tip()
            print("p10s drop")
            p200.pick_up_tip()
            print("p200 pick")
        elif last_pipette == "neither":
            p200.pick_up_tip()
            print("p200 pick")

        print("Adding {}ul to well {} with the p200".format(vol, well))

        if current_vol - vol < 100:
            tube_count += 1
            current_vol = 1200
        current_vol -= vol
        print("currently {}ul in tube {}".format(current_vol,tubes[tube_count]))

        p200.transfer(vol, centrifuge_tube[tubes[tube_count]].bottom(), resuspend_plate.wells(well),touch_tip=True, blow_out=True, new_tip='never')
        last_pipette = "p200"

if last_pipette == "p10s":
    p10s.drop_tip()
    print("p10s drop")
elif last_pipette == "p200":
    p200.drop_tip()
    print("p200 drop")

stop = datetime.now()
print(stop)
runtime = stop - start
print("Total runtime is: ", runtime)
