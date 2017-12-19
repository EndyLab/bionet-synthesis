## FOR GIT REPOSITORY -- Queries the database for build ready constructs and then sets up the plates

from opentrons import robot, containers, instruments
import argparse
import sys

import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math
import datetime

## Take in required information

# Load files
parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
args = parser.parse_args()

# User input
#num_reactions = input("How many reactions: ")
#max_plates = input("Number of plates to pull from: ")
#max_frag = input("Max number of fragments to assemble: ")
#manual = input("Would you like to specify the plates: ")
#if manual != "":
#    plates = raw_input("Specify the plates: ")
#else:
#    plates = []

## Create a list of all well indexes to add in later
letter = ["A","B","C","D","E","F","G","H"]
number = ["1","2","3","4","5","6","7","8","9","10","11","12"]
target_well = []
temp_well = 0
for n in number:
    for l in letter:
        temp_well = l + n
        target_well.append(temp_well)

## Define initial conditions
targets = []
counter = 0
frag_list = []
master_well = []

max_plates = 2
plates = ["pSHPs0807B412037MU", "pSHPs0807B412038MU"]
#plates = []

gene_list = []

num_reactions = 20
max_frag = 2

# Query the database and iterate through each json file
for file in glob.glob("../data/BBF10K*/*.json"):
    print(file)

    # Open and store the data within the json file
    with open(file,"r") as json_file:
        data = json.load(json_file)

    # Determine if it is build ready
    if data["status"]["build_ready"] != "TRUE":
        continue
    print("build ready")

    #Set a limit for how many you want to run
    if len(gene_list) == num_reactions:
        break

    # Pull general information about the gene
    gene = data["gene_id"]
    locations = data["location"]["fragments"]
    frag_num = len(locations)

    # Set the limit on how many fragments you want to build in one reaction
    if frag_num > max_frag:
        continue
    print("number of fragments: ", frag_num)

    # Check if we'll need too many source plates if we add these gene to the build.
    new_plates = [loc.split("_")[0] for loc in locations.values()]
    new_plate_count = len(pd.unique(plates + new_plates))
    if new_plate_count > max_plates:
        # We'll have too many source plates if we add this gene.
        continue

    gene_list.append(gene)
    plates += new_plates
    dest_well = target_well[len(gene_list) - 1]

    print("number of unique plates:", len(pd.unique(plates)))

    # Iterate through all of the fragments
    for frag_count, (fragment, frag_loc) in enumerate(locations.items()):
        print(fragment)

        plate_loc, well = frag_loc.split("_")
        row = [gene, plate_loc, well, dest_well]
        targets.append(row)

        print()

    # If all of the fragments are cleared the number of reactions needed for the master mix is added along with the well it corresponds to
    frag_list.append(frag_num)
    print("num frag list: ", len(frag_list))
    master_well.append(dest_well)
    print("num wells: ", len(master_well))
    print()
    print()
print("unique plates: ", pd.unique(plates))

# Creates a dataframe with the well and
master_plan = pd.DataFrame({
    "Well" : master_well,
    "Fragments" : frag_list
})
master_plan = master_plan[["Well","Fragments"]]

targets = np.array(targets)
plan = pd.DataFrame({
    "Plate" : targets[:,1],
    "Gene" : targets[:,0],
    "Well" : targets[:,2],
    "Destination" : targets[:,3]
    })
plan = plan[["Gene","Plate","Well","Destination"]]
print(plan)
print(len(target_well))
print(master_plan)
print()
input("Press enter to continue")

## Setting up the OT-1 deck

# Configuration
SOURCE_SLOTS = ['D2','D3','D1','B2']


# Specify locations, note that locations are indexed by the spot in the array
locations = np.array([["tiprack-200", "A3"],
                    ["tiprack-10", "E2"],
                    ["tiprack-10s1", "E3"],
                    ["tiprack-10s2", "E1"],
                    ["trash", "B1"],
                    ["PCR-strip-tall", "C3"],
                    ["DEST_PLATE", "C2"],
                    ["Tube Rack","A1"]])

# Attach a location to each of the source plates
layout = list(zip(pd.unique(plates),SOURCE_SLOTS[:len(plates)]))
locations = np.append(locations,layout, axis=0)

# Make the dataframe to represent the OT-1 deck
deck = ['A1','B2','C3','D2','E1']
slots = pd.Series(deck)
columns = sorted(slots.str[0].unique())
rows = sorted(slots.str[1].unique(), reverse=True)
layout_table = pd.DataFrame(index=rows, columns=columns)
layout_table.fillna("---", inplace=True)

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

# Determine the number of rows to aliquot
num_rows = num_reactions // 8
print("Building {} reactions in {} rows".format(num_reactions, num_rows))

total_num = 0

for index, row in master_plan.iterrows():
    rxn_needed = int(row['Fragments'])
    total_num += rxn_needed

extra_master = 1.25

master_reactions = total_num * extra_master

print("You need {} rxns of master mix".format(master_reactions))

master_mix = pd.DataFrame({
            'Component':['Cutsmart','ATP','Vector','T4 Ligase','BbsI','H2O','Total'],
            'Amount':[master_reactions,master_reactions,(0.25*master_reactions),(master_reactions*(50/96)),(master_reactions*(6/96)),(5.166*master_reactions),(master_reactions*8)]
                        })
master_mix = master_mix[['Component','Amount']]
print("Use the table below to create the master mix")
print()
print(master_mix)
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
]

p10s_tipracks = [
    containers.load('tiprack-10ul', locations[2,1]),
    containers.load('tiprack-10ul', locations[3,1])
]

trash = containers.load('point', locations[4,1], 'holywastedplasticbatman')
centrifuge_tube = containers.load('tube-rack-2ml',locations[7,1])
master = containers.load('PCR-strip-tall', locations[5,1])

dest_plate = containers.load('96-PCR-tall', locations[6,1])

source_plates = {}
for plate, slot in layout:
    source_plates[plate] = containers.load('96-flat', slot)
    #print("source_plates", source_plates[plate])

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

## Aliquot the master mix into the PCR tube strip

vol_per_tube = num_rows * 8 * extra_master
vol_A1 = (total_num - num_reactions) * 8 * extra_master

print("{}ul into each PCR tube and {}ul extra into A1".format(vol_per_tube,vol_A1))

p200.pick_up_tip()
for well in range(8):
    print("Transferring {}ul to well {}".format(vol_per_tube,well))
    p200.transfer(vol_per_tube, centrifuge_tube['A1'].bottom(),master.wells(well).bottom(),touch_tip=True, mix_before=(3,50),new_tip='never')

print("Transferring an extra {}ul to well A1".format(vol_A1))
p200.transfer(vol_A1, centrifuge_tube['A1'].bottom(),master['A1'].bottom(), touch_tip=True, mix_before=(3,50),new_tip='never')
p200.drop_tip()

## Aliquot the master mix into all of the desired wells
p10.pick_up_tip()
for row in range(num_rows):
    print("Transferring master mix to row {}".format(row))
    p10.transfer(8, master['A1'].bottom(), dest_plate.rows(row).bottom(), touch_tip=True, mix_before=(3,8), new_tip='never')
p10.drop_tip()

p10s.pick_up_tip()
for index, row in master_plan.iterrows():
    if int(row['Fragments']) > 1:
        extra_volume = int(row['Fragments'] - 1) * 8
        current_well = str(row['Well'])
        print("Transferring {}ul of extra MM to {}".format(extra_volume,current_well))
        p10s.transfer(extra_volume,master['A1'].bottom(),dest_plate.wells(current_well).bottom(), touch_tip=True, mix_before=(3,8), new_tip='never')
p10s.drop_tip()

## Add the fragments from the source plates to the destination plate

# Sets the volume to pipet of each fragment to 2uL
frag_vol = 2

for index, row in plan.iterrows():
    plate = row['Plate']
    start_well = row['Well']
    dest_well = row['Destination']
    gene = row['Gene']
    print("Transferring {} from plate {} well {} to well {} of the dest plate".format(gene,plate,start_well,dest_well))
    p10s.pick_up_tip()
    p10s.transfer(frag_vol,source_plates[plate].wells(start_well).bottom(),dest_plate.wells(dest_well).bottom(), mix_before=(3,8), touch_tip=True)

# Record the current date and time
now, seconds = str(datetime.datetime.now()).split(".")

build_num = 0
if glob.glob("../builds/build*.csv"):
    print("previous builds")
    for file in glob.glob("../builds/build*.csv"):
        if "bad" in file:
            continue
        build_num = int(file[15]) + 1
else:
    print("no previous builds")
    build_num = '1'

outcome = int(input("1 = Good run, 2 = Bad run: "))
build_name = "build{}".format(build_num)

if outcome != 2:
    file_name = "../builds/{}_{}.csv".format(build_name,now)
else:
    file_name = "../builds/bad-{}_{}.csv".format(build_name, now)

plate_map = plan[["Gene","Destination"]]
plate_map = plate_map.drop_duplicates(subset=['Gene'])
plate_map.set_index("Gene")
print()
print(plate_map)

#plate_map.to_csv(file_name)

## Udate the json file of all of the attempted genes

for index, row in plate_map.iterrows():
    file = glob.glob("../data/{}/{}.json".format(row["Gene"],row["Gene"]))
    print(file)

    # Open and store the data within the json file
    with open(file,"r") as json_file:
        data = json.load(json_file)

    data["status"]["build_attempts"][build_name] = {
        "build_well" : row["Destination"],
        "build_date" : now,
        "build_outcome" : "In process"
    }
    data["status"]["building"] = True



print()
