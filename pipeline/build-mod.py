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
from datetime import datetime
import getch
from config import *

## Take in required information
## ============================================

# Set starting paths
# BASE_PATH is defined in the config file
PIPELINE_PATH = BASE_PATH + "/pipeline"
BUILDS_PATH = BASE_PATH + "/builds"
DATA_PATH = BASE_PATH + "/data"

# Load files
parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
parser.add_argument('-m', '--manual', required=False, action="store_true", help="Maunal entry of parameters.")
args = parser.parse_args()

plates = []
names = []
cand = []
plate = []
well = []
frag_nums = []
dest_wells = []

# Establish initial functions
def change_height(container,target):
    x = 0
    print("Change height - s-g:up h-l:down x:exit")
    while x == 0:
        c = getch.getch()
        if c == "s":
            p10.robot._driver.move(z=20,mode="relative")
        elif c == "d":
            p10.robot._driver.move(z=5,mode="relative")
        elif c == "f":
            p10.robot._driver.move(z=0.5,mode="relative")
        elif c == "g":
            p10.robot._driver.move(z=0.1,mode="relative")
        elif c == "h":
            p10.robot._driver.move(z=-0.1,mode="relative")
        elif c == "j":
            p10.robot._driver.move(z=-0.5,mode="relative")
        elif c == "k":
            p10.robot._driver.move(z=-5,mode="relative")
        elif c == "l":
            p10.robot._driver.move(z=-20,mode="relative")
        elif c == "x":
            x = 1
    p10s.calibrate_position((container,target.from_center(x=0, y=0, z=-1,reference=container)))

# Verify that the correct robot is being used
if args.run:
    robot_name = str(os.environ["ROBOT_DEV"][-5:])
    robot_number = str(input("Run on this robot: {} ? 1-Yes, 2-No ".format(robot_name)))
    if robot_number == "1":
        print("Proceeding with run")
    else:
        sys.exit("Run . robot.sh while in the /opentrons/robots directory to change the robot")

# Allows the user to specify different aspects of the run
if args.manual:
    max_reactions = input("How many reactions: ")
    max_plates = input("Number of plates to pull from: ")
    max_frag = input("Max number of fragments to assemble: ")
    manual = raw_input("Would you like to specify the plates: ")
    if manual != "":
        for num in range(max_plates):
            single_plate = raw_input("Specify plate: ")
            plates.append(single_plate)
    else:
        plates = []
else:
    max_plates = 3
    #plates = []

    print("Max plates: ", max_plates)
    print("Pulling from plates: ", plates)
    max_reactions = 96
    print("Number of reactions: ", max_reactions)
    max_frag = 4
    print("Max number of fragments: ", max_frag)
    input("Press enter to continue ")


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
gene_list = []

## Takes in all of the previously attempted genes
previous_genes = []
for file in glob.glob(BUILDS_PATH+"/*/*20*.csv"):
    build = pd.read_csv(file)
    previous_genes += list(build['Gene'])

## IF THERE ARE SPECIFIC GENES REQUIRED
## ============================================

#required_genes = pd.read_csv("./testing/select_plates.csv")
#required_genes = list(required_genes["Gene ID"])
#print(required_genes)

#completed_genes = []
#complete = False

#input("continue?")

# Get the specified set first:
#for required_gene in required_genes:
#    for file in glob.glob(DATA_PATH + "/{}/{}.json".format(required_gene,required_gene)):
#        print(file)
#        with open(file,"r") as json_file:
#            data = json.load(json_file)
#        gene = data["gene_id"]
#        locations = data["location"]["fragments"]
#        frag_num = len(locations)
#
#        new_plates = [loc.split("_")[0] for loc in locations.values()]
#        new_plate_count = len(pd.unique(plates + new_plates))
#
#        gene_list.append(gene)
#        plates += new_plates
#        dest_well = target_well[len(gene_list) - 1]

#        for frag_count, (fragment, frag_loc) in enumerate(locations.items()):
#            print(fragment)

#            plate_loc, well = frag_loc.split("_")
#            row = [gene, plate_loc, well, dest_well]
#            targets.append(row)

#        frag_list.append(frag_num)
#        master_well.append(dest_well)


#print("remaining_reactions",remaining_reactions)
#print("\n\n")


# FIND GENES TO BUILD AND THEIR CORRESPONDING FRAGMENTS
## ============================================

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv(BASE_PATH + "/pipeline/testing/data_testing/10K_CDS.csv")
dictionary = dict(zip(data['gene_name'], data['idnum']))

# Pulls in all of the previously attempted genes
previous_genes = []
for file in glob.glob(BASE_PATH + "/builds/*/build*_20*.csv"):
    #print(file)
    build = pd.read_csv(file)
    previous_genes += list(build['Gene'])

# Extracts all of the gene names that are currently buildable and have not been attempted
for file in glob.glob(BASE_PATH + "/data/*/*.json"):
    #print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    if data["status"]["build_ready"] == True and data["gene_id"] not in previous_genes:
        gene_name = data["gene_name"]
        names.append(gene_name)

print(previous_genes)
print("Names: ", names)
print(len(names))

input("continue?")

unknown_names = []

# Pulls all of the plate_maps and stores the locations of all of the flagged fragments
for plate_map in glob.glob(BASE_PATH + "/plate_maps/*.csv"):
    details = pd.read_csv(plate_map)
    print("plate_map", plate_map)
    for index,row in details.iterrows():
        if re.match(r'.+_[0-9]',row["customer_line_item_id"]):
            gene_name = row["customer_line_item_id"][:-2].strip()
        else:
            unknown_names.append(row["customer_line_item_id"])
            continue
        print(gene_name)
        if gene_name in names and gene_name not in previous_genes:
            print("match")
            cand.append(gene_name)
            plate.append(row["Plate"])
            well.append(row["Well"])

plan = pd.DataFrame({
        "Gene" : cand,
        "Plate" : plate,
        "Well" : well
    })

for index, row in plan.iterrows():
    for file in glob.glob(BASE_PATH + "/data/{}/{}.json".format(dictionary[row["Gene"]],dictionary[row["Gene"]])):
        frag_num = 0
        #print(file)
        with open(file,"r") as json_file:
            data = json.load(json_file)
        #print(data["location"]["fragments"])
        for fragment in data["location"]["fragments"]:
            frag_num += 1
        frag_nums.append(str(frag_num))

plan["Fragments"] = frag_nums

print("length", len(list(pd.unique(plan["Gene"]))))
num_wells = len(list(pd.unique(plan["Gene"])))
print("unique", pd.unique(plan["Gene"]))
print("target_well", target_well)
print("target_well_index", target_well[0])


well_dict = dict(zip(list(pd.unique(plan["Gene"])),target_well[:num_wells]))
for index,row in plan.iterrows():
    dest_wells.append(well_dict[row["Gene"]])

plan["Destination"] = dest_wells

plan = plan[["Gene","Plate","Well","Fragments","Destination"]]

# Creates a dataframe with the well and
master_plan = pd.DataFrame({
    "Well" : plan["Destination"],
    "Fragments" : plan["Fragments"]
    })
master_plan = master_plan[["Well","Fragments"]].drop_duplicates()
num_reactions = len(plan)


print(plan,"\n")
print("Total reactions: ",len(target_well),"\n")
print(master_plan,"\n")

unique_plates = list(pd.unique(plan["Plate"]))

group_plates = [unique_plates[n:n+3] for n in range(0, len(unique_plates), 3)]
print(group_plates)

for group in group_plates:
    print(group)

input("Press enter to continue")


#for plate in pd.unique(plates):
#    constructs_from_plate = len(plan[plan.Plate == plate])
#    print("Using {} fragments from plate {}".format(constructs_from_plate,plate))

#input("Press enter to continue")

## Setting up the OT-1 deck
## ============================================
SOURCE_SLOTS = ['D2','D3','B2']

# Specify locations, note that locations are indexed by the spot in the array
locations = np.array([["tiprack-200", "A3"],
                    ["tiprack-10", "E2"],
                    ["tiprack-10s1", "E3"],
                    ["tiprack-10s2", "E1"],
                    ["trash", "D1"],
                    ["PCR-strip-tall", "C3"],
                    ["DEST_PLATE", "C2"],
                    ["Tube Rack","B1"]])

# Attach a location to each of the source plates
#layout = list(zip(pd.unique(plates),SOURCE_SLOTS[:len(plates)]))
#locations = np.append(locations,layout, axis=0)

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
print("\n Please arrange the plates in the following configuration: \n")
print(layout_table,"\n")
input("Press enter to continue")


## SETUP THE MASTER MIX
## ============================================
# Determine the number of rows to aliquot
num_rows = num_reactions // 8
print("Building {} reactions in {} rows".format(num_reactions, num_rows))
total_num = 0
for index, row in master_plan.iterrows():
    rxn_needed = int(row['Fragments'])
    total_num += rxn_needed

## ONLY IF YOU ARE ONLY RUNNING TWO FRAGMENT REACTIONS
#total_num = total_num / 2

extra_master = 1.3
master_reactions = total_num * extra_master
print("You need {} rxns of master mix".format(master_reactions))

# REVIEW: FIX THE MASTER MIX CALCULATION
master_volume = 8

cutsmart = 1
atp = 1
vector = 0.25
ligase = 0.5
enzyme = 0.25
water = master_volume - (cutsmart + atp + vector + ligase + enzyme)

master_mix = pd.DataFrame({
            'Component':['Cutsmart','ATP','Vector','T4 Ligase','BbsI','H2O','Total'],
            'Amount':[(cutsmart*master_reactions),(atp*master_reactions),(vector*master_reactions),(ligase*master_reactions),(enzyme*master_reactions),(water*master_reactions),(master_volume*master_reactions)]
                        })
master_mix = master_mix[['Component','Amount']]
print("Use the table below to create the master mix")
print()
print(master_mix)
print()
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
centrifuge_tube = containers.load('tube-rack-2ml',locations[7,1])
master = containers.load('PCR-strip-tall', locations[5,1])

dest_plate = containers.load('96-PCR-tall', locations[6,1])

#source_plates = {}
#for plate, slot in layout:
#    source_plates[plate] = containers.load('96-flat', slot)

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

p10s = instruments.Pipette(
    axis='a',
    max_volume=10,
    min_volume=0.5,
    tip_racks=p10s_tipracks,
    trash_container=trash,
    channels=1,
    name='p10-8s',
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

## OT-ONE PROTOCOL
## ============================================
# Aliquot the master mix into the PCR tube strip

master_volume = 6
vol_per_tube = num_rows * master_volume * extra_master

print("{}ul into each PCR tube".format(vol_per_tube))

p200.pick_up_tip()
for well in range(8):
    print("Transferring {}ul to well {}".format(vol_per_tube,well))
    p200.transfer(vol_per_tube, centrifuge_tube['A1'].bottom(),master.wells(well).bottom(), mix_before=(3,50),new_tip='never')
p200.drop_tip()

# Aliquot the master mix into all of the desired wells
p10.pick_up_tip()
for row in range(num_rows):
    print("Transferring {}ul of master mix to row {}".format(master_volume,row))
    p10.transfer(master_volume, master['A1'].bottom(), dest_plate.rows(row).bottom(), mix_before=(1,8), blow_out=True, new_tip='never')
p10.drop_tip()

# Aliquot master mix into the last row if not a complete row
if num_reactions % 8 > 0:
    p10.pick_up_tip()
    print("need single channel for {}".format(num_reactions % 8))
    for missing in range(num_reactions % 8):
        extra_volume = 8
        current_well = (8 * num_rows) + (missing)
        print("Transferring {}ul of extra MM to {}".format(extra_volume,current_well))
        p10s.transfer(extra_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
    p10s.drop_tip()

# Aliquot extra master mix into wells with multiple fragments
p10s.pick_up_tip()
for index, row in master_plan.iterrows():
    if int(row['Fragments']) > 1:
        extra_volume = (int(row['Fragments']) - 1) * 8
        current_well = str(row['Well'])
        print("Transferring {}ul of extra MM to {}".format(extra_volume,current_well))
        p10s.transfer(extra_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
p10s.drop_tip()

## Add the fragments from the source plates to the destination plate
## ============================================
# Sets the volume to pipet of each fragment to 2uL
frag_vol = 2
dil_vol = 5
used_plates = []
plate_counter = 0
current_group = group_plates[plate_counter]
print(current_group)

def change_plates(current_plates):
    source_plates = {}
    plate_locations = list(zip(pd.unique(current_plates),SOURCE_SLOTS[:len(current_plates)]))
    print("plate_locations","\n", plate_locations)
    input("Switch out plates for those listed. Press enter when ready.")
    for plate, slot in plate_locations:
        source_plates[plate] = containers.load('96-flat', slot)
    return source_plates

source_plates = change_plates(current_group)

for index, row in plan.iterrows():
    plate = row['Plate']
    if plate not in current_group:
        plate_counter += 1
        current_group = group_plates[plate_counter]
        source_plates = change_plates(current_group)
    start_well = row['Well']
    dest_well = row['Destination']
    gene = row['Gene']
    p10s.pick_up_tip()

    ## ADD BACK IN IF VOLUME IN PLATES ARE LOW
    if plate == "pSHPs0807B412037MU" or plate == "pSHPs0807B412038MU" or plate == "pSHPs0807B412039MU":
        print("Diluting sample in plate {} well {} with {}uL of water".format(plate,start_well,dil_vol))
        p10s.transfer(dil_vol,centrifuge_tube['B1'].bottom(),source_plates[plate].wells(start_well).bottom(),new_tip='never')

    print("Transferring {} from plate {} well {} to well {} of the dest plate".format(gene,plate,start_well,dest_well))
    p10s.mix(3, 8, source_plates[plate].wells(start_well).bottom())
    if plate not in used_plates:
        change_height(source_plates[plate],source_plates[plate].wells(start_well))
    p10s.transfer(frag_vol,source_plates[plate].wells(start_well).bottom(),dest_plate.wells(dest_well).bottom(),blow_out=True)

    used_plates.append(plate)

# Record the current date and time
now, seconds = str(datetime.now()).split(".")

build_num = 0

## GENERATE THE BUILD MAP
## ============================================
# Assigns this build a unique number following the most recent build
max_num = 0
if glob.glob(BUILDS_PATH + "/*/build*_20*.csv"):
    for build_map in glob.glob(BUILDS_PATH + "/*/build*_20*.csv"):
        if "bad" in build_map:
            continue
        build_num = re.match(
            r'.+build([0-9]+)_20.+',
            build_map).groups()
        #current_build_num = str(int(build_num[0]) + 1).zfill(3)
        current_build_num = int(build_num[0]) + 1
        print("current_build_num", current_build_num)
        if int(current_build_num) > max_num:
            max_num = current_build_num
            build_str = str(max_num).zfill(3)
else:
    print("no previous builds")
    build_num = '001'

print("build number: ", build_str)

# Asks if it should record the results
outcome = int(input("1 = Good run, 2 = Bad run: "))
build_name = "build{}".format(build_str)

print("build_name: ", build_name)

if outcome != 2:
    file_name = BUILDS_PATH + "/{}/{}_{}.csv".format(build_name,build_name,now)
else:
    file_name = BUILDS_PATH + "/{}/bad-{}_{}.csv".format(build_name,build_name, now)

path = BUILDS_PATH + "/" + build_name

if os.path.exists(path):
    print("Directory {} already exists".format(build_name))
else:
    # Generates a new directory with the ID# as its name
    os.makedirs(path)


plate_map = plan[["Gene","Destination"]]
plate_map = plate_map.drop_duplicates(subset=['Gene'])
plate_map.set_index("Gene")
print("\n",plate_map)

# Export the build map
plate_map.to_csv(file_name)

## Update the json file of all of the attempted genes
## ============================================
for index, row in plate_map.iterrows():
    if outcome == 2:
        break
    gene = row["Gene"]
    for file in glob.glob(DATA_PATH + "/{}/{}.json".format(gene,gene)):
        print(file)

        # Open and store the data within the json file
        with open(file,"r") as json_file:
            data = json.load(json_file)

        # Adds a new set of build descriptions
        if data["status"]["build_attempts"][0]["build_well"] != "":
            data["status"]["build_attempts"].append({})

        attempt_num = len(data["status"]["build_attempts"]) - 1
        data["status"]["build_attempts"][attempt_num] = {}
        data["status"]["build_attempts"][attempt_num]["build_well"] = row["Destination"]
        data["status"]["build_attempts"][attempt_num]["build_date"] = now
        data["status"]["build_attempts"][attempt_num]["build_number"] = build_name
        data["status"]["build_attempts"][attempt_num]["build_outcome"] = "In Process"
        data["status"]["building"] = True

        with open(DATA_PATH + "/{}/{}.json".format(gene,gene),"w+") as json_file:
            json.dump(data,json_file,indent=2)

print()
stop = datetime.now()
print(stop)
runtime = stop - start
print("Total runtime is: ", runtime)
robot.home()
