from opentrons import robot, containers, instruments
import argparse
import sys

import pandas as pd
import json
import glob
from datetime import datetime
import getch
import re

from config import *
import interact_db

## ============================================
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
args = parser.parse_args()

## ============================================
## ESTABLISH INITIAL FUNCTIONS
## ============================================
def change_height(container,target):
    '''Allows for real-time calibration of the p10 single height'''
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

def well_addresses():
    '''Generates a list of well address A1-H12'''
    letter = ["A","B","C","D","E","F","G","H"]
    number = ["1","2","3","4","5","6","7","8","9","10","11","12"]
    target_well = []
    temp_well = 0
    for n in number:
        for l in letter:
            temp_well = l + n
            target_well.append(temp_well)
    return target_well

def change_plates(current_plates):
    '''Allows the user to swap plates in the middle of a protocol'''
    source_plates = {}
    plate_locations = list(zip(pd.unique(current_plates),SOURCE_SLOTS[:len(current_plates)]))
    print("plate_locations","\n", plate_locations)
    input("Switch out plates for those listed. Press enter when ready.")
    for plate, slot in plate_locations:
        source_plates[plate] = containers.load('96-flat', slot)
    return source_plates

def make_gg_rxns(num_rxns,rxn_vol):
    '''Calculates the amount of each reagent to add to reach the desired master mix'''
    cutsmart = 1 * num_rxns
    atp = 1 * num_rxns
    vector = 0.25 * num_rxns
    ligase = 0.5 * num_rxns
    enzyme = 0.25 * num_rxns
    water = (rxn_vol - ((cutsmart + atp + vector + ligase + enzyme)/num_rxns)) * num_rxns
    master_mix = pd.DataFrame(
        {'Component':['H2O','Cutsmart','ATP','Vector','T4 Ligase','Restriction Enzyme','Total'],
        'Amount':[water,cutsmart,atp,vector,ligase,enzyme,rxn_vol*num_rxns]},
        columns=["Component","Amount"]
    )
    return master_mix

# Verify that the correct robot is being used
if args.run:
    robot_name = str(os.environ["ROBOT_DEV"][-5:])
    robot_number = str(input("Run on this robot: {} ? 1-Yes, 2-No ".format(robot_name)))
    if robot_number == "1":
        print("Proceeding with run")
    else:
        sys.exit("Run roboswitch 'ROBOT_NAME' to change to the correct robot")

target_well = well_addresses()
db = interact_db.query_db()

## =============================================
## CREATE LIST OF DESIRED GENES
## =============================================
# Pull out the desirable genes
build_db = db[['data.gene_id','data.status.build_ready','data.status.abandoned','data.status.build_complete','data.location.fragments']]
build_db = build_db[build_db['data.status.build_ready'] == True]
build_db = build_db[build_db['data.status.abandoned'] == False]
build_db = build_db[build_db['data.status.build_complete'] != "Good_sequence"]

# Remove the extra fragments
build_db = build_db.reset_index()
build_db = build_db.set_index(["data.gene_id"])
build_db = build_db[build_db['index'].str[-1] == "1"]

# Set the max number of reactions that you want to run
max_genes = 94

# Pull the subset of genes that will actually be built in this run
building = build_db.iloc[:max_genes]
building = building.reset_index()
gene_list = [g for g in building['data.gene_id']]

## =============================================
## CREATE BUILD PLAN
## =============================================
# Create the dataframe containing information for the build and
# for aliquoting master mix
gene_ids = []
frag_names = []
locations = []
dest_wells = []
master_wells = []
frag_nums = []
well_counter = 0
for gene in gene_list:
    with open("{}/data/{}/{}.json".format(BASE_PATH,gene,gene),"r") as json_file:
        data = json.load(json_file)
    dest_well = target_well[well_counter]
    frag_count = 0
    for fragment in data["location"]["fragments"]:
        gene_ids.append(data["gene_id"])
        frag_names.append(fragment)
        locations.append(data["location"]["fragments"][fragment])
        dest_wells.append(dest_well)
        frag_count += 1
    master_wells.append(dest_well)
    frag_nums.append(frag_count)
    well_counter += 1

plates = []
wells = []
for loc in locations:
    plate,well = loc.split("_")
    plates.append(plate)
    wells.append(well)

# Plan tells OT where to transfer the source DNA to
plan = pd.DataFrame({
    "Gene" : gene_ids,
    "Fragment" : frag_names,
    "Plate" : plates,
    "Well" : wells,
    "Destination" : dest_wells
},columns=["Gene","Fragment","Plate","Well","Destination"])

# Master tells OT how much master mix to add to each tube
master = pd.DataFrame({
    "Well" : master_wells,
    "Fragments" : frag_nums
},columns=["Well","Fragments"])
plan = plan.sort_values("Plate")

## =============================================
## CREATE PLATE GROUPS
## =============================================
# Currently available spots on the OT-one deck
SOURCE_SLOTS = ['D2','D3','B2']

# Group the unique plates by the number of available slots
unique_plates = list(pd.unique(plan["Plate"]))
group_plates = [unique_plates[n:n+len(SOURCE_SLOTS)] for n in range(0, len(unique_plates), len(SOURCE_SLOTS))]
for num,group in enumerate(group_plates):
    print("Group{}: {}".format(num+1,group))

input("Press enter to continue")

# Create a means to index to the plates
plate_dict = dict(enumerate(unique_plates))

# Present the plates that will be used to select if any need to be diluted
for num,plate in enumerate(unique_plates):
    print("{}: {}".format(num,plate))
print("{}: No plates".format(num+1))

# Add plates that need to be diluted into a list
more = True
dil_plates = []
plate_num = int(input("Choose the number of any plate that is low: "))
if plate_num == int(num)+1:
    print("No plates were added")
    more = False
else:
    dil_plates.append(plate_dict[int(plate_num)])
    while more:
        print("Plates to be diluted: {}".format(dil_plates))
        plate_num = int(input("Another plate? "))
        if plate_num == (int(num)+1):
            more = False
        else:
            dil_plates.append(plate_dict[int(plate_num)])
            print("No more plates were added")


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
            "PCR-strip-tall" : "C3",
            "DEST_PLATE" : "C2",
            "Tube_rack" : "B1"
        }

# Make the dataframe to represent the OT-1 deck
deck = ['A1','B2','C3','D2','E1']
slots = pd.Series(deck)
columns = sorted(slots.str[0].unique())
rows = sorted(slots.str[1].unique(), reverse=True)
layout_table = pd.DataFrame(index=rows, columns=columns)
layout_table.fillna("---", inplace=True)

# Fill in the data frame with the locations
for obj in locations:
        layout_table.loc[locations[obj][1], locations[obj][0]] = obj

# Displays the required plate map and waits to proceed
print("\n Please arrange the items in the following configuration: \n")
print(layout_table,"\n")
input("Press enter to continue")

## =============================================
## SETUP THE MASTER MIX
## =============================================
# Calculate how many reactions worth of MM to make
total_num = 0
for index, row in master.iterrows():
    rxn_needed = int(row['Fragments'])
    total_num += rxn_needed
num_wells = len(master)
print("wells",num_wells)
num_rows = num_wells // 8
print('rows',num_rows)

# Set a multiplier to account for pipetting error and evaporation
extra_master = 1.3
master_reactions = total_num * extra_master

# Set the volume of master mix to use per reaction
master_volume = 8

# Generate the dataframe to present the master mix composition
master_mix = make_gg_rxns(master_reactions,master_volume)
print("Use the table below to create the master mix")
print()
print(master_mix)
print()
input("Press enter to continue")

## =============================================
## INITIALIZE THE OT-1
## =============================================
# Determine whether to simulate or run the protocol
if args.run:
    port = os.environ["ROBOT_DEV"]
    print("Connecting robot to port {}".format(port))
    robot.connect(port)
else:
    print("Simulating protcol run")
    robot.connect()

# Declare components on the deck
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
centrifuge_tube = containers.load('tube-rack-2ml',locations["Tube Rack"])
master = containers.load('PCR-strip-tall', locations["PCR-strip-tall"])
dest_plate = containers.load('96-PCR-tall', locations["DEST_PLATE"])

# Declare all of the pipettes
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

## =============================================
## OT-1 PROTOCOL
## =============================================

# Start timer
start = datetime.now()
print("Starting run at: ",start)

# Home the robot to start
robot.home()

# Aliquot the master mix into the PCR tube strip
vol_per_tube = num_rows * master_volume * extra_master
print("Aliquoting MM into PCR tubes")
print("{}ul into each tube".format(vol_per_tube))
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
    p10s.pick_up_tip()
    print("need single channel for {}".format(num_reactions % 8))
    for missing in range(num_reactions % 8):
        current_well = (8 * num_rows) + (missing)
        print("Transferring {}ul of extra MM to {}".format(master_volume,current_well))
        p10s.transfer(master_volume, centrifuge_tube['A1'].bottom(),dest_plate.wells(current_well).bottom(),blow_out=True, mix_before=(1,8), new_tip='never')
    p10s.drop_tip()

# Aliquot extra master mix into wells with multiple fragments
p10s.pick_up_tip()
for index, row in master.iterrows():
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

# Sets the volume of water to dilute with, if needed
dil_vol = 5

# Sets the first group of plates
used_plates = []
plate_counter = 0
current_group = group_plates[plate_counter]
source_plates = change_plates(current_group)

# Trasfer DNA based on each item in plan
for index, row in plan.iterrows():
    plate = row['Plate']
    # Switches to the next group when the plate is not longer in it
    if plate not in current_group:
        plate_counter += 1
        current_group = group_plates[plate_counter]
        source_plates = change_plates(current_group)
    start_well = row['Well']
    dest_well = row['Destination']
    gene = row['Gene']

    p10s.pick_up_tip()

    # Only dilutes wells in plates marked to have low volumes
    if plate in dil_plates:
        print("Diluting sample in plate {} well {} with {}uL of water".format(plate,start_well,dil_vol))
        p10s.transfer(dil_vol,centrifuge_tube['B1'].bottom(),source_plates[plate].wells(start_well).bottom(),new_tip='never')

    print("Transferring {} from plate {} well {} to well {} of the dest plate".format(gene,plate,start_well,dest_well))
    p10s.mix(3, 8, source_plates[plate].wells(start_well).bottom())

    # Checks the calibration to make sure that it can aspirate correctly
    p10s.aspirate(frag_vol,source_plates[plate].wells(start_well).bottom())
    if plate not in used_plates:
        change_height(source_plates[plate],source_plates[plate].wells(start_well))
    p10s.dispense(frag_vol,dest_plate.wells(dest_well).bottom(),blow_out=True)
    used_plates.append(plate)


robot.home()

# Record the time that the protocol ended
now, seconds = str(datetime.now()).split(".")

## =============================================
## GENERATE THE BUILD MAP
## =============================================

previous_builds = (sorted(glob.glob(BASE_PATH + "/builds/*_20*.csv")))
last = previous_builds[-1]
build_name = "build{}".format(str(int(last[-3:])+1).zfill(3))
print("build number: ",build_name)

path_name = "{}/{}".format(BUILDS_PATH,build_name)
file_name = "{}/{}_{}.csv".format(path_name,build_name,now)

if os.path.exists(path):
    print("Directory {} already exists".format(build_name))
else:
    # Generates a new directory with the ID# as its name
    os.makedirs(path)

# Pull the relevant information from the build plan
plate_map = plan[["Gene","Destination"]]
plate_map = plate_map.drop_duplicates(subset=['Gene'])
plate_map = plate_map.sort_index()

# Fill in the remaining wells with 'Empty'
last_well = plate_map.iloc[-1]["Destination"]
remaining_wells = target_well[(target_well.index(last_well)+1):]
empty = ["Empty" for num in range(len(remaining_wells))]
empty = pd.DataFrame({
    "Gene" : empty,
    "Destination" : remaining_wells
})
plate_map = pd.concat([plate_map,empty])

# Export the build map
plate_map.to_csv(file_name,index=False)

## =============================================
## UPDATE THE DATABASE
## =============================================

for index, row in plate_map.iterrows():
    gene = row["Gene ID"]
    for file in glob.glob("{}/{}/{}.json".format(DATA_PATH,gene,gene)):
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












#
