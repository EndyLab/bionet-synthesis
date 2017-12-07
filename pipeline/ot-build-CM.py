from opentrons import robot, containers, instruments
import pandas as pd
import argparse
import sys

# Configuration
SOURCE_SLOTS = ['D2','D3','D1','B2']

# Load our files
parser = argparse.ArgumentParser(description="Run a DNA build on an Opentrons OT-1 robot.")
parser.add_argument('-l', '--layout', required=False, help="A CSV file describing the layout of the sourcep plates.")
parser.add_argument('-b', '--build-plan', required=True, help="A CSV file describing the build plan.")
parser.add_argument('-s', '--simulate', required=False, action="store_true", help="Simulate the robot run and print command output.")
args = parser.parse_args()

plan = pd.read_csv(args.build_plan, usecols=['Gene','Wells'])



if args.layout:
    # We were given an explicit layout
    layout = pd.read_csv(args.layout)
    layout = layout.set_index('Position').to_dict()['Name'] # Turn into a location->name dict
else:
    plate_names = plan['Wells'].str.split('-').str[0].unique()

    if len(plate_names) > len(SOURCE_SLOTS):
        print("Error: This build plan requires too many source plates.")
        sys.exit(1)

    layout = list(zip(SOURCE_SLOTS[:len(plate_names)], plate_names))

    slots = pd.Series(SOURCE_SLOTS)
    columns = sorted(slots.str[0].unique())
    rows = sorted(slots.str[1].unique(), reverse=True)

    layout_table = pd.DataFrame(index=rows, columns=columns)
    layout_table.fillna("", inplace=True)

    for slot, plate in layout:
        layout_table.loc[slot[1], slot[0]] = plate

    print("Please arrange the plates in the following configuration:")
    print()
    print(layout_table)
    print()
    input("Press enter to continue")

    layout = dict(layout)


# Configure the robot

#  Layout:
#    A     B       C      D      E
#  3 p200  master  master source p10
#  2       dest    dest   source p10
#  1       trash          source p10
#

if args.simulate:
    print("Simulating protcol run")
    robot.connect()
else:
    port = robot.get_serial_ports_list()[0]
    print("Connecting robot to port {}".format(port))
    robot.connect(port)

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

# Run the protocol

# Load dest plates

# Distribute down the plates

max_reactions = 96

if len(plan) <= max_reactions:
    num_reactions = len(plan)
    print("Will run {} reactions".format(num_reactions))
else:
    num_reactions = max_reactions
    print("Too many reactions, will only run first {} reactions".format(num_reactions))

#num_reactions = len(plan)
num_rows = num_reactions // 8
print(num_rows)
all_wells = dest_plates[0].wells() + dest_plates[1].wells()

print("Building {} reactions in {} rows".format(num_reactions, num_rows))

p10.pick_up_tip()

for i in range(num_rows):
    p = i // 12
    r = i % 12

    print("Transferring master mix to plate {} row {}".format(p, r))
    p10.transfer(8, master['A1'], dest_plates[p].rows(r).bottom(), blow_out=True, touch_tip=True, new_tip='never')

p10.drop_tip()

j = 0

# Add multiples of mastermix to plates with multiple fragments
p10s.pick_up_tip()
for i, construct in plan.iterrows():
    vol = 8 * (len(construct['Wells'].split(',')) - 1)
    j = j + (len(construct['Wells'].split(',')) - 1)
    if vol > 0:
        print("Adding {} ul to well {} for multifragment assembly".format(vol, i))
        p10s.transfer(vol, master['A1'], all_wells[int(i)].bottom(), blow_out=True, touch_tip=True, new_tip='never')
    if i == num_reactions-1:
        break
p10s.drop_tip()

master_reactions = (num_reactions + j) * 1.1
print("Make {} rxns of master mix".format(master_reactions))

master_mix = pd.DataFrame({
        'Component':['Cutsmart','ATP','Vector','T4 Ligase','BbsI','H2O','Total'],
        'Amount':[master_reactions,master_reactions,(0.25*master_reactions),(master_reactions*(50/96)),(master_reactions*(6/96)),(5.166*master_reactions),(master_reactions*8)]
            })
master_mix = master_mix[['Component','Amount']]
print(master_mix)

input("Press enter to continue")


# Move source DNA into dest mastermixes
#for i,construct in plan.iterrows():
#    print("Building gene {} {}".format(i, construct['Gene']))
#    fragments = construct['Wells'].split(',')
#    for fragment in fragments:
#        plate, well = fragment.split('-')
#        print("    Adding fragment from plate {} well {}".format(plate,well))
#        p10s.transfer(2, source_plates[plate].wells(well).bottom(), all_wells[int(i)].bottom(), blow_out=True, touch_tip=True, mix_before=(3,5))
#    if i == num_reactions-1:
#        break
#if args.simulate:
#    print()
#    print("Ran commands:")
#    for c in robot.commands():
#        print(c)

plan[:max_reactions].to_csv('../synth1/round1-assemblies.csv')
