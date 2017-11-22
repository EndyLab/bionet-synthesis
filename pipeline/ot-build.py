from opentrons import robot, containers, instruments
import pandas as pd
import argparse
import sys

# Configuration
SOURCE_SLOTS = ['D1','D2','D3', 'E1','C2']

# Load our files
parser = argparse.ArgumentParser(description="Run a DNA build on an Opentrons OT-1 robot.")
parser.add_argument('-l', '--layout', required=False, help="A CSV file describing the layout of the sourcep plates.")
parser.add_argument('-b', '--build-plan', required=True, help="A CSV file describing the build plan.")
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
        print(slot[1], slot[0])
        layout_table.loc[slot[1], slot[0]] = plate

    print("Please arrange the plates in the following configuration:")
    print(layout_table)
    print()
    input("Press enter to continue")
    sys.exit(0)

    layout = dict(layout)

# Configure the robot

#  Layout:
#    A     B       C      D      E
#  3 trash master  source source p10
#  2 p200  dest    source source p10
#  1 p200  dest    source source p10
#

robot.connect() #robot.get_serial_ports_list()[0])

p200_tipracks = [
    containers.load('tiprack-200ul', 'A1'),
    containers.load('tiprack-200ul', 'A2')
]

p10_tipracks = [
    containers.load('tiprack-10ul', 'E1'),
    containers.load('tiprack-200ul', 'E2'),
    containers.load('tiprack-200ul', 'E3')
]

trash = containers.load('point', 'A3', 'trash')
master = containers.load('96-flat', 'B3', 'mastermix')

dest_plates = [
    containers.load('96-flat', 'B1'),
    containers.load('96-flat', 'B2')
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
num_reactions = len(plan)
num_rows = num_reactions // 8 + 1
all_wells = dest_plates[0].wells() + dest_plates[1].wells()

# for i in range(1,num_rows):
p10.transfer(8, master, all_wells[:num_reactions], blow_out=True, touch_tip=True)

# Add multiples of mastermix to plates with multiple fragments
# for i, construct in plan.iterrows():
#     vol = 8 * (len(construct['Wells'].split(',')) - 1)
#     if vol > 0:
#         p200.transfer(vol, master, all_wells[int(i)], blow_out=True, touch_tip=True)
#
# # Move source DNA into dest mastermixes
# for i,construct in plan.iterrows():
#     fragments = construct['Wells'].split(',')
#     for fragment in fragments:
#         plate, well = fragment.split('-')
#         p200.transfer(2, source_plates[plate].wells(well), all_wells[int(i)], blow_out=True, touch_tip=True, mix_before=(3,10), mix_after=(3,10))

print(robot.commands())
# robot.run()
