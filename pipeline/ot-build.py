from opentrons import robot, containers, instruments
import pandas as pd
import argparse

# Load our files
parser = argparse.ArgumentParser(description="Run a DNA build on an Opentrons OT-1 robot.")
parser.add_argument('-l', '--layout', required=True, help="A CSV file describing the layout of the sourcep plates.")
parser.add_argument('-b', '--build-plan', required=True, help="A CSV file describing the build plan.")
args = parser.parse_args()

layout = pd.read_csv(args.layout)
layout = layout.set_index('Position').to_dict()['Name'] # Turn into a location->name dict

plan = pd.read_csv(args.build_plan, usecols=['Gene','Wells'])

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

source_plates = {
    layout['C1']: containers.load('96-flat', 'C1'),
    layout['C2']: containers.load('96-flat', 'C2'),
    layout['C3']: containers.load('96-flat', 'C3'),
    layout['D1']: containers.load('96-flat', 'D1'),
    layout['D2']: containers.load('96-flat', 'D2'),
    layout['D3']: containers.load('96-flat', 'D3'),
}

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
