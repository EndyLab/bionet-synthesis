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
#  3 p200  master  master source p10
#  2       dest    dest   source p10
#  1       trash          source p10
#

robot.connect(robot.get_serial_ports_list()[0])
robot.home()

p200_tipracks = [
    containers.load('tiprack-200ul', 'A3'),
]

p10_tipracks = [
    containers.load('tiprack-10ul', 'E2'),
]

p10s_tipracks = [
    containers.load('tiprack-10ul', 'E3')
]

trash = containers.load('point', 'B1', 'holywastedplasticbatman')
master = containers.load('PCR-strip-tall', 'C3')

dest_plates = [
    containers.load('96-PCR-tall', 'C2'),
    containers.load('96-PCR-tall', 'B2')
]

source_plates = {
    layout['D1']: containers.load('96-flat', 'D1'),
    layout['D2']: containers.load('96-flat', 'D2'),
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
num_reactions = len(plan)
num_rows = num_reactions // 8 + 1
all_wells = dest_plates[0].wells() + dest_plates[1].wells()

print("Building {} reactions in {} rows".format(num_reactions, num_rows))

p10.pick_up_tip()

for i in range(num_rows):
    p = i // 12
    r = i % 12

    print("Transferring master mix to plate {} row {}".format(p, r))
    p10.transfer(8, master['A1'], dest_plates[p].rows(r).bottom(), blow_out=True, touch_tip=True, new_tip='never')

p10.drop_tip()

# Add multiples of mastermix to plates with multiple fragments
p10s.pick_up_tip()
for i, construct in plan.iterrows():
    vol = 8 * (len(construct['Wells'].split(',')) - 1)
    if vol > 0:
        print("Adding {} to well {} for multifragment assembly".format(vol, i))
        p10s.transfer(vol, master['A1'], all_wells[int(i)].bottom(), blow_out=True, touch_tip=True, new_tip='never')
p10s.drop_tip()

# Move source DNA into dest mastermixes
for i,construct in plan.iterrows():
    print("Building gene {} {}".format(i, construct['Gene']))
    fragments = construct['Wells'].split(',')
    for fragment in fragments:
        plate, well = fragment.split('-')
        print("    Adding fragment from plate {} well {}".format(plate,well))
        p10s.transfer(2, source_plates[plate].wells(well).bottom(), all_wells[int(i)].bottom(), blow_out=True, touch_tip=True, mix_before=(3,5))
