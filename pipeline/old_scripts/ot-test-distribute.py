from opentrons import robot, containers, instruments
import pandas as pd
import argparse

# Load our files
#parser = argparse.ArgumentParser(description="Run a DNA build on an Opentrons OT-1 robot.")
#parser.add_argument('-l', '--layout', required=True, help="A CSV file describing the layout of the sourcep plates.")
#parser.add_argument('-b', '--build-plan', required=True, help="A CSV file describing the build plan.")
#args = parser.parse_args()

#layout = pd.read_csv(args.layout) # Use plate-loc-resuspend-sm.csv
#layout = layout.set_index('Position').to_dict()['Name'] # Turn into a location->name dict

#plan = pd.read_csv(args.build_plan, usecols=['Plate','Well','Volume']) # Used testing/plate-yield.csv

# Configure the robot

#  Layout:
#    A     B       C      D      E
#  3 trash master  source source p10
#  2 p200  dest    source source p10
#  1 p200  dest    source source p10
#

# Connect to OT-1

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
]

trash = containers.load('point', 'A3', 'holywastedplasticbatman')
master = containers.load('PCR-strip-tall', 'C3')

source_plates = {
    layout['D2']: containers.load('96-flat', 'D2'),
}

dest_plates = [
    containers.load('96-flat', 'C2')
]

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
    tip_racks=p10_tipracks,
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

# Load plates to be resuspended

# Add water to all of the wells

p200.pick_up_tip()
p200.distribute(
    20,
    source_plates.wells('A1'),
    dest_plates.wells(['B2','B3','B4']),
    disposal_vol=10
)
#for i in range(0,len(plan)):
#    print("Resuspending well {} on plate {} with {}ul".format(plan['Well'][i], plan['Plate'][i], plan['Volume'][i]))
#
#    pipette.distribute(
#        plan['Volume'][i],
#        dest_plates.wells(plan['Well'][5]),
#        plan['Plate'][i].wells(plan['Well'][i]),
#        disposal_vol=10
#    )
    
print(robot.commands())
