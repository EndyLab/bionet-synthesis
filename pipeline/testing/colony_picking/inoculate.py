from opentrons import robot, containers, instruments
from opentrons.helpers import helpers
import argparse
import os

def mix_in_well(pipette,depth=-0.75,location=None,radius=0.7):

    well_edges = [
        location.from_center(x=radius, y=0, z=depth),       # right edge
        location.from_center(x=0, y=radius, z=depth),       # back edge
        location.from_center(x=radius * -1, y=0, z=depth),  # left edge
        location.from_center(x=0, y=radius * -1, z=depth)   # front edge
    ]

    [pipette.move_to((location, e), strategy='direct') for e in well_edges]

    return pipette

def inoculate(pipette,location=None,depth=-0.75,radius=0.7,mix=3):
    
    _description = 'Inoculating'
    pipette.robot.add_command(_description)
    
    if location:
        pipette.move_to(location, strategy='arc')
    else:
        location = pipette.previous_placeable
    
    for num in range(mix):
        mix_in_well(pipette,location=location,depth=depth,radius=radius)
    
    pipette.move_to((location,location.from_center(x=0, y=0, z=depth)),strategy='direct')

    pipette.motor.move(pipette._get_plunger_position('drop_tip'))
    pipette.motor.move(pipette._get_plunger_position('bottom'))

    pipette.current_volume = 0
    pipette.current_tip(None)

    return

if __name__ == "__main__":

    # Load files
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    args = parser.parse_args()

    if args.run:
        port = os.environ["ROBOT_DEV"]
        print("Connecting robot to port {}".format(port))
        robot.connect(port)
    else:
        print("Simulating protcol run")
        robot.connect()

    p10s_tipracks = [
        containers.load('tiprack-10ul', 'E2'),
    ]
    trash = containers.load('point', 'D1', 'holywastedplasticbatman')
    culture = containers.load('96-deep-well', 'C2')

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

    robot.home()
    
    for num in range(3):
        p10s.pick_up_tip()
        inoculate(p10s,location=culture.wells(num))








#
