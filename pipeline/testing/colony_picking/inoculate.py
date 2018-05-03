from opentrons import robot, containers, instruments
import argparse


def inoculate(pipette,location=None,radius=1.0):

    _description = 'Inoculating'
    pipette.robot.add_command(_description)

    height_offset = 0

    if helpers.is_number(location):
        height_offset = location
        location = None

    # if no location specified, use the previously
    # associated placeable to get Well dimensions
    if location:
        pipette.move_to(location, strategy='arc')
    else:
        location = pipette.previous_placeable

    v_offset = (0, 0, height_offset)

    well_edges = [
        location.from_center(x=radius, y=0, z=0.5),       # right edge
        location.from_center(x=radius * -1, y=0, z=0.5),  # left edge
        location.from_center(x=0, y=radius, z=0.5),       # back edge
        location.from_center(x=0, y=radius * -1, z=0.5)   # front edge
    ]

    # Apply vertical offset to well edges
    well_edges = map(lambda x: x + v_offset, well_edges)

    [pipette.move_to((location, e), strategy='direct') for e in well_edges]
    [pipette.move_to((location, e), strategy='direct') for e in well_edges]
    [pipette.move_to((location, e), strategy='direct') for e in well_edges]

    pipette.motor.move(pipette._get_plunger_position('drop_tip'))
    pipette.motor.move(pipette._get_plunger_position('bottom'))

    pipette.current_volume = 0
    pipette.current_tip(None)

    return pipette


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
    trash = containers.load('point', locations["trash"], 'holywastedplasticbatman')
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

    p10s.pick_up_tip()

    p10s.move_to(culture.wells('B1').bottom())

    p10s = inoculate(p10s)









#
