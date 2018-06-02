from opentrons import robot, containers, instruments
import argparse

import numpy as np
import pandas as pd
import getch
import shutil
import os
import sys
import cv2
import ot_functions as ot

# Take in the 'run' argument from the command line
parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
args = parser.parse_args()

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
#dest_plate = containers.load('96-PCR-tall', locations["DEST_PLATE"])


p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)

pipette_dict = {'p10':p10,'p10s':p10s,'p200':p200}

if args.run:
    port = os.environ["ROBOT_DEV"]
    print("Connecting robot to port {}".format(port))
    robot.connect(port)
else:
    print("Simulating protcol run")
    robot.connect()


def move_ot(pipette):
    '''Allows for easy manual driving of the OT-One'''
    step = 1
    print("Initiating manual control of OT")
    print("'w'-'s' = back - forward\n'a'-'d' = left - right\n'W'-'S' = up - down")
    print("Use 'r' and 'f' to increase/decrease step size")
    print("Use 'x' to exit")
    print("Current step size is: {}".format(step))
    while True:
        directions = {
            'w' : [0,step,0],
            's' : [0,-step,0],
            'a' : [step,0,0],
            'd' : [-step,0,0],
            'W' : [0,0,step],
            'S' : [0,0,-step],
        }
        c = getch.getch()
        if c == 'r':
            step = step * 2
            print("Current step size is: {}".format(step))
        elif c == 'f':
            step = step / 2
            print("Current step size is: {}".format(step))
        elif c in directions.keys():
            print(c,directions[c])
            coords = directions[c]
            p10.robot._driver.move(x=coords[0],y=coords[1],z=coords[2],mode="relative")
        elif c == 'x':
            print('Exit')
            return
        else:
            print('invalid input')

def check_calibration(robot,pipette,container,sub_location):
    pipette.move_to(sub_location.top(z=20))
    ans = ot.request_info('Does it line up? (y/n): ',type='string')
    if ans != 'n':
        pipette.move_to(sub_location.top(z=5))
        ans = ot.request_info('Still good? (y/n): ',type='string')
        if ans != 'n':
            pipette.move_to(sub_location.bottom())
    else:
        ans = ot.request_info('Need to rehome? (y/n): ',type='string')
        if ans != 'n':
            robot.home()
            ans = ot.request_info('Try a little higher? (y/n): ',type='string')
            if ans != 'n':
                pipette.move_to(sub_location.top(z=40))

    print("Adjust calibration")
    move_ot(pipette)
    #pipette.calibrate_position((container,sub_location.from_center(x=0, y=0, z=-1,reference=container)))


layout = pd.read_csv('./all_deck_layout.csv')

#check_calibration(robot,pipette_dict['p10'],dest_plate,dest_plate.rows(0))

for i,row in layout.iterrows():
    container_type = row['container_type']
    slot = row['slot']
    current_container = containers.load(container_type,slot)
    p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)
    pipette_dict = {'p10':p10,'p10s':p10s,'p200':p200}
    pipette = pipette_dict[row['pipette']]

    print(row['pipette'])
    print(container_type)
    print(slot)
    print(current_container)
    input('Make sure {} is in slot {} {}'.format(container_type,slot,row['pipette']))
    #import ipdb; ipdb.set_trace()
    #check_calibration(robot,pipette,current_container,current_container.rows(0))
    if container_type == 'point':
        check_calibration(robot,pipette,current_container,current_container)
    elif row['pipette'] == 'p10':
        print('p10 elif')
        check_calibration(robot,pipette,current_container,current_container.rows(0))
    else:
        check_calibration(robot,pipette,current_container,current_container.wells(0))



#check_calibration(robot,p10,dest_plate,dest_plate.rows(1))







#
