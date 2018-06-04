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

# Declare components on the deck
p200_tipracks = []
p10_tipracks = []
p10s_tipracks = []
trash = containers.load('point', 'D1', 'holywastedplasticbatman')

if args.run:
    port = os.environ["ROBOT_DEV"]
    ot.print_center("...Connecting robot to port {}...".format(port))
    robot.connect(port)
else:
    ot.print_center("...Simulating protcol run...")
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

def check_tips(pipette,sub_location):
    pipette.drop_tip(location=sub_location,home_after=False)
    pipette.pick_up_tip(location=sub_location)
    pipette.move_to(sub_location.top(z=30))
    ans = ot.request_info('Are tips sufficiently attached? (y/n): ',type='string')
    if ans == 'n':
        pipette.drop_tip(location=sub_location,home_after=False)
        move_ot(pipette)
        return check_tips(pipette,sub_location)
    else:
        return

def check_calibration(robot,pipette,container,height=20):
    if 's' in pipette.name or '200' in pipette.name:
        sub_location = container.wells(0)
    elif container.get_type() == 'point':
        sub_location = container
    else:
        sub_location = container.rows(0)
    pipette.move_to(sub_location.top(z=height),strategy='arc')
    ans = ot.request_info('Does it line up? (y/n): ',type='string')
    if ans != 'n':
        if 'tiprack' in str(container.get_type()):
            pipette.move_to(sub_location)
        else:
            pipette.move_to(sub_location.top())
        ans = ot.request_info('Still aligned? (y/n): ',type='string')
        if ans != 'n':
            pipette.move_to(sub_location.bottom())
        else:
            print('Drive OT from here')
    else:
        ans = ot.request_info('Need to rehome? (y/n): ',type='string')
        if ans != 'n':
            robot.home()
            ans = ot.request_info('Try a little higher? (y/n): ',type='string')
            if ans != 'n':
                return check_calibration(robot,pipette,container,height=height+20)
    print("Adjust {} calibration".format(container.get_type()))
    move_ot(pipette)
    if 'tiprack' in str(container.get_type()):
        print("Adjust tiprack calibration")
        check_tips(pipette,sub_location)
    pipette.calibrate_position((container,sub_location.from_center(x=0, y=0, z=-1,reference=container)))
    return sub_location

layout = pd.read_csv('./all_deck_layout.csv')
layout['type'] = layout.container_type.apply(lambda x: 0 if 'tiprack' in x else 1)
sections = layout.groupby(['pipette','type'])
for (pip,type),data in sections:
    data = data.reset_index()
    for i,row in data.iterrows():

        input('Make sure {} is in slot {}'.format(row['container_type'],row['slot']))
        current_container = containers.load(row['container_type'],row['slot'])

        p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)
        pipette_dict = {'p10':p10,'p10s':p10s,'p200':p200}
        pipette = pipette_dict[row['pipette']]
        sub_location = check_calibration(robot,pipette,current_container)
        if type == 0:
            print('Tipracks',len(data))

            if i == len(data)-1:
                print('Last rack')
                last_rack = sub_location
            else:
                pipette.drop_tip(sub_location)
                print('On to next rack')

        elif type == 1:
            if i == len(data)-1:
                print('Last container')
                pipette.drop_tip(last_rack)
            else:
                print('On to next container')
            # check_calibration(robot,pipette,current_container)


    # print(pip,'======',type)
    # print(data)
# print(layout)
# input('check')

#check_calibration(robot,pipette_dict['p10'],dest_plate,dest_plate.rows(0))

# for i,row in layout.iterrows():
#     container_type = row['container_type']
#     slot = row['slot']
#     current_container = containers.load(container_type,slot)
#     p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)
#     pipette_dict = {'p10':p10,'p10s':p10s,'p200':p200}
#     pipette = pipette_dict[row['pipette']]
#
#     print(row['pipette'])
#     print(container_type)
#     print(slot)
#     print(current_container)
#     input('Make sure {} is in slot {} {}'.format(container_type,slot,row['pipette']))
#     #import ipdb; ipdb.set_trace()
#     #check_calibration(robot,pipette,current_container,current_container.rows(0))
#     if container_type == 'point':
#         check_calibration(robot,pipette,current_container,current_container)
#     elif row['pipette'] == 'p10':
#         print('p10 elif')
#         check_calibration(robot,pipette,current_container,current_container.rows(0))
#     else:
#         check_calibration(robot,pipette,current_container,current_container.wells(0))



#check_calibration(robot,p10,dest_plate,dest_plate.rows(1))







#
