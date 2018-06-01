from opentrons import robot, containers, instruments
import numpy as np
import pandas as pd
import getch
import shutil
import os
import sys
import cv2
import ot_functions as ot

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

p10,p10s,p200 = ot.initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash)

a = 1
acceptable_directions = ['w','s','a','d','W','S']
print('started')
print("Current step size is: {}".format(a))
while True:
    directions = {
        'w' : [a,0,0],
        's' : [-a,0,0],
        'a' : [0,a,0],
        'd' : [0,-a,0],
        'W' : [0,0,a],
        'S' : [0,0,-a],
    }
    c = getch.getch()
    if c == 'r':
        a = a * 2
        print("Current step size is: {}".format(a))
    elif c == 'f':
        a = a / 2
        print("Current step size is: {}".format(a))
    elif c in acceptable_directions:
        print(c,directions[c])
        coords = directions[c]
        p10.robot._driver.move(x=coords[0],y=coords[1],z=coords[2],mode="relative")
    elif c == 'x':
        print('Exit')
        break
    else:
        print('invalid input')







#
