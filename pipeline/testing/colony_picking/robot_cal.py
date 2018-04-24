from matplotlib import pyplot as plt
from scipy.signal import butter, lfilter, freqz, argrelextrema

from imutils import contours
from skimage import measure
import numpy as np
import pandas as pd
import itertools
import argparse
import imutils
import cv2
import glob
import os
import getch
import math

from opentrons import robot, containers, instruments
import sys

from cam_calibrate import *


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

robot.home()

p10s_tipracks = [
    containers.load('tiprack-10ul', 'E1'),
    containers.load('tiprack-10ul', 'E2')
]
trash = containers.load('point', 'D1', 'holywastedplasticbatman')
well_plate = containers.load('96-deep-well', 'D2')
trans_plate = containers.load('point','B2')


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

def run_ot(coords):
	'''
	Pass the coordinates to the robot to pick the colony
	'''
	for i,[x,y] in enumerate(coords):
		print(i,":",x,y)
		p10s.move_to((trans_plate,[x,y,0]))
		input("click to move to next colony")

y_max = 120
x_max = 78

coords = [[0,0],[0,y_max],[x_max,y_max],[x_max,0]]
cols = [10,20,30,40,50,60,70]
y = range(10,120,10)
for col in cols:
    x = [col for num in range(12)]
    coords += zip(x,y)

run_ot(coords)











#
