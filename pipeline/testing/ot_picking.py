from opentrons import robot, containers, instruments
import argparse
import sys

import os
import re
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime
import getch


port = os.environ["ROBOT_DEV"]
print("Connecting robot to port {}".format(port))
robot.connect(port)

robot.home()

p10s_tipracks = [
    containers.load('tiprack-10ul', 'E1'),
    containers.load('tiprack-10ul', 'E2')
]
trash = containers.load('point', 'D1', 'holywastedplasticbatman')
agar_plate = containers.load('96-deep-well', 'D2')
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

y_max = 118
x_max = 78
for num in range(2):
    print(num)
    p10s.move_to(trans_plate)
    p10s.move_to((trans_plate,[0,y_max,0]))
    p10s.move_to((trans_plate,[x_max,y_max,0]))
    p10s.move_to((trans_plate,[x_max,0,0]))

print('end')
