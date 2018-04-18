from opentrons import robot, containers, instruments
import argparse
import sys

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import re
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime
import getch

from config import *
from db_config import *

port = os.environ["ROBOT_DEV"]
print("Connecting robot to port {}".format(port))
robot.connect(port)

robot.home()

p10s_tipracks = [
    containers.load('tiprack-10ul', locations["tiprack-10s1"]),
    containers.load('tiprack-10ul', locations["tiprack-10s2"])
]
trash = containers.load('point', locations["trash"], 'holywastedplasticbatman')
agar_plate = containers.load('96-deep-well', 'D2')

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
for num in range(3):
    p10s.move_to(agar_plate)
    p10s.move_to((agar_plate,[0,5,0])








#
