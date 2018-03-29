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
from db_gen import *
from sql_frag_loc import *
from sql_resuspension import *
from sql_build import *
from sql_plating import *

session = build_db()
session = frag_assign(session)
session = resuspension(session)
session = resuspension(session)
session = run_build(session)
session = plate(session)

assemblies = []
print("Choose which plate you would like to assess:")
for index,assembly in enumerate(session.query(Plate).filter(Plate.assessed == 'assessed').order_by(Plate.id)):
    print("{}. {}-{}".format(index,assembly.builds.id,assembly.plate_name))
    assemblies.append(assembly)

plate_num = int(input("Enter plate here: "))
# target_plate = assemblies[plate_num]

print("Enter well addresses that didn't work. e.g. ['A1','H6']")
failed = input("List: ")
print(failed)

## =============================================
## QUERY DATABASE FOR REPLACEMENTS
## =============================================

# for well in session.query(Well).join(Part,Well.parts










#
