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


parts = [part for part in session.query(Part).limit(50)]

new_plate = Plate(parts,'seq_plate','test')

# for part in parts:
#     new_plate.add_item(part)
    # input()

new_plate.add_item(parts[1],address='A1')
new_plate.add_item(parts[2],address='H10')
new_plate.add_item(parts[3],address='D7')
new_plate.add_item(parts[4])
new_plate.add_item(parts[5])
new_plate.add_item(parts[6])













#
