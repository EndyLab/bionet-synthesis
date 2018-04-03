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

from config import *
from db_config import *

for instance in session.query(Plate):#.filter(Plate.plate_type == 'syn_plate'):
    print(instance.plate_name)
    if '38' in instance.plate_name:
        for well in instance.wells:
            print(well.address,well.syn_yield,well.volume)
