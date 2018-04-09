import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import numpy as np
import pandas as pd

from config import *
from db_config import *

status = []
for part in session.query(Part):
    print(part.part_id,part.status)
    status.append(part.status)

status = pd.Series(status)
print(status.value_counts())
