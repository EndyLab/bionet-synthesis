import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import numpy as np
import pandas as pd

from config import *
from db_config import *

status = []
multi = []
# ignore = ['sequence_confirmed','cloning_mutation','building']
for part in session.query(Part).order_by(Part.id):
    print(part.part_id,part.status)
    # if part.status in ignore:
    #     continue
    status.append(part.status)
    tries = [(well.plates.builds.build_name,well.seq_outcome) for well in part.wells if well.plates.plate_type == 'seq_plate']

    # if len(tries) > 1:
    # print(part.part_id,tries)
        # multi.append(tries)


# for build in session.query(Build).order_by(Build.id):
#     for plate in build.plates:
#         print(build.build_name)
#         if plate.plate_type != 'seq_plate':
#             continue
#         status = [well.seq_outcome for well in plate.wells]
#         status = pd.Series(status)
#         print(status.value_counts())


status = pd.Series(status)
print(status.value_counts())
# multi = pd.Series(multi)
# print(multi.value_counts())
