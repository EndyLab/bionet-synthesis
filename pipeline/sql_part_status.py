import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import numpy as np
import pandas as pd

from config import *
from db_config import *

status = []
not_attempted = []
# ignore = ['sequence_confirmed','cloning_mutation','building']
counter = 0
for part in session.query(Part).order_by(Part.id):
    if counter == 750:
        break
    counter += 1
    tries = len([well for well in part.wells if well.plates.plate_type == 'seq_plate' and well.misplaced != 'True'])
    if tries == 0 and part.status != 'ordered' and part.status != 'received':
        array = [tries,part.part_id,part.status,[[well.misplaced,well.plates.builds.build_name] for well in part.wells if well.plates.plate_type == 'seq_plate']]
        print(array)
        not_attempted.append(array)
    print(part.part_id,part.status,str(tries))
    frags = len([frag for frag in part.fragments])
    state = part.status+" - "+str(tries)+" - "+part.cloning_enzyme+" - "+str(frags)
    status.append(state)



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
for gene in not_attempted:
    print(gene)
# multi = pd.Series(multi)
# print(multi.value_counts())
