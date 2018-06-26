import sys
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
import ot_functions as ot

session,engine = connect_db()

builds = []
print('Builds awaiting a manual check:')
for i,build in enumerate(session.query(Build).filter(Build.status == 'manual_check_pending').order_by(Build.build_name)):
    print('{}: {}'.format(i,build.build_name))
    builds.append(build)

ans = ot.request_info('Which build would you like to check: ',type='int',select_from=[x for x in range(len(builds))])
tar_build = builds[ans]

data = pd.read_csv(glob.glob('{}/builds/{}/{}_manual_check.csv'.format(BASE_PATH,tar_build.build_name,tar_build.build_name))[0])
outcome_dict = {
    'k':'keep',
    'p':'perfect',
    'm':'mutation',
    'f':'cloning_failure'
}
data['manual'] = data.manual.apply(lambda x: outcome_dict[x])

ot.print_center('...Updating the database...')

for i,row in data.iterrows():
    tar_well = session.query(Well).join(Plate,Well.plates).join(Build,Plate.builds)\
        .join(Part,Well.parts).filter(Build.build_name == tar_build.build_name)\
        .filter(Part.part_id == row['part_id']).filter(Plate.plate_type == 'seq_plate')\
        .filter(Well.address == row['address']).one()
    tar_well.seq_outcome = row['manual']

session.commit()

query = "SELECT wells.seq_outcome,wells.misplaced FROM wells\
            INNER JOIN plates ON wells.plate_id = plates.id\
            INNER JOIN builds ON plates.build_id = builds.id\
            WHERE builds.build_name = '{}'\
                AND plates.plate_type = 'seq_plate'\
                AND wells.misplaced IS NULL".format(tar_build.build_name)

results = pd.read_sql_query(query,con=engine)

print(results.seq_outcome.value_counts())





#
