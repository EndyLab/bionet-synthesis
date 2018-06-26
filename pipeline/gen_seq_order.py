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

import ot_functions as ot
from config import *
from db_config import *
session,engine = connect_db()


def gen_seq_order():
    # Determine which build/plate to assess
    assemblies = []
    print("Choose which builds you would like to send for sequencing:")
    for index,assembly in enumerate(session.query(Build).filter(Build.status == 'sequencing').order_by(Build.build_name)):
        print("{}. {}".format(index,assembly.build_name))
        assemblies.append(assembly)
    ans = ot.request_info('Enter numbers separated by a space (i.e. 10 13): ',type='list')
    builds = [assemblies[num] for num in ans]
    print('Will add the following builds to the order:')
    print([build.build_name for build in builds])

    ans = ot.request_info('Are these correct? ',type='string',select_from=['y','n'])
    if ans == 'n':
        ans = ot.request_info('Enter numbers separated by a space (i.e. 10 13): ',type='list')
        builds = [assemblies[num] for num in ans]
        print('Will add the following builds to the order:')
        print([build.build_name for build in builds])

    for build in builds:
        query = "SELECT parts.part_id,wells.address,parts.seq,builds.build_name from builds\
                    INNER JOIN plates ON builds.id = plates.build_id\
                    INNER JOIN wells ON plates.id = wells.plate_id\
                    INNER JOIN parts ON wells.part_id = parts.id\
                    WHERE builds.build_name = '{}'\
                        AND plates.plate_type = 'seq_plate'".format(build.build_name)

        data = pd.read_sql_query(query,con=engine)

        data['Template Type'] = 'Plasmid'
        data['Template Conc (ng/ul)'] = 50
        data['Template Size (bp)'] = data.seq.apply(len)
        data = data.rename(columns={'part_id':'Template Name','address':'Well'})
        data = data[['Well','Template Name','Template Type','Template Conc (ng/ul)','Template Size (bp)']]

        path = '{}/builds/{}'.format(BASE_PATH,build.build_name)
        ot.make_directory(path)
        data.to_csv('{}/{}_seq_order.csv'.format(path,build.build_name),index=False)
        input("Build complete")






if __name__ == "__main__":
    gen_seq_order()










#
