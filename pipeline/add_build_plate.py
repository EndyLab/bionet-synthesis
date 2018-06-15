import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import pandas as pd
import sys
import glob
from config import *
import ot_functions as ot
from db_config import *
session,engine = connect_db()

def add_build_plate():
    build = ot.request_info("Please enter the build number you would like to enter into the database: ",type='int')

    build = 'build'+str(build).zfill(3)

    if len([build for build in session.query(Build).filter(Build.build_name == build)]) != 0:
        print("Build already exists")
        add_build_plate()
        sys.exit()

    build_path = '{}/builds/{}/{}_trans_map.csv'.format(BASE_PATH,build,build)
    if len(glob.glob(build_path)) == 0:
        print("No trans_map exists for this build")
        add_build_plate()
        sys.exit()
    data = pd.read_csv('{}/builds/{}/{}_trans_map.csv'.format(BASE_PATH,build,build))
    parts = [part for part in session.query(Part).filter(Part.part_id.in_(data.Gene.tolist()))]
    target_build = Build(parts,build_name=build)
    for i,(gene,well) in data.iterrows():
        print("working on: ", gene,well)
        gene_obj = session.query(Part).filter(Part.part_id == gene).one()
        target_build.add_item(gene_obj,well)

    ot.print_center("Build plate has been added")

    commit = int(ot.request_info("Commit changes (1-yes, 2-no): ",type='int'))
    if commit == 1:
        session.commit()

if __name__ == "__main__":
    add_build_plate()
