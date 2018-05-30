import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime
import getch
import shutil
import sys

from config import *
import ot_functions as ot

from db_config import *
session,engine = connect_db()

def create_build_plans(num_parts,enzyme='BbsI',max_rxns=96):
    ot.print_center("============ Creating Plan ============")
    num_builds = math.ceil(num_parts / 96)
    past_builds = pd.read_sql_query("SELECT builds.build_name FROM builds", con=engine).build_name.tolist()
    last_build = max([int(build[-3:]) for build in past_builds])
    new_builds = ["build"+str(last_build+num).zfill(3) for num in range(1,num_builds+1)]
    print(new_builds)

    acceptable_status = ['received'] # List with all of the status that you want to pull from
    rework = ['cloning_error','cloning_failure','trans_failure']

    used_plates = []

    for build in new_builds:
        ## Pulls in parts that have the desirable status and then sorts them by the plates their fragments are in,
        ## sorts them by their plate numbers and returns the earliest set of parts
        ot.print_center("...Finding genes for {}...".format(build))

        to_build = [part for part in session.query(Part)\
                    .join(Fragment,Part.fragments)\
                    .join(Well,Fragment.wells)\
                    .join(Plate,Well.plates)\
                    .filter(Part.status.in_(acceptable_status))\
                    .filter(Part.cloning_enzyme == enzyme)\
                    .filter(Plate.plate_id.notin_(used_plates))\
                    .order_by(Plate.plate_id)]

        if len(to_build) < 96:
            print('Not enough remaining parts to fill a build')
            rework_ans = ot.request_info("Include previously failed genes? (y/n): ")
            if rework_ans == 'y':
                to_build += [part for part in session.query(Part)\
                            .join(Fragment,Part.fragments)\
                            .join(Well,Fragment.wells)\
                            .join(Plate,Well.plates)\
                            .filter(Part.status.in_(rework))\
                            .filter(Part.cloning_enzyme == enzyme)\
                            .filter(Plate.plate_id.notin_(used_plates))\
                            .order_by(Plate.id)]
        to_build = to_build[:max_rxns]
        print("{} includes {} parts".format(build,len(to_build)))

        ot.print_center("...Finding the required plates...")
        unique_plates = []
        for part in to_build:
            for frag in part.fragments:
                unique_plates.append(frag.wells[0].plates.plate_id)
        unique_plates = pd.Series(unique_plates).unique().tolist()
        print(unique_plates)
        used_plates += unique_plates



if __name__ == "__main__":
    print('1 build = 96 parts\n2 builds = 192 parts\n3 builds = 288 parts')
    num_builds = ot.request_info("Enter the desired number of builds: ")
    max_rxns = int(num_builds) * 96

    enzyme_num = ot.request_info("Which enzyme (1 - BbsI or 2 - BtgZI): ")
    if int(enzyme_num) == 1:
        enzyme = 'BbsI'
    else:
        enzyme = 'BtgZI'
    print("Using enzyme: ",enzyme)

    create_build_plans(max_rxns,enzyme=enzyme)
