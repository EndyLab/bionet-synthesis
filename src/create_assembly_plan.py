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
import argparse

from config import *
import ot_functions as ot

from db_config import *

def create_build_plans(session,engine,num_parts,frag_nums,enzyme='BbsI',max_rxns=96):
    # Take in the 'run' argument from the command line
    parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
    parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
    args = parser.parse_args()

    ot.print_center("============ Creating Plan ============")
    num_builds = math.ceil(num_parts / 96)
    print(num_builds)
    used_plates = []

    if num_parts < max_rxns:
        print('Generating a partial build')
        max_rxns = num_parts

    previous_plans = [build.build_name for build in session.query(Build).order_by(Build.build_name).filter(Build.status == 'planning')]
    if len(previous_plans) > 0:
        print("Are you sure you don't want to use these plans instead:\n",previous_plans)
        ans = ot.request_info("Use one of these plans instead (y/n):\n{}".format(previous_plans))
        if ans == 'y':
            print('Enter the plan into the build script')
            sys.exit()
    print('Will continue to exclude parts in previous plans')

    past_builds = pd.read_sql_query("SELECT builds.build_name FROM builds", con=engine).build_name.tolist()
    if len(past_builds) == 0:
        last_build = 0
    else:
        last_build = [int(build[-3:]) for build in past_builds][-1]
    print('last_build',last_build)
    new_builds = ["build"+str(last_build+num).zfill(3) for num in range(1,num_builds+1)]
    print(new_builds)

    acceptable_status = ['received'] # List with all of the status that you want to pull from
    rework = ['cloning_error','cloning_failure','trans_failure']

    def list_to_string(ls):
        string = ''
        for l in ls:
            string += "'{}',".format(l)
        return string[:-1]

    def query_for_parts(status,enzyme,engine):
        query_parts = "SELECT parts.part_id,parts.status,fragments.fragment_name,plates.plate_id,wells.address,wells.volume,plates.id FROM parts\
                INNER JOIN part_frag ON parts.id = part_frag.part_id\
                INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
                INNER JOIN wells ON fragments.id = wells.fragment_id\
                INNER JOIN plates on wells.plate_id = plates.id\
                WHERE parts.status IN ({})\
                    AND parts.cloning_enzyme = '{}'".format(list_to_string(status),enzyme)

        return pd.read_sql_query(query_parts, con=engine)

    for build,frag_num in zip(new_builds,frag_nums):
        ## Pulls in parts that have the desirable status and then sorts them by the plates their fragments are in,
        ## sorts them by their plate numbers and returns the earliest set of parts
        ot.print_center("...Finding genes for {}...".format(build))

        to_build = query_for_parts(acceptable_status,enzyme,engine)

        to_build.id = to_build.plate_id.apply(lambda x: 0 if x in used_plates else 1)

        grouped = to_build.groupby('part_id').filter(lambda x: len(x) == x['id'].sum())
        # print(grouped)
        if frag_num != 0:
            grouped = grouped.groupby('part_id').filter(lambda x: len(x) == frag_num)

        sort_group = grouped.sort_values('plate_id')
        # print(sort_group)
        sub = sort_group.part_id.unique().tolist()

        if len(sub) < max_rxns:
            print('Not enough remaining parts to fill a build')
            rework_ans = ot.request_info("Include previously failed genes? (y/n): ",type='string')
            if rework_ans == 'y':
                rework_df = query_for_parts(rework,enzyme,engine).sort_values('plate_id')
                sort_group = pd.concat([sort_group,rework_df])
                sub = sort_group.part_id.unique().tolist()[:max_rxns]
        else:
            print('truncating parts')
            sub = sub[:max_rxns]
            print(len(sub))

        to_build = sort_group[sort_group.part_id.isin(sub)]
        print(len(to_build))
        # print('To Build:\n',to_build)
        num_reactions = len(to_build)

        print("{} includes {} parts".format(build,len(to_build.part_id.unique().tolist())))

        print("Fragments sourced from each plate: ")
        print(to_build.groupby('part_id').agg(len).status.value_counts())

        used_plates += to_build.plate_id.unique().tolist()
        print(to_build.plate_id.unique().tolist())

        master_mix = ot.make_gg_rxns(num_reactions,10)
        print("Below is the required master mix for {}\n{}".format(build,master_mix))

        build_parts = to_build.part_id.unique().tolist()

        if args.run:

            ot.print_center('...Updating the database...')

            build_objs = [part for part in session.query(Part).filter(Part.part_id.in_(build_parts))]

            target_build = Build(build_objs,build_name=build)
            target_build.status = 'planning'
            session.add(target_build)

            for part in build_objs:
                part.change_status('planning')

            session.commit()
        else:
            print('Did not write plan. Add -r to the end of the function to write them')






if __name__ == "__main__":
    session,engine = connect_db()

    df = ot.query_everything(engine)
    df_rec = df[df.status == 'received']
    print('Below is the breakdown of remaining constructs by fragment number:')
    print(df_rec.Fragments.value_counts())


    print('1 build = 96 parts\n2 builds = 192 parts\n3 builds = 288 parts')
    num_parts = int(ot.request_info("Enter the desired number of parts to build: ",type='int'))
    num_builds = math.ceil(num_parts / 96)
    frag_nums = ot.request_info("Enter the number of fragments desired for the {} builds \n(i.e. '1 2 0' where 0 represents any number): ".format(num_builds), type='list',length=num_builds)

    enzyme_num = ot.request_info("Which enzyme (1 - BbsI or 2 - BtgZI): ",type='int')
    if int(enzyme_num) == 1:
        enzyme = 'BbsI'
    else:
        enzyme = 'BtgZI'
    print("Using enzyme: ",enzyme)

    create_build_plans(session,engine,num_parts,frag_nums,enzyme=enzyme)









#
