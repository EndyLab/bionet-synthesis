import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import json
from shutil import copyfile
import os
import glob
import db_config
from config import *


def author_status():
    session,engine = db_config.connect_db()

    print(datetime.now(),'Began run')

    # Creates a dataframe with all of the parts that have been attempted and what their outcomes were
    query_outcomes = "SELECT parts.part_id,parts.status,wells.seq_outcome,wells.plate_type,builds.build_name,wells.misplaced FROM parts \
            INNER JOIN wells ON parts.id = wells.part_id\
            INNER JOIN plates ON wells.plate_id = plates.id\
            INNER JOIN builds ON plates.build_id = builds.id"

    # Creates a dataframe of all parts and their fragments and which orders they were submitted in
    query_frag = "SELECT parts.part_id,fragments.fragment_name,twist_orders.sub_name FROM parts\
            INNER JOIN part_frag ON parts.id = part_frag.part_id\
            INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
            INNER JOIN frag_order ON fragments.id = frag_order.frag_id\
            INNER JOIN twist_orders ON twist_orders.id = frag_order.twist_order_id"

    # Builds a table of all the current parts in the database
    query_parts = "SELECT * FROM parts"

    # Creates a dictionary to link the part_ids to the # of fragments and the order #
    df_frag = pd.read_sql_query(query_frag, con=engine)
    frags = df_frag.groupby('part_id')['fragment_name'].agg(len)
    frags.name = 'Count'
    frags = pd.DataFrame(frags).reset_index()
    frags_dict = dict(zip(frags.part_id.tolist(),frags.Count.tolist()))
    subs_dict = dict(zip(df_frag.part_id.tolist(),df_frag.sub_name.tolist()))
    print(datetime.now(),'Finished frags')

    # Creates a dictionary to link all of the parts with their author names
    author_dict = []
    for file in sorted(glob.glob('../data/*/*.json')):
        with open(file,"r") as json_file:
            data = json.load(json_file)
        author_dict.append([data['gene_id'],data['author']['name']])
    author_dict = dict(author_dict)

    def multiple(x):
        '''Adds an N/A for parts that have less than 2 attempts'''
        if len(x) == 1:
            x.append('N/A')
        return x

    def find_outcome(x):
        if x in df_out_dict.keys():
            return df_out_dict[x]
        else:
            return ['N/A','N/A']

    def find_build(x):
        if x in df_build_dict.keys():
            return df_build_dict[x]
        else:
            return ['N/A','N/A']

    def simplify_outcome(x):
        if "mutation" in x:
            return 'cloning_mutation'
        elif "bad" in x:
            return 'sequence_failure'
        else:
            return x

    # Takes in all of the build outcome information
    df_res = pd.read_sql_query(query_outcomes, con=engine)
    df_res = df_res[df_res.plate_type == 'seq_plate']

    df_out = df_res.groupby('part_id')['seq_outcome'].apply(list)
    df_out.name = 'Outcomes'
    df_out = pd.DataFrame(df_out).reset_index()
    df_out.Outcomes = df_out.Outcomes.apply(multiple)
    df_out_dict = dict(zip(df_out.part_id.tolist(),df_out.Outcomes.tolist()))

    df_build = df_res.groupby('part_id')['build_name'].apply(list)
    df_build.name = 'Builds'
    df_build = pd.DataFrame(df_build).reset_index()
    df_build.Builds = df_build.Builds.apply(multiple)
    df_build_dict = dict(zip(df_build.part_id.tolist(),df_build.Builds.tolist()))
    print(datetime.now(),'Finished outcomes')

    # Creates a comprehensive tabel that incorporates all of the gathered information
    df_parts = pd.read_sql_query(query_parts, con=engine)
    print('finished part query')
    df_parts['Fragments'] = df_parts.part_id.apply(lambda x: frags_dict[x])
    df_parts['Submission'] = df_parts.part_id.apply(lambda x: subs_dict[x])
    df_parts['Order_number'] = df_parts.Submission.apply(lambda name: int(name[-3:]))
    df_parts['Outcomes'] = df_parts.part_id.apply(find_outcome)
    df_parts['Builds'] = df_parts.part_id.apply(find_build)
    print('finished outcome and builds')
    df_parts['Attempt_1_Outcome'] = df_parts.Outcomes.apply(lambda x: x[0])
    df_parts['Attempt_1_Outcome_G'] = df_parts.Attempt_1_Outcome.apply(simplify_outcome)
    df_parts['Attempt_1_Build'] = df_parts.Builds.apply(lambda x: x[0])
    df_parts['Attempt_2_Outcome'] = df_parts.Outcomes.apply(lambda x: x[1])
    df_parts['Attempt_2_Outcome_G'] = df_parts.Attempt_2_Outcome.apply(simplify_outcome)
    df_parts['Attempt_2_Build'] = df_parts.Builds.apply(lambda x: x[1])
    df_parts['Length'] = df_parts.seq.apply(len)
    df_parts['Author'] = df_parts.part_id.apply(lambda x: author_dict[x])

    print(datetime.now(),'Finished building dataframe')

    print(df_parts.Author.value_counts())

    author_name = input('Enter author name: ')
    author_df = df_parts[df_parts.Author == author_name]
    limited_df = author_df[['part_name','part_id','status','Attempt_1_Outcome_G','Attempt_2_Outcome_G']]
    print(limited_df)

    current = str(datetime.now()).split(" ")[0]

    # Generates a directory and csv file of all of the cloning information
    file_name = "{}_status_{}.csv".format(author_name,current)
    path = '{}/authors/{}'.format(BASE_PATH,author_name)
    if os.path.exists(path):
        print("Directory for {} already exists".format(author_name))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)
        print("Making directory for {}".format(author_name))
    file_path = path + "/" + file_name
    limited_df.to_csv(file_path,index=False)

if __name__ == "__main__":
    author_status()











#
