import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime

import ot_functions as ot
from config import *
from db_config import *
session,engine = connect_db()

print(datetime.now(),'Began run')

query_outcomes = "SELECT parts.part_id,parts.status,wells.seq_outcome,wells.plate_type,builds.build_name,wells.misplaced FROM parts \
        INNER JOIN wells ON parts.id = wells.part_id\
        INNER JOIN plates ON wells.plate_id = plates.id\
        INNER JOIN builds ON plates.build_id = builds.id"

query_parts = "SELECT * FROM parts"

author_dict = []
for file in sorted(glob.glob('../data/*/*.json')):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    author_dict.append([data['gene_id'],data['author']['name']])
author_dict = dict(author_dict)

def multiple(x):
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
#     elif x == 'cloning_error':
#         return 'cloning_failure'
    else:
        return x

def find_author(x):
    return author_dict[x]

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

df_parts = pd.read_sql_query(query_parts, con=engine)
print('finished part query')

df_parts = df_parts[df_parts.part_id != 'BBF10K_000745']

df_parts['Outcomes'] = df_parts.part_id.apply(find_outcome)
df_parts['Builds'] = df_parts.part_id.apply(find_build)
print('finished outcome and builds')
df_parts['Attempt_1_Outcome'] = df_parts.Outcomes.apply(lambda x: x[0])
df_parts['Attempt_1_Outcome_G'] = df_parts.Attempt_1_Outcome.apply(simplify_outcome)
df_parts['Attempt_1_Build'] = df_parts.Builds.apply(lambda x: x[0])
df_parts['Attempt_2_Outcome'] = df_parts.Outcomes.apply(lambda x: x[1])
df_parts['Attempt_2_Outcome_G'] = df_parts.Attempt_2_Outcome.apply(simplify_outcome)
df_parts['Attempt_2_Build'] = df_parts.Builds.apply(lambda x: x[1])
df_parts['Author'] = df_parts.part_id.apply(find_author)
print(datetime.now(),'Finished building dataframe')

df_authors = df_parts[['part_name','part_id','status','Attempt_1_Outcome_G','Attempt_2_Outcome_G','Author']]

author_groups = df_authors.groupby('Author')
for author,parts in author_groups:
    print(author)
    print(parts)
    path = '{}/authors/{}'.format(BASE_PATH,author)
    ot.make_directory(path)
    parts.to_csv('{}/{}_status_{}'.format(path,author,str(datetime.now()).split(" ")[0]))
