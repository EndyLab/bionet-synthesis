import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import numpy as np
import pandas as pd
from datetime import datetime

from config import *
import ot_functions as ot
from db_config import *
session,engine = connect_db()
conn = engine.connect()

from moclopy import fixer
from synbiolib import codon


query_seqs = "SELECT parts.part_id,parts.part_name,parts.part_type,parts.original_seq,parts.seq,parts.organism FROM parts\
                WHERE parts.seq IS NULL"

df = pd.read_sql_query(query_seqs,con=engine)

def optimize_sequences(row):
    if row['part_type'] != 'cds':
        return row['original_seq']
    if row['organism'] == None:
        table = codon.load_codon_table(species='ecoli')
    else:
        table = codon.load_codon_table(species=row['organism'])
    protein_seq = fixer.translate(row['original_seq'])
    optimized = codon.optimize_protein(table,protein_seq)

    fixed = fixer.fix_sequence(row['part_id'],optimized)

    return fixed

ot.print_center("...Optimizing sequences...")
print("Starting at: {}".format(datetime.now()))
df = df.loc[0:99,:]
df['seq'] = df.apply(optimize_sequences, axis=1)
print("Finished at: {}".format(datetime.now()))

seq_dict = dict(zip(df.part_id.tolist(),df.seq.tolist()))

for part in session.query(Part).filter(Part.part_id.in_(df.part_id.tolist())):
    part.seq = seq_dict[part.part_id]

session.commit()
