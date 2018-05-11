import numpy as np
import pandas as pd
import glob

from config import *
import db_config
session,engine = db_config.connect_db()

def add_order(file,no_frag=[],prev_subs=[]):
    order = pd.read_csv(file)
    sub_name = file.split('/')[-1][:-4]
    prev_frags = []
    if sub_name not in prev_subs:
        print('making a new order')
        current_order = Twist_order(sub_name=sub_name)
        session.add(current_order)
    else:
        print('order already exists')
        current_order = session.query(Twist_order).filter(Twist_order.sub_name == sub_name).one()
        for frag in current_order.fragments:
            prev_frags.append(frag.seq)
    added = []
    for i,row in order.iterrows():
        print(i,row['gene name'])
        seq = row['FASTA_seq'].upper()
        if seq in prev_frags:
            print('Previously added')
            continue
        frag = session.query(Fragment).filter(Fragment.seq == seq).first()
        if type(frag) == type(None):
            print("No fragment found")
            no_frag.append(row['gene name'])
        else:
            added.append(frag.fragment_name)
            current_order.fragments.append(frag)
    for frag in current_order.fragments:
        print(frag.fragment_name)
    print(added)
    return no_frag



if __name__ == "__main__":
    no_frag = []
    prev_subs = [sub.sub_name for sub in session.query(Twist_order).order_by(Twist_order.id)]
    for file in sorted(glob.glob('{}/submissions/submission005.csv'.format(BASE_PATH))):
        print(file)
        no_frag = add_order(file,no_frag=no_frag,prev_subs=prev_subs)
    print("Missing fragments:")
    print(no_frag)
    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()











#
