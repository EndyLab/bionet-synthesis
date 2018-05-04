import numpy as np
import pandas as pd
import glob

from config import *
from db_config import *

def add_order(file,no_frag=[],prev_subs=[]):
    order = pd.read_csv(file)
    sub_name = file.split('/')[-1][:-4]
    if sub_name in prev_subs:
        return no_frag
    new_order = Twist_order(sub_name=sub_name)
    session.add(new_order)
    for i,row in order.iterrows():
        print(i,row['gene name'])
        seq = row['FASTA_seq'].upper()
        frag = session.query(Fragment).filter(Fragment.seq == seq).first()
        if type(frag) == type(None):
            print("No fragment found")
            no_frag.append(row['gene name'])
        else:
            new_order.fragments.append(frag)
    for frag in new_order.fragments:
        print(frag.fragment_name)
    return no_frag



if __name__ == "__main__":
    no_frag = []
    prev_subs = [sub.sub_name for sub in session.query(Twist_order).order_by(Twist_order.id)]
    for file in sorted(glob.glob('{}/submissions/submission*.csv'.format(BASE_PATH))):
        print(file)
        no_frag = add_order(order,no_frag=no_frag,prev_subs=prev_subs)
    print("Missing fragments:")
    print(no_frag)
    commit = int(input("Commit changes (1-yes, 2-no): "))
    if commit == 1:
        session.commit()











#
