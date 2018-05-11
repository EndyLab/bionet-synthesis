import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import json
import os
import numpy as np
import pandas as pd
import glob
from datetime import datetime

from config import *
from db_config import *
session,engine = connect_db()

data = []
mis = []
counter = 1
no_order = []
max = int(input("Max number to query: "))
for part in session.query(Part).order_by(Part.part_id):
    counter += 1
    if counter == max:
        break
    print(part.part_id)
    try:
        order_num = int(part.fragments[0].twist_orders[0].sub_name[-3:])
    except:
        order_num = 12
        no_order.append(part.part_id)

    frags = len([frag for frag in part.fragments])
    attempts = [well for well in part.wells if well.plates.plate_type == 'seq_plate' and well.misplaced != 'True']
    tups = []
    for well in attempts:
        build = int(well.plates.builds.build_name[-3:])
        tups.append([well.seq_outcome,build])
    tups = sorted(tups, key=lambda tup: tup[1])
    attempts = [out for out,num in tups]
    misplaced = [well.seq_outcome for well in part.wells if well.plates.plate_type == 'seq_plate' and well.misplaced == 'True']
    if len(misplaced) == 0:
        misplaced = ['N/A']
    mis.append(misplaced)
    if len(attempts) == 0:
        attempts += ['N/A','N/A']
    elif len(attempts) == 1:
        attempts += ['N/A']
    elif len(attempts) == 2:
        print('two attemps')
    else:
        print(len(attempts))
        input('More than two attempts')
    row = [part.part_id,part.status,len(part.seq),frags,attempts[0],attempts[1],part.cloning_enzyme,order_num]
    data.append(row)


data = np.array(data)
data_df = pd.DataFrame({
    'Part':data[:,0],
    "Status":data[:,1],
    'Length':data[:,2],
    'Frags':data[:,3],
    'Enzyme':data[:,6],
    'Order':data[:,7],
    'Attempt_1':data[:,4],
    'Attempt_2':data[:,5],
    'Misplaced':mis
})
simple_status = []
for attempts in zip(data_df['Attempt_1'].tolist(),data_df['Attempt_2'].tolist()):
    row = []
    for att in attempts:
        if "mutation" in att:
            row.append('cloning_mutation')
        elif "bad" in att:
            row.append('sequence_failure')
        elif att == 'cloning_error':
            row.append('cloning_failure')
        else:
            row.append(att)
    simple_status.append(row)
simple_status = np.array(simple_status)
data_df['Attempt_1_G'] = simple_status[:,0]
data_df['Attempt_2_G'] = simple_status[:,1]

data_df = data_df[['Part','Status','Length','Frags','Enzyme','Order','Attempt_1_G','Attempt_2_G','Attempt_1','Attempt_2','Misplaced']]
data_df['Length'] = pd.to_numeric(data_df['Length'])
data_df['Frags'] = pd.to_numeric(data_df['Frags'])

data_df.sort_values(by=['Length'])

data_b = data_df[data_df['Enzyme'] == 'BbsI']



# sankey = {}
#
# nodes = ['Total_ordered']
# links = []
# sankey['nodes'] = []
# sankey['links'] = []
#
# for order in order_raw.groupby('Order'):
#     order_num = int(order[0])
#     order_name = 'Order_' + str(order_num).zfill(3)
#
#     order_total = sum(order[1].Count.tolist())
#     links.append(['Total_ordered',order_name,order_total])
#
#     cur_nodes = order[1].Status.tolist()
#     cur_nodes = [node+'_'+str(order_num) for node in cur_nodes]
#     nodes += cur_nodes
#     nodes.append(order_name)
#
#     for i,row in order[1].iterrows():
#         links.append([order_name,row.Status+'_'+str(order_num),row.Count])
#
# node_dict = dict([[y,x] for x,y in enumerate(pd.Series(nodes).unique())])
#
# for name in node_dict.keys():
#     sankey['nodes'].append({'name' : name})
#
# for source,target,value in links:
#     sankey['links'].append({
#         "source":node_dict[source],
#         "target":node_dict[target],
#         "value":value
#     })
# with open("../docs/sankey/sankey.json","w+") as json_file:
#     json.dump(sankey,json_file,indent=2)



def add_branch(source,target,amount,nodes,links):
    if amount == 0:
        return nodes,links
    nodes += [source,target]
    links.append([source,target,amount])
    return nodes,links

def gen_sankey(nodes,links,sankey={'nodes':[],'links':[]}):
    node_dict = dict([[y,x] for x,y in enumerate(pd.Series(nodes).unique())])
    for name in node_dict.keys():
        sankey['nodes'].append({'name' : name})
    for source,target,value in links:
        sankey['links'].append({
            "source":node_dict[source],
            "target":node_dict[target],
            "value":value
        })
    return sankey

sankey = {}

nodes = ['Total_ordered']
links = []
sankey['nodes'] = []
sankey['links'] = []

# Orders
total = len(data_b)
orders = pd.DataFrame(data_b.Order.value_counts())
order_names = ['Order_'+str(order).zfill(3) for order in orders.index.tolist()]
order_counts = [int(count) for count in orders.Order.tolist()]

for order,count in zip(order_names,order_counts):
    nodes,links = add_branch(order,'Total_ordered',count,nodes,links)

# Abandoned
abandoned = len(data_b[data_b.Status == 'ordered'])
nodes,links = add_branch('Total_ordered','Abandoned',abandoned,nodes,links)

# Received
received = total-abandoned
nodes,links = add_branch('Total_ordered','Received',received,nodes,links)

# Attempted
not_attempted = len(data_b[data_b.Status == 'received'])
nodes,links = add_branch('Received','Not_attempted',not_attempted,nodes,links)
attempted = received-not_attempted
nodes,links = add_branch('Received','Attempted',attempted,nodes,links)

# Outcomes
data_att = data_b[data_b.Status != 'ordered']
outcomes = pd.DataFrame(data_att.Status.value_counts())
out = outcomes.index.tolist()
count = outcomes.Status.tolist()
for out,count in zip(out,count):
    print(out,count)
    nodes,links = add_branch('Attempted',out,count,nodes,links)


sankey = gen_sankey(nodes,links)

with open("../docs/sankey/sankey.json","w+") as json_file:
    json.dump(sankey,json_file,indent=2)
print("Check sankey")
outcomes



#
