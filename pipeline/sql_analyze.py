import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import json
from shutil import copyfile
import os

import matplotlib
font = {'size'   : 15}
matplotlib.rc('font', **font)

from db_config import *
session,engine = connect_db()

print(datetime.now(),'Began run')

query_outcomes = "SELECT parts.part_id,parts.status,wells.seq_outcome,wells.plate_type,builds.build_name,wells.misplaced FROM parts \
        INNER JOIN wells ON parts.id = wells.part_id\
        INNER JOIN plates ON wells.plate_id = plates.id\
        INNER JOIN builds ON plates.build_id = builds.id"

query_frag = "SELECT parts.part_id,fragments.fragment_name,twist_orders.sub_name FROM parts\
        INNER JOIN part_frag ON parts.id = part_frag.part_id\
        INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
        INNER JOIN frag_order ON fragments.id = frag_order.frag_id\
        INNER JOIN twist_orders ON twist_orders.id = frag_order.twist_order_id"

query_parts = "SELECT * FROM parts"

df_frag = pd.read_sql_query(query_frag, con=engine)
frags = df_frag.groupby('part_id')['fragment_name'].agg(len)
frags.name = 'Count'
frags = pd.DataFrame(frags).reset_index()
frags_dict = dict(zip(frags.part_id.tolist(),frags.Count.tolist()))
subs_dict = dict(zip(df_frag.part_id.tolist(),df_frag.sub_name.tolist()))
print(datetime.now(),'Finished frags')

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
    else:
        return x

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
df_parts = df_parts[df_parts.part_id != 'BBF10K_000745']
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
print(datetime.now(),'Finished building dataframe')
df_parts

## BBREAKDOWN OF OUTCOMES BY NUMBER OF FRAGMENTS

data_b = df_parts[df_parts.cloning_enzyme == 'BbsI']

# data_fail = data_b[data_b.status != 'sequence_confirmed']
data_att = data_b[data_b.status != 'ordered']

frag_norm = (data_att.groupby('Fragments')['status'].value_counts() / data_att.groupby('Fragments')['status'].agg(len))
frag_norm.name = 'Percent'
norm_frag_df = pd.DataFrame(frag_norm).reset_index()
norm_frag_df.Percent = norm_frag_df.Percent * 100

dims = (16, 6)

fig1, ax1 = plt.subplots(figsize=dims)
g = sns.barplot(ax=ax1,data=norm_frag_df, x='Fragments', y='Percent', hue='status',palette="Spectral")
g.set_yscale('log')
plt.title("Normalized Stat Percentage vs. Number of Fragments")
fig1.show()
input("move to next image")
plt.close()

total_frag_raw = (data_att.groupby('Fragments')['status'].value_counts())
total_frag_raw.name = 'Count'
total_raw = pd.DataFrame(total_frag_raw).reset_index()

fig2, ax2 = plt.subplots(figsize=dims)
t = sns.barplot(ax=ax2,data=total_raw, x='Fragments', y='Count', hue='status',palette='Spectral')
t.set_yscale('log')
plt.title("Raw Status Count vs. Number of Fragments")
fig2.show()
input("move to next image")
plt.close()

print(datetime.now(),'Finished outcome vs. fragments\n')

data_att.groupby(data_att.status).Fragments.describe()

## CLONING FAILURE VERSUS LENGTH

data_fail = data_att[data_att.status != 'sequence_confirmed']
data_fail = data_fail[data_fail.status != 'abandoned']
data_fail = data_fail[data_fail.status != 'received']


fig3, ax3 = plt.subplots()

data_fail.groupby(data_fail.status).Length.plot.hist(ax=ax3,alpha=0.5,legend=True,figsize=(16,6),bins=10)
data_fail.groupby(data_fail.status).Length.describe()

plt.xlabel('Seqence Length (bp)')
plt.ylabel('Counts')
plt.title('Cloning Failure vs. Sequence Length')

fig3.savefig('../docs/Overall/raw_length.png')

## SYNTHESIS FAILURES PLOTTED AGAINST SEQUENCE LENGTH

data_order = data_b[data_b.status == 'abandoned']
fig4, ax4 = plt.subplots()

ax4.set_xticks([1000,2000,3000,4000,5000])

data_order.groupby(data_order.status).Length.plot.hist(ax=ax4,alpha=0.5,legend=True,figsize=(16,6),bins=10)

plt.xlabel('Seqence Length (bp)')
plt.ylabel('Counts')
plt.title('Synthesis Abandonment Versus Sequence Length')
fig4.show()
input("move to next image")
plt.close()
fig4.savefig('../docs/Overall/syn_fail.png')

## SECOND ATTEMPT BREAKDOWN

# Elimination of irrelevant information for the plots
data_att = data_b[data_b.status != 'ordered']
data_att = data_att[data_att.status != 'received']
data_att = data_att[data_att.Attempt_1_Outcome_G != 'sequence_confirmed']
data_attempt = data_att[data_att.Attempt_2_Outcome_G != 'N/A']

# Normalizing the data by dividing the counts of each 2nd outcome by the total counts of each group of 1st outcomes
attempts_norm = (data_attempt.groupby('Attempt_1_Outcome_G')['Attempt_2_Outcome_G'].value_counts() / data_attempt.groupby('Attempt_1_Outcome_G')['Attempt_2_Outcome_G'].agg(len))
attempts_norm.name = 'Percent'
attempts_norm = pd.DataFrame(attempts_norm).reset_index()
attempts_norm.Percent = attempts_norm.Percent * 100

# Plot normalized data
dims = (16, 6)
fig1, ax = plt.subplots(figsize=dims)
p = sns.barplot(data=attempts_norm, x='Attempt_1_Outcome_G', y='Percent', hue='Attempt_2_Outcome_G',palette="Spectral").set_title("Normalized 2nd Attempt Outcomes")
plt.xlabel("1st Attempt Outcome")
fig1.show()
input("move to next image")
plt.close()
fig1.savefig('../docs/Overall/norm_attempt.png')

# Caluclate the raw counts
attempts_raw = (data_attempt.groupby('Attempt_1_Outcome_G')['Attempt_2_Outcome_G'].value_counts())
attempts_raw.name = 'Count'
attempts_raw = pd.DataFrame(attempts_raw).reset_index()

# Plot the raw counts
fig2, ax = plt.subplots(figsize=dims)
g = sns.barplot(data=attempts_raw, x='Attempt_1_Outcome_G', y='Count', hue='Attempt_2_Outcome_G',palette="Spectral").set_title("2nd Attempt Outcomes")
plt.xlabel("1st Attempt Outcome")
fig2.show()
input("move to next image")
plt.close()
fig2.savefig('../docs/Overall/raw_attempt.png')

## BREAKING DOWN CLONING OUTCOMES BY ORDER NUMBER

# Normalizing the data by dividing the counts of each 2nd outcome by the total counts of each group of 1st outcomes
order_norm = (data_b.groupby('Order_number')['status'].value_counts() / data_b.groupby('Order_number')['status'].agg(len))
order_norm.name = 'Percent'
order_norm = pd.DataFrame(order_norm).reset_index()
order_norm.Percent = order_norm.Percent * 100

# Plot normalized data
dims = (16, 6)
fig1, ax = plt.subplots(figsize=dims)
p = sns.barplot(data=order_norm, x='Order_number', y='Percent', hue='status',palette="Spectral").set_title("Normalized Outcomes by Order")
plt.xlabel("Order Number")
fig1.show()
input("move to next image")
plt.close()
fig1.savefig('../docs/Overall/norm_order.png')

# Caluclate the raw counts
order_raw = (data_b.groupby('Order_number')['status'].value_counts())
order_raw.name = 'Count'
order_raw = pd.DataFrame(order_raw).reset_index()

# Plot the raw counts
fig2, ax = plt.subplots(figsize=dims)
g = sns.barplot(data=order_raw, x='Order_number', y='Count', hue='status',palette="Spectral")

g.set_yscale('log')
plt.title("Outcome Counts by Order")
plt.xlabel("Order Number")
fig2.show()
input("move to next image")
plt.close()
fig2.savefig('../docs/Overall/raw_order.png')

## BUILD ATTEMPT BREAKDOWN

df_int = df_res
df_int['Outcome'] = df_int.seq_outcome.apply(simplify_outcome)

# Normalizing the data by dividing the counts of each 2nd outcome by the total counts of each group of 1st outcomes
build_norm = (df_int.groupby('build_name')['Outcome'].value_counts() / df_int.groupby('build_name')['Outcome'].agg(len))
build_norm.name = 'Percent'
build_norm = pd.DataFrame(build_norm).reset_index()
build_norm.Percent = build_norm.Percent * 100
build_norm['build_name'] = build_norm['build_name'].apply(lambda x: int(x[-3:]))

# Plot normalized data
dims = (16, 6)
fig1, ax = plt.subplots(figsize=dims)
p = sns.barplot(data=build_norm, x='build_name', y='Percent', hue='Outcome',palette="Spectral").set_title("Normalized Outcomes by Build")
plt.xlabel("Build Number")
fig1.show()
input("move to next image")
plt.close()
fig1.savefig('../docs/Overall/norm_build.png')

# Caluclate the raw counts
build_raw = (df_int.groupby('build_name')['status'].value_counts())
build_raw.name = 'Count'
build_raw = pd.DataFrame(build_raw).reset_index()
build_raw['build_name'] = build_raw['build_name'].apply(lambda x: int(x[-3:]))

# Plot the raw counts
fig2, ax = plt.subplots(figsize=dims)
g = sns.barplot(data=build_raw, x='build_name', y='Count', hue='status',palette="Spectral")

# g.set_yscale('log')
plt.title("Outcome Counts by Build")
plt.xlabel("Build Number")
fig2.show()
input("move to next image")
plt.close()
fig2.savefig('../docs/Overall/raw_build.png')

# for i,build in build_norm.groupby('build_name'):
#     print(build,'\n')

## GENERATE SANKEY DIAGRAMS

def add_branch(source,target,amount,nodes,links):
    if amount == 0:
        return nodes,links
    nodes += [source,target]
    links.append([source,target,amount])
    return nodes,links

def gen_sankey(nodes,links,sankey={'nodes':[],'links':[]}):
    node_dict = dict([[y,x] for x,y in enumerate(pd.Series(nodes).unique())])
    print(node_dict)
    for name in node_dict.keys():
        sankey['nodes'].append({'name' : name})
    for source,target,value in links:
        sankey['links'].append({
            "source":node_dict[source],
            "target":node_dict[target],
            "value":value
        })
    return sankey

desired_dfs = [df[1] for df in data_b.groupby('Order_number')]
df_names = ["Order_"+str(num+1).zfill(3) for num in range(len(data_b.groupby('Order_number')))]
print(df_names)

desired_dfs.append(data_b)
df_names.append('Overall')

for name,df in zip(df_names,desired_dfs):
    sankey = {}

    nodes = ['Total_ordered']
    links = []
    sankey['nodes'] = []
    sankey['links'] = []

    # Orders
    total = len(df)
    orders = pd.DataFrame(df.Order_number.value_counts())
    order_names = ['Order_'+str(order).zfill(3) for order in orders.index.tolist()]
    order_counts = [int(count) for count in orders.Order_number.tolist()]

    for order,count in zip(order_names,order_counts):
        nodes,links = add_branch(order,'Total_ordered',count,nodes,links)

    # Synthesizing
    syn = len(df[df.status == 'ordered'])
    nodes,links = add_branch('Total_ordered','Synthesizing',syn,nodes,links)


    # Abandoned
    abandoned = len(df[df.status == 'abandoned'])
    nodes,links = add_branch('Total_ordered','Abandoned',abandoned,nodes,links)

    # Received
    received = total-abandoned-syn
    nodes,links = add_branch('Total_ordered','Received',received,nodes,links)

    # Attempted
    not_attempted = len(df[df.status == 'received'])
    nodes,links = add_branch('Received','Not_attempted',not_attempted,nodes,links)
    attempted = received-not_attempted
    nodes,links = add_branch('Received','Attempted',attempted,nodes,links)

    # Outcomes
    data_att = df[df.status != 'ordered']
    data_att = data_att[data_att.status != 'abandoned']
    data_att = data_att[data_att.status != 'received']
    outcomes = pd.DataFrame(data_att.status.value_counts())
    out = outcomes.index.tolist()
    count = outcomes.status.tolist()
    for out,count in zip(out,count):
        print(out,count)
        nodes,links = add_branch('Attempted',out,count,nodes,links)

    sankey = gen_sankey(nodes,links,sankey=sankey)

    path = '{}/docs/{}'.format(BASE_PATH,name)
    if os.path.exists(path):
        print("Directory for {} already exists".format(name))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)
        print("Making directory for {}".format(name))

    with open("{}/sankey.json".format(path),"w+") as json_file:
        json.dump(sankey,json_file,indent=2)

    copyfile('../docs/sankey/sankey.html','{}/sankey.html'.format(path))

    date = str(datetime.now()).split(" ")[0]
    if name == 'Overall':
        with open('../docs/sankey/index.md','r') as md_file:
            data = md_file.read()
            print(data)
            data = data.replace('[DATE]',date)
            print(data)
        with open('{}/index.md'.format(path),'w') as md_file:
            md_file.write(data)
    else:
        with open('../docs/sankey/order.md','r') as md_file:
            data = md_file.read()
            print(data)
            data = data.replace('[ORDER]',name)
            data = data.replace('[DATE]',date)
            print(data)
        with open('{}/index.md'.format(path),'w') as md_file:
            md_file.write(data)

## GENERATE THE FAILURES TABLE

fail2 = df_parts[df_parts.Attempt_1_Outcome_G == 'cloning_failure']
fail2 = fail2[fail2.Attempt_2_Outcome_G == 'cloning_failure']
fail2_limit = fail2[['part_id','part_name','Length','Fragments','part_type']]
fail2_limit

string = '| Part ID | Gene Name | Sequence Length | # of Fragments | Part Type |\n'
string += '| ------------- | ------------- | :-------------: | :-------------: | ------------- |\n'
for i,row in fail2_limit.iterrows():
    string += '| {} | {} | {} | {} | {} |\n'.format(row.part_id,row.part_name,row.Length,row.Fragments,row.part_type)
total_ordered_bp = df_parts.Length.sum()
total_ordered_g = len(df_parts)

rec = ['ordered','abandoned']
total_received = df_parts[df_parts.apply(lambda row: row.status not in rec, axis=1)]
total_received_bp = total_received.Length.sum()
total_received_g = len(total_received)

total_built = df_parts[df_parts.status == 'sequence_confirmed']
total_built_bp = total_built.Length.sum()
total_built_g = len(total_built)

rework = ['build008','build010']
normal_builds = df_parts[df_parts.apply(lambda row: row.Attempt_1_Build not in rework, axis=1)]
total_attempted = normal_builds[normal_builds.Attempt_1_Outcome_G != 'N/A']
total_built_first = normal_builds[normal_builds.Attempt_1_Outcome_G == 'sequence_confirmed']
total_attempted_g = len(total_attempted)
total_built_first_g = len(total_built_first)
success_rate = round((total_built_first_g / total_attempted_g)*100)


overall = '| Amount Ordered | Amount Received  | Amount Built | Cloning Success Rate |\n'
overall += '| :-------------: | :-------------: | :-------------: | :-------------: |\n'
overall += '| {}bp / {}genes | {}bp / {}genes | {}bp / {}genes | {}% |'.format(total_ordered_bp,total_ordered_g,total_received_bp,total_received_g,total_built_bp,total_built_g,success_rate)

with open('../docs/Overall/index.md','r') as md_file:
    data = md_file.read()
    print(data)
    data = data.replace('[Failure Table]',string)
    data = data.replace('[General Table]',overall)
    print(data)
with open('../docs/Overall/index.md'.format(path),'w') as md_file:
    md_file.write(data)










#
