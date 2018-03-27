import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math

multi_attempt = []
outcomes = []
empty = []
status = []
original = []

with open("./json/template.json","r") as temp_file:
    temp = json.load(temp_file)

for file in glob.glob('../../data/*/*.json'):
    with open(file,'r') as json_file:
        data = json.load(json_file)

    try:
        data['version'] == '2.0.0'
        continue
    except:
        print(data['gene_id'])

    temp['gene_id'] = data['gene_id']
    temp['author']['name'] = data['author']['name']
    temp['author']['email'] = data['author']['email']
    temp['author']['affiliation'] = data['author']['affiliation']
    temp['author']['orcid'] = data['author']['orcid']

    temp['info']['documentation']['gene_name'] = data['gene_name']
    temp['info']['documentation']['what'] = data['project_description']

    if data["info"]["database_links"] == [""]:
        temp['info']['documentation']['database_links'] = []
    else:
        temp['info']['documentation']['database_links'] = data["info"]["database_links"]

    if data['info']['type']['part_type'] == "CDS":
        data['info']['type']['part_type'] = "cds"
    temp['info']['gene_metadata']["cloning"]['part_type'] = data['info']['type']['part_type']


    if data['info']['type']['build_type'] == "10K_MoClo-EntryCDS-BbsI":
        temp['info']['gene_metadata']["cloning"]['cloning_enzyme'] = "BbsI"
    elif data['info']['type']['build_type'] == "10K_MoClo-EntryCDS-BtgZI":
        temp['info']['gene_metadata']["cloning"]['cloning_enzyme'] = "BtgZI"
    else:
        temp['info']['gene_metadata']["cloning"]['cloning_enzyme'] = data['info']['type']['build_type']

    temp['info']['gene_metadata']["cloning"]['retrieval_enzyme'] = "BsaI"

    temp['info']['gene_metadata']["cloning"]['optimize'] = True

    # Starts the state at submitted
    temp['status']['current_status'] = 'submitted'


    temp['info']['gene_metadata']['IP']['submission_number'] = data['info']['IP']['submission_number']
    temp['info']['gene_metadata']['IP']['results'] = data['info']['IP']['results']

    temp['info']['safety'] = ""

    temp['info']['order_number'] = data['info']['order_number']
    if temp['info']['order_number'] != "":
        temp['status']['current_status'] = 'ordered'
    if data['status']['abandoned'] == True:
        temp['status']['current_status'] = 'synthesis_abandoned'
    if data['status']['build_ready'] == True:
        temp['status']['current_status'] = 'received'


    temp['sequence']['original_sequence'] = data['sequence']['original_sequence'].upper()
    temp['sequence']['optimized_sequence'] = data['sequence']['optimized_sequence'].upper()
    temp['sequence']['fragment_sequences'] = data['sequence']['fragment_sequences']

    # Reset the build attempts
    temp['status']['build_attempts'] = [{
        "building": "",
        "build_well": "",
        "build_date": "",
        "build_number": "",
        "build_outcome": "",
        "forward_read": "",
        "reverse_read": ""
    }]

    for num,attempt in enumerate(reversed(data['status']['build_attempts'])):
        multi = False
        if num > 0:
            temp['status']['build_attempts'].append({
    			"build_number": "",
    			"build_well": "",
    			"build_date": "",
    			"build_outcome": "",
    			"forward_read": "",
    			"reverse_read": ""
            })
            multi = True
        temp['status']['build_attempts'][num]['build_number'] = attempt['build_number']
        temp['status']['build_attempts'][num]['build_well'] = attempt['build_well']

        try:
            temp['status']['build_attempts'][num]['build_date'] = attempt['build_date']
        except:
            temp['status']['build_attempts'][num]['build_date'] = '2018-02-14 17:30:11'

        temp['status']['build_attempts'][num]['build_number'] = attempt['build_number']
        temp['status']['build_attempts'][num]['build_outcome'] = attempt['build_outcome']


        if attempt['build_outcome'] == 'Good_sequence' or attempt['build_outcome'] == "Perfect":
            temp['status']['build_attempts'][num]['build_outcome'] = 'sequence_confirmed'
            temp['status']['current_status'] = 'sequence_confirmed'
        elif 'utat' in attempt['build_outcome']:
            temp['status']['build_attempts'][num]['build_outcome'] = 'cloning_mutation'
            temp['status']['current_status'] = 'cloning_mutation'
        elif attempt['build_outcome'] == "Split Reads" or attempt['build_outcome'] == "Partial" or attempt['build_outcome'] == "Incomplete":
            temp['status']['build_attempts'][num]['build_outcome'] = 'cloning_recombination'
            temp['status']['current_status'] = 'cloning_failure'
        elif attempt['build_outcome'] == 'Failed':
            temp['status']['build_attempts'][num]['build_outcome'] = 'robot_failure'
            temp['status']['current_status'] = 'cloning_failure'
        elif attempt['build_outcome'] == 'Original Vector Sequence' or attempt['build_outcome'] == 'Original_vector_sequence':
            temp['status']['build_attempts'][num]['build_outcome'] ='original_vector'
            temp['status']['current_status'] = 'cloning_failure'
        elif attempt['build_outcome'] == "":
            temp['status']['build_attempts'][num]['build_outcome'] ='empty'
            empty.append(temp['gene_id'])
        elif attempt['build_outcome'] == 'In_process':
            temp['status']['build_attempts'][num]['build_outcome'] = 'building'
            temp['status']['current_status'] = 'building'
        elif attempt['build_outcome'] == 'Bad_reads':
            temp['status']['build_attempts'][num]['build_outcome'] = 'sequence_failure'
            temp['status']['current_status'] = 'sequence_failure'
        elif "nknown" in attempt['build_outcome'] or attempt['build_outcome'] == "CHECK":
            temp['status']['build_attempts'][num]['build_outcome'] = 'unknown_sequence'
            temp['status']['current_status'] = 'unknown_sequence'
        else:
            temp['status']['build_attempts'][num]['build_outcome'] = attempt['build_outcome']

        outcomes.append(temp['status']['build_attempts'][num]['build_outcome'])

    if data['status']['build_complete'] == True or data['status']["build_complete"] == "Good_sequence":
        temp['status']['current_status'] = 'sequence_confirmed'

    temp['location']['fragments'] = data['location']['fragments']
    # temp['location']['cloned'] = data['cloned']

    temp['dates']['submitted'] = data['dates']['submitted']
    temp['dates']['ordered'] = data['dates']['ordered']
    temp['dates']['received'] = data['dates']['build_ready']
    temp['dates']['complete'] = data['dates']['complete']

    status.append(temp['status']['current_status'])
    original.append(data["status"]['build_complete'])

    with open(file,'w') as json_file:
        json.dump(temp,json_file,indent=2)

    # print(temp)
    # input()

# print(empty)
# print(pd.Series(status).value_counts())
# print("previous")
# print(pd.Series(original).value_counts())
