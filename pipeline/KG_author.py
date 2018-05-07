import json
import glob

author = input("What author? : ").lower()
for file in glob.glob('../data/*/*.json'):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    file_author = data['author']['name'].lower()
    if author in file_author:
        print(data['gene_id'])


    
