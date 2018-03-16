import glob
import json

for file in glob.glob('../../data/*/*.json'):
    with open(file,"r") as json_file:
        data = json.load(json_file)
    print(data['gene_id'])
    print(data["status"]["build_attempts"][0]["build_outcome"])
    if data["status"]["build_attempts"][0]["build_outcome"] == 'In_process':
        data["status"]["build_complete"] = "In_process"
    else:
        print("good")
    with open(file,"w+") as json_file:
        json.dump(data,json_file,indent=2)
