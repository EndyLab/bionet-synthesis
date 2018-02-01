import numpy as np
import pandas as pd

# Create a dictionary to link the gene name to the corresponding id number
data = pd.read_csv("./data_testing/10K_CDS.csv")
dictionary = dict(zip(data['gene_name'], data['idnum']))

id_nums = []

PATH = "../../raw_files/distribution_maps/JCVI_MMSYN1_1/MMSYN1_plate1.csv"

data = pd.read_csv(PATH)

for index, row in data.iterrows():
    print(row["Gene Name"])
    if row["Gene Name"] == "Empty":
        id_nums.append("Empty")
        continue
    id_num = dictionary[row["Gene Name"]]
    id_nums.append(id_num)


data["Gene ID"] = id_nums

print(data)

data = data[["Gene ID","Gene Name","Well"]]

data.to_csv("./build000.csv")
