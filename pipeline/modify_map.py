import numpy as np
import pandas as pd

import json
import os
import glob
import re
import math
import datetime
from datetime import datetime

gene_names = []
wells = []

dist_map = pd.read_csv('JCVI_plate1.csv')

for index,row in dist_map.iterrows():
    full_name = str(row['Gene Name']).strip()
    gene_name = full_name[:-2]
    gene_names.append(gene_name)
    well = row['Well']
    wells.append(well)

new_data = pd.DataFrame({
    "Gene Name" : gene_names,
    "Well" : wells
})

new_data = new_data[["Gene Name","Well"]]
print(new_data)

new_data.to_csv('MMSYN1_plate1.csv')
