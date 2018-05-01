import numpy as np
import pandas as pd

import json
import os
import glob
import re

orders = []

for file in glob.glob("../../data/*/*.json"):
    print(file)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    orders.append(str(data['info']['order_number']))

order_nums = pd.Series(orders)

order_counts = order_nums.value_counts()

print(order_counts)
