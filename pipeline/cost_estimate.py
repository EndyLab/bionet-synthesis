import pandas as pd
import glob
from config import *

file = glob.glob('{}/**/openfoundry_material_list.csv'.format(BASE_PATH),recursive=True)[0]
data = pd.read_csv(file)
cost = data.cost_per_build.sum()
number = int(input('Enter how many builds: '))
total = (number * cost)
print("The approximate cost is: ${}".format(total))
