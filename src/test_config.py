import os,sys,inspect
import pandas as pd
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)

print(currentdir,parentdir)
#
# from config import *
#
# print(BASE_PATH)
#

# import sys
# sys.path[0] = '/Users/conarymeyer/bionet-synthesis/pipeline/configuration/'
# from config import *
# print(BASE_PATH)
#
# sys.path[0] = '/Users/conarymeyer/bionet-synthesis/pipeline/'
#
# import ot_functions as ot
#
# data = pd.read_csv('../all_deck_layout.csv')
# print(data)

# import sys
# sys.path[0] = '/Users/conarymeyer/bionet-synthesis/pipeline/configuration/'
# print(sys.path[0])
# from config import *
# print(BASE_PATH)
