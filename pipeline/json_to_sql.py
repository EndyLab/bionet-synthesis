import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import numpy as np
import pandas as pd
import glob
import json
import re

from config import *
from db_config import *

enzyme_dict = {}
counter = 0
for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
    counter += 1
    if counter == 800:
        break
    print(counter)
    with open(file,"r") as json_file:
        data = json.load(json_file)
    part = data['gene_id']
    enz = data['info']['gene_metadata']['cloning']['cloning_enzyme']
    part_dict = {part:enz}
    enzyme_dict.update(part_dict)

for part in session.query(Part).order_by(Part.part_id):
    part.cloning_enzyme = enzyme_dict[part.part_id]
    print(part.part_id,part.cloning_enzyme)

commit = int(input("Commit changes (1-yes, 2-no): "))
if commit == 1:
    session.commit()












#
