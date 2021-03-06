import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

from config import *
from db_config import *
session,engine = connect_db()


for part in session.query(Part).order_by(Part.part_id):
    print(part.part_id)
    part.eval_status()

commit = int(input("Commit changes (1-yes, 2-no): "))
if commit == 1:
    session.commit()
