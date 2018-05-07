import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

from config import *
import db_config
session,engine = db_config.connect_db()


for part in session.query(Part).order_by(Part.id):
    print(part.part_id)
    part.eval_status()

commit = int(input("Commit changes (1-yes, 2-no): "))
if commit == 1:
    session.commit()
