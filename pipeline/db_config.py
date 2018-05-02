import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import pandas as pd
from config import *

## Connect to the AWS server running the openfoundry database
conn_str = 'postgresql+psycopg2://openfoundry:freegenestomakegenesfree@freegenes-openfoundry.cwtlxuukykrr.us-east-1.rds.amazonaws.com:5432/openfoundry'
engine = sqlalchemy.create_engine(conn_str, echo=False)

## Begin using sqlite as a local database
# engine = create_engine('sqlite:///:memory:', echo=False)
# engine = create_engine('sqlite:///free_genes.db')

# inspector = inspect(engine)

Base = declarative_base()

## Builds a many to many relationship between parts and fragments
part_frag = Table('part_frag', Base.metadata,
    Column('part_id', ForeignKey('parts.id'), primary_key=True),
    Column('fragment_id', ForeignKey('fragments.id'), primary_key=True)
)
## Builds a many to many relationship between fragments and orders
frag_order = Table('frag_order', Base.metadata,
    Column('frag_id', ForeignKey('fragments.id'), primary_key=True),
    Column('twist_order_id', ForeignKey('twist_orders.id'), primary_key=True)
)

class Part(Base):
    '''Descibes complete sequences to be synthesized and cloned'''
    __tablename__ = 'parts'

    id = Column(Integer, primary_key=True)
    part_id = Column(String) # Holds the id previously given to it
    part_name = Column(String)
    part_type = Column(String) # For fragmentation
    cloning_enzyme = Column(String) # Enzyme used in initial cloning
    organism = Column(String) # For codon optimization
    seq = Column(String)
    status = Column(String) # Gives the most updated status of this part

    # A part can have many fragments and one fragment can have many parts
    fragments = relationship('Fragment',
                            secondary=part_frag,
                            back_populates='parts')

    # Allows the wells within different plates to link to the part inside
    wells = relationship("Well",back_populates='parts')

    def change_status(self,status):
        possible = ['submitted','optimized','ordered','synthesis_abandoned','received','trans_failure',
                   'cloning_failure','cloning_error','sequence_failure','cloning_mutation','building',
                   'sequencing','cloning_abandoned','sequence_confirmed']
        rank = dict(zip(possible,range(len(possible))))
        if status not in possible:
            print("Not a possible status")
        elif rank[status] > rank[self.status]:
            # print(status,rank[status],'is greater than',self.status,rank[self.status])
            self.status = status
        # else:
            # print("no change")
    def eval_status(self):
        status = [well.seq_outcome for well in self.wells if well.plates.plate_type == 'seq_plate']
        print(self.part_id,status)
        simple_status = []
        for s in status:
            if "mutation" in s:
                simple_status.append('cloning_mutation')
            elif "bad" in s:
                simple_status.append('sequence_failure')
            else:
                simple_status.append(s)

        possible = ['cloning_failure','cloning_error','sequence_failure','cloning_mutation','sequence_confirmed']
        rank = dict(zip(possible,range(len(possible))))
        try:
            new_status = possible[max([rank[r] for r in simple_status])]
            print(new_status)
            self.status = new_status
        except:
            print("skipped")

class Fragment(Base):
    '''Describes a single fragment to be synthesized'''
    __tablename__ = 'fragments'

    id = Column(Integer, primary_key=True)
    fragment_name = Column(String)
    retrieval_enzyme = Column(String) # Enzyme used to pull it out of dest vector
    syn_yield = Column(Integer) # Yield of DNA from synthesis (in ng)
    seq = Column(String) # Contains the complete sequence including the overhangs for cloning
    cloning_method = Column(String) # Details which vector class it must be cloned into

    # A part can have many fragments and one fragment can have many parts
    parts = relationship('Part',
                            secondary=part_frag,
                            back_populates='fragments')

    # Many fragments are included in an order and fragments could be reordered
    twist_orders = relationship('Twist_order',
                            secondary=frag_order,
                            back_populates='fragments')

    # Allows the wells within different plates to link to the fragment inside
    wells = relationship("Well",back_populates='fragments')



class Twist_order(Base):
    '''Describes an order sent out to Twist for synthesis'''
    __tablename__ = 'twist_orders'

    id = Column(Integer, primary_key=True)
    date = Column(String) # Date that the order was sent out
    invoice = Column(String) # Invoice number for order

    # Many fragments are included in an order and fragments could be reordered
    fragments = relationship('Fragment',
                            secondary=frag_order,
                            back_populates='twist_orders')

class Seq_order(Base):
    '''Describes an order sent out to Elim for sequencing'''
    __tablename__ = 'seq_orders'

    id = Column(Integer, primary_key=True)
    order = Column(String) # Order number assigned by Elim

    # A single sequencing order could include several builds
    builds = relationship("Build", back_populates="seq_orders")

    def __init__(self,builds):
        self.builds.append(builds) # Creates the relationship between seq_order and builds

class Well(Base):
    '''Describes a well within a plate'''
    __tablename__ = 'wells'

    id = Column(Integer,primary_key=True)
    address = Column(String) # Takes an address like 'A1'
    plate_type = Column(String) # Uses the same types listed in the Plate class

    # Columns for "syn_plate" wells
    syn_yield = Column(Integer) # The yield specified in the plate map from twist
    concentration = Column(Integer) # Number of fmoles in 2µL
    volume = Column(Integer) # Volume in the well in µL

    # Columns for "assembly_plate" wells
    vector = Column(String) # The vector used in the cloning reaction
    trans_outcome = Column(String)

    # Columns for "seq_plate" wells
    for_read = Column(String)
    rev_read = Column(String)
    seq_outcome = Column(String)
    misplaced = Column(String) # To indicate if it is the result of the BLAST alignment

    # One plate can have many wells but a well can only have one plate
    plate_id = Column(Integer,ForeignKey('plates.id'))
    plates = relationship("Plate", back_populates="wells")

    # One well can only have a single part inside of it if assembly/seq_plate
    part_id = Column(Integer,ForeignKey('parts.id'))
    parts = relationship("Part",back_populates='wells')

    # One well can only have a single fragment inside of it if syn_plate
    fragment_id = Column(Integer,ForeignKey('fragments.id'))
    fragments = relationship("Fragment",back_populates='wells')


    def __init__(self,plate_type,item,address,syn_yield='',vector='',\
                trans_outcome='',for_read='',rev_read='',seq_outcome='',misplaced=''):
        self.plate_type = plate_type
        self.address = address
        if self.plate_type == 'syn_plate':
            self.fragments = item # Links the well to the fragment inside
            self.syn_yield = syn_yield
        elif self.plate_type == 'assembly_plate':
            self.parts = item # Links the well to the part inside
            self.vector = vector
            self.trans_outcome = trans_outcome
        elif self.plate_type == 'seq_plate':
            self.parts = item # Links the well to the part inside
            self.for_read = for_read
            self.rev_read = rev_read
            self.seq_outcome = seq_outcome
            self.misplaced = misplaced
        else:
            print(plate_type)
            input("plate_type didn't match")

def well_addresses():
    '''Generates a list of well address A1-H12'''
    letter = ["A","B","C","D","E","F","G","H"]
    number = ["1","2","3","4","5","6","7","8","9","10","11","12"]
    target_well = []
    temp_well = 0
    for n in number:
        for l in letter:
            temp_well = l + n
            target_well.append(temp_well)
    return target_well


class Plate(Base):
    '''Describes a 96 well plate'''
    __tablename__ = 'plates'

    id = Column(Integer, primary_key=True)
    plate_type = Column(String) # Takes in syn_plate, assembly_plate, or seq_plate
    plate_name = Column(String)
    resuspended = Column(String) # 'resuspended' or 'not resuspended'
    plated = Column(String) # 'plated' or 'not_plated'
    assessed = Column(String) # 'assessed' or 'not_assessed'

    # One build can have many plates, but one plate can only belong to a single build
    build_id = Column(Integer,ForeignKey('builds.id'))
    builds = relationship("Build", back_populates="plates")

    # One plate can have many wells but a well can only have one plate
    wells = relationship("Well", back_populates="plates")

    def __init__(self,items,plate_type,plate_name):
        # Increment through all of the possible well addresses
        self.plate_name = plate_name
        self.plate_type = plate_type
        if self.plate_type == 'syn_plate':
            self.resuspended = 'not_resuspended'
        if self.plate_type == 'assembly_plate':
            self.plated = 'not_plated'
        # Generate a well for every part provided
        for item in items:
            self.add_item(item)

    def add_item(self,item,address='',syn_yield='',vector='',trans_outcome='',for_read='',rev_read='',seq_outcome=''):
        '''Allows the user to set the specific wells to
        associate with each part'''
        used_wells = [well.address for well in self.wells]
        all_wells = well_addresses()
        empty_wells = [well for well in all_wells if well not in used_wells]
        if len(used_wells) == len(all_wells):
            print("Plate is full")
            return
        elif address == '':
            address = empty_wells[0]
        elif address not in empty_wells:
            print("Address is already taken")
            return

        self.wells.append(Well(self.plate_type,item,address,syn_yield=syn_yield,\
                            vector=vector,trans_outcome=trans_outcome,\
                            for_read=for_read,rev_read=rev_read,seq_outcome=seq_outcome))
    def resuspend(self):
        self.resuspended = 'resuspended'

    def plate(self):
        self.plated = 'plated'

    def assess(self):
        self.assessed = 'assessed'


class Build(Base):
    '''Describes a complete build'''
    __tablename__ = 'builds'

    id = Column(Integer, primary_key=True)
    build_name = Column(String)
    status = Column(String) # Possible states: 'building','sequencing','complete'
    date = Column(String) # Date that the build was conducted
    master_mix = Column(String) # The master mix that was used

    # One build can have many plates, but one plate can only belong to a single build
    plates = relationship("Plate", back_populates="builds")

    # A single sequencing order could include several builds
    seq_orders_id = Column(Integer,ForeignKey('seq_orders.id'))
    seq_orders = relationship("Seq_order", back_populates="builds")

    def __init__(self,items,build_name='',rxn_per_plate=96):
        # Generate an assembly plate
        group_parts = [items[n:n+rxn_per_plate] for n in range(0, len(items), rxn_per_plate)]
        self.build_name = build_name
        self.status = 'building'
        for num,group in enumerate(group_parts):
            plate_name = 'Assembly_' + str(num)
            self.plates.append(Plate(group,'assembly_plate',plate_name))

    def add_item(self,item,address,vector='',trans_outcome=''):
        '''Allows the user to set the specific wells to
        associate with each part'''
        self.plates[0].wells.append(Well(self.plates[0].plate_type,item,address=address,\
                    vector=vector,trans_outcome=trans_outcome))

## Create and commit all of the tables
Base.metadata.create_all(engine)

Session = sessionmaker(bind=engine)
Session.configure(bind=engine)
session = Session()



#
