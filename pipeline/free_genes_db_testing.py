import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import glob
import json
import pandas as pd
from config import *

## Begin using sqlite as a local database
engine = create_engine('sqlite:///:memory:', echo=False)

inspector = inspect(engine)

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
## Builds a many to many relationship between builds and parts
build_part = Table('build_part', Base.metadata,
    Column('build_id', ForeignKey('builds.id'), primary_key=True),
    Column('part_id', ForeignKey('parts.id'), primary_key=True)
)

class Part(Base):
    '''Descibes complete sequences to be synthesized and cloned'''
    __tablename__ = 'parts'

    id = Column(Integer, primary_key=True)
    part_id = Column(String) # Holds the id previously given to it
    part_type = Column(String) # For fragmentation
    organism = Column(String) # For codon optimization
    seq = Column(String)
    status = Column(String) # Gives the most updated status of this part

    # A part can have many fragments and one fragment can have many parts
    fragments = relationship('Fragment',
                            secondary=part_frag,
                            back_populates='parts')

    # A build can include many different parts and parts can be attempted many times
    builds = relationship('Build',
                            secondary=build_part,
                            back_populates='parts')

    # Allows the wells within different plates to link to the part inside
    wells = relationship("Well",back_populates='parts')

    def change_status(self,status):
        possible = ['submitted','optimized','ordered','synthesis_abandoned','received','building',
                   'sequencing','sequence_confirmed','cloning_mutation','sequence_failure','cloning_failure',
                   'cloning_abandoned']
        if status not in possible:
            print("Not a possible status")
        else:
            self.status = status

class Fragment(Base):
    '''Describes a single fragment to be synthesized'''
    __tablename__ = 'fragments'

    id = Column(Integer, primary_key=True)
    fragment_name = Column(String)
    location = Column(String)
    cloning_enzyme = Column(String) # Enzyme used in initial cloning
    retrieval_enzyme = Column(String) # Enzyme used to pull it out of dest vector
    syn_yield = Column(Integer) # Yield of DNA from synthesis (in ng)
    seq = Column(String) # Contains the complete sequence including the overhangs for cloning

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
    # Takes an address like 'A1' to indicate its position in the specified plate
    address = Column(String)

    # One plate can have many wells but a well can only have one plate
    plate_id = Column(Integer,ForeignKey('plates.id'))
    plates = relationship("Plate", back_populates="wells")

    # One well can only have a single part inside of it if assembly/seq_plate
    part_id = Column(Integer,ForeignKey('parts.id'))
    parts = relationship("Part",back_populates='wells')

    # One well can only have a single fragment inside of it if syn_plate
    fragment_id = Column(Integer,ForeignKey('fragments.id'))
    fragments = relationship("Fragment",back_populates='wells')


    def __init__(self,plate_type,item,address,syn_yield='',vector='',trans_outcome='',for_read='',rev_read='',seq_outcome=''):
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

    # One build can have many plates, but one plate can only belong to a single build
    build_id = Column(Integer,ForeignKey('builds.id'))
    builds = relationship("Build", back_populates="plates")

    # One plate can have many wells but a well can only have one plate
    wells = relationship("Well", back_populates="plates")

    def __init__(self,items,plate_type,plate_name):
        # Increment through all of the possible well addresses
        well_list = well_addresses()
        self.plate_name = plate_name
        self.plate_type = plate_type
        self.counter = 0
        self.next_well = well_list[self.counter]

        # Generate a well for every part provided
        for item in items:
            print("made well: ",self.next_well)
            self.wells.append(Well(self.plate_type,item,self.next_well))
            self.counter += 1
            self.next_well = well_list[self.counter]

    def add_item(self,item,address,syn_yield='',vector='',trans_outcome='',for_read='',rev_read='',seq_outcome=''):
        '''Allows the user to set the specific wells to
        associate with each part'''
        self.wells.append(Well(self.plate_type,item,address,syn_yield=syn_yield,vector=vector,trans_outcome=trans_outcome,\
                              for_read=for_read,rev_read=rev_read,seq_outcome=seq_outcome))

class Build(Base):
    '''Describes a complete build'''
    __tablename__ = 'builds'

    id = Column(Integer, primary_key=True)
    date = Column(String) # Date that the build was conducted
    master_mix = Column(String) # The master mix that was used

    # A build can include many different parts and parts can be attempted many times
    parts = relationship('Part',
                            secondary=build_part,
                            back_populates='builds')

    # One build can have many plates, but one plate can only belong to a single build
    plates = relationship("Plate", back_populates="builds")

    # A single sequencing order could include several builds
    seq_orders_id = Column(Integer,ForeignKey('seq_orders.id'))
    seq_orders = relationship("Seq_order", back_populates="builds")

    def __init__(self,items):
        # Generate an assembly plate
        self.plates = [Plate(items,'assembly_plate','test_name')]
        print("made a plate")

    def add_item(self,item,address,vector='',trans_outcome=''):
        '''Allows the user to set the specific wells to
        associate with each part'''
        self.plates[0].wells.append(Well(self.plates[0].plate_type,item,address,vector=vector,trans_outcome=trans_outcome))


## Create and commit all of the tables
Base.metadata.create_all(engine)

Session = sessionmaker(bind=engine)
Session.configure(bind=engine)
session = Session()

## Take in plate maps from twist and generate fragment objects
twist1 = pd.read_csv("{}/plate_maps/O-001_A-001_0000.00.00.csv".format(BASE_PATH))
twist1['customer_line_item_id'] = twist1['customer_line_item_id'].str.strip()
plates_made = []
for index,row in twist1.iterrows():
    if row['Plate'] not in plates_made:
        if plates_made != []:
            session.add(current_plate)
        current_plate = Plate([],'syn_plate',plate_name=row['Plate'])
        plates_made.append(row['Plate'])
    current_plate.add_item(Fragment(fragment_name=row['customer_line_item_id'],seq=row['Insert Sequence']),row['Well'],syn_yield=row['Yield (ng)'])

# Verify that the synthesis plates objects work
for instance in session.query(Plate):
    print()
    print(instance.plate_name)
    for well in instance.wells:
        print(well.address,":",well.fragments.fragment_name,"-",well.syn_yield)

## Generate part and fragment objects from JSON database
j_counter = 0
for file in sorted(glob.glob("{}/data/*/*.json".format(BASE_PATH))):
    with open(file,"r") as json_file:
        data = json.load(json_file)
#     print(data['gene_id'])
    new = Part(part_id=data['gene_id'],
        part_type=data['info']['gene_metadata']['cloning']['part_type'],
        seq=data['sequence']['optimized_sequence'])
    frags = []
    for fragment in data['sequence']['fragment_sequences']:
        session.add(Fragment(
            fragment_name=fragment,
            location=data['location']['fragments'][fragment],
            seq=data['sequence']['fragment_sequences'][fragment],
            parts=[new]
            ))
#     print(new.part_id)
    session.add(new)

    if j_counter > 650:
        break
    j_counter += 1

# Test that the json file information is transferred
test_part = session.query(Part).filter(Part.part_id == 'BBF10K_000397').first()
print("Test part name: ",test_part.part_id,"Seq: \n",test_part.seq)

## Convert the build maps from csv's to SQL objects
build1 = pd.read_csv("{}/builds/build006/build006_2018-02-09 11:53:22-1.csv".format(BASE_PATH))

test_build = Build([])
for index,row in build1.iterrows():
    part_obj = session.query(Part).filter(Part.part_id == row['Gene']).first()
    test_build.add_item(part_obj,row['Destination'],vector='pOpen_v.1.1',trans_outcome='good')

session.add(test_build)

# Test that the json file information is transferred
test_part = session.query(Part).filter(Part.part_id == 'BBF10K_000397').first()
print("Test part in build {} well {}".format(test_part.wells[0].plates.builds.id,test_part.wells[0].address))

# Test the build plate object works:
for instance in session.query(Build):
    for plate in instance.plates:
        for well in plate.wells:
            print(well.parts.part_id,well.address,well.vector,well.trans_outcome)







#
