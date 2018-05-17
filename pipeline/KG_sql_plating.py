from opentrons import robot, containers, instruments
import argparse
import sys

import sqlalchemy
from sqlalchemy import create_engine,Column,Integer,String,ForeignKey,Table,Text,inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker,relationship

import os
import re
import math
import glob
import json
import numpy as np
import pandas as pd
from datetime import datetime


from config import *
from db_config import *
session,engine = connect_db()

from opentrons_functions import locations as ot_locations , recalibration as recal, freegenes_funcs as ff

PIPELINE_PATH = BASE_PATH + "/pipeline"
BUILDS_PATH = BASE_PATH + "/builds"
DATA_PATH = BASE_PATH + "/data"

BACKBONE_PATH = BASE_PATH + "/sequencing_files/popen_v1-1_backbone.fasta"
DICTIONARY_PATH = PIPELINE_PATH + "/testing/data_testing/10K_CDS.csv"

def db_choose_plate():
    assembly_options = []
    assemblies = []
    print("Choose which plate you would like to transform/plate:")
    for index,assembly in enumerate(session.query(Plate).join(Build,Plate.builds).filter(Plate.plated == 'not_plated').order_by(Build.build_name)):
        assembly_options.append("{}. {}".format(index,assembly.builds.build_name))
        assemblies.append((assembly.builds.build_name, assembly))
    target_plate = ff.option_list(assemblies)[1]
    build_map = target_plate.wells
    return len(build_map)

def agar_plate_preparation(num_rows, locations, annotations, portion=1, trans_per_plate=3):
    print("num rows: ",num_rows)
    num_plates = num_rows // trans_per_plate

    agar_plate_names = []
    for plate_num in range(num_plates):
        if portion == 2:
            plate_num += 2
        current_plate = "Buildxxx" + "_p" + str(plate_num + 1)
        agar_plate_names.append(current_plate)
    print(agar_plate_names)
    print("You will need {} agar plates".format(len(agar_plate_names)))
    ot_locations.show_opentrons_locations(locations, annotations)
    return agar_plate_names

def OT1_plate(num_reactions, simulate=True, num_dilutions=4, plate_vol=7.5, dilution_vol=9, waste_vol=2.5, plating_row=0, media_per_tube = 150):
    # Setup parameters
    locations = {
        "A3":"tiprack-200ul",
        "E3":"tiprack-10ul",
        "E2":"tiprack-10ul",
        "E1":"tiprack-10ul",
        "D1":"point",
        "C3":"PCR-strip-tall",
        "C2":"96-PCR-tall",
        "B1":"tube-rack-2ml",
        "D2":"96-deep-well",
        "D3":"96-deep-well",
        }
    annotations = {
        "D1":"Trash",
        "C2":"transformation_plate",
        "D2":"Agar_plate_1",
        "D3":"Agar_plate_2"
        }

    # Setup robot
    location_data = ot_locations.load_robot_from_locations(locations, annotations)

    p10 = location_data["p10"]
    p200 = location_data["p200"]
    p200_tipracks = location_data["p200_tipracks"]
    p10_tipracks = location_data["p10_tipracks"]

    transformation_plate = location_data["container_dict"]["transformation_plate"]
    trash = location_data["container_dict"]["Trash"]
    centrifuge_tube = location_data["container_dict"]["B1"]
    master = location_data["container_dict"]["C3"]

    # Activate simulation
    recal.simulation(simulate=True)

    # Setup agar
    num_rows = math.ceil(num_reactions / 8)

    def aliquot_LB(lb_location='A1'):
        robot.home()
        # Aliquot the LB into the PCR tube strip for the dilutions
        p200.pick_up_tip()
        for well in range(8):
            print("Transferring {}ul to tube {}".format(media_per_tube,well))
            p200.transfer(media_per_tube, centrifuge_tube[lb_location].bottom(),master.wells(well).bottom(),new_tip='never')
        p200.drop_tip()
        if lb_location == 'A1':
            return plating_manager()

    # Iterate through each row of the transformation plate
    def plate_cells(half="First",plating_row=plating_row):
        redo = True
        robot.home()
        current_plate = location_data["container_dict"]["Agar_plate_1"]
        agar_plate_names = agar_plate_preparation(num_rows, locations, annotations)
        plate = agar_plate_names[0]
        for trans_row in range(num_rows):
            if trans_row == 0:
                continue
            # Iterates through each dilution
            for dilution in range(num_dilutions):
                # Resets to switch to the next plate
                if plating_row == 12:
                    p200.pick_up_tip()
                    for well in range(8):
                        aliquot_LB(lb_location='B1')
                    p200.drop_tip()
                    current_plate = location_data["container_dict"]["Agar_plate_2"]
                    plate = agar_plate_names[1]
                    print("changing to :", plate)
                    plating_row = 0
                    redo = True

                print("trans_row",trans_row, "dilution",dilution,"plate",plate,"plate row",plating_row)
                p10.pick_up_tip()
                print("Diluting cells in row {} with {}ul".format(trans_row, dilution_vol))
                p10.transfer(dilution_vol, master['A1'].bottom(), transformation_plate.rows(trans_row).bottom(),new_tip='never',mix_before=(1,9))
                print("Plating {}ul from transformation row {} onto {} in row {}".format(plate_vol,trans_row,plate,plating_row))
                p10.aspirate(plate_vol,transformation_plate.rows(trans_row).bottom())
                p10.dispense(plate_vol, current_plate.rows(plating_row).bottom())

                if redo:
                    redo = recal.change_height(current_plate,current_plate.rows(plating_row)[0],location_data)
                    p10.blow_out()

                print("Discard {}ul from transformation row {} into waste tube".format(waste_vol,trans_row))
                p10.aspirate(waste_vol,transformation_plate.rows(trans_row).bottom())
                p10.drop_tip()
                plating_row += 1
        print("Rehoming")
        robot.home()
        return plating_manager()

    def plating_manager():
        options = ("Aliquot LB","Plate first half","Plate second half","End simulation","Exit")
        choice = ff.option_list(options)

        if choice == "Aliquot LB":
            return aliquot_LB()
        elif choice == "Plate first half":
            return plate_cells()
        elif choice == "Plate second half":
            return plate_cells()
        elif choice == "End simulation":
            recal.simulation(simulate=False)
            return plating_manager()
        elif choice == "Exit":
            sys.exit()
    plating_manager()



OT1_plate(int(input("How many reactions do you wish to run? :" )))


#
