from opentrons import robot, containers, instruments
import numpy as np
import pandas as pd
import getch
import shutil
import os
import sys

def initialize_pipettes(p10_tipracks,p10s_tipracks,p200_tipracks,trash):
    # Declare all of the pipettes
    p10 = instruments.Pipette(
        axis='a',
        max_volume=10,
        min_volume=0.5,
        tip_racks=p10_tipracks,
        trash_container=trash,
        channels=8,
        name='p10-8',
        aspirate_speed=400,
        dispense_speed=800
    )

    p10s = instruments.Pipette(
        axis='a',
        max_volume=10,
        min_volume=0.5,
        tip_racks=p10s_tipracks,
        trash_container=trash,
        channels=1,
        name='p10-8s',
        aspirate_speed=400,
        dispense_speed=800
    )

    p200 = instruments.Pipette(
        axis='b',
        max_volume=200,
        min_volume=20,
        tip_racks=p200_tipracks,
        trash_container=trash,
        channels=1,
        name='p200-1',
        aspirate_speed=400,
        dispense_speed=800
    )
    return p10,p10s,p200

def display_deck(robot):
    df = pd.DataFrame(np.zeros((3,5)), columns=['A','B','C','D','E'], index=['3','2','1'])
    df.loc[:,:] = "---"

    for slot in robot.deck:
        for child in slot.get_children_list():
            print(slot.get_name()[0],slot.get_name()[1],child.get_name())
            df.loc[slot.get_name()[1],slot.get_name()[0]] = child.get_name()
    print(df)

def print_layout(locations):

    # Generate an empty dataframe with the right shape
    layout_table = pd.DataFrame(np.zeros((3,5)), columns=['A','B','C','D','E'], index=['3','2','1'])
    layout_table.loc[:,:] = "---"

    # Fill in the data frame with the locations
    for obj in locations:
            layout_table.loc[locations[obj][1], locations[obj][0]] = obj

    # Displays the required plate map and waits to proceed
    print("\n Please arrange the items in the following configuration: \n")
    print(layout_table,"\n")
    input("Press enter to continue")

def change_speed(robot):
    robot.head_speed(5000)

def change_height(pipette,container,target,recalibrate=False):
    counter = 0
    z = 0
    print("Change height - s-g:up h-l:down x:exit")
    while True:
        c = getch.getch()
        if c == "s":
            print("Up 20mm")
            pipette.robot._driver.move(z=20,mode="relative")
            z += 20
        elif c == "d":
            print("Up 5mm")
            pipette.robot._driver.move(z=5,mode="relative")
            z += 5
        elif c == "f":
            print("Up 0.5mm")
            pipette.robot._driver.move(z=0.5,mode="relative")
            z += 0.5
        elif c == "g":
            print("Up 0.1mm")
            pipette.robot._driver.move(z=0.1,mode="relative")
            z += 0.1
        elif c == "h":
            print("Down 0.1mm")
            pipette.robot._driver.move(z=-0.1,mode="relative")
            z += -0.1
        elif c == "j":
            print("Down 0.5mm")
            pipette.robot._driver.move(z=-0.5,mode="relative")
            z += -0.5
        elif c == "k":
            print("Down 5mm")
            pipette.robot._driver.move(z=-5,mode="relative")
            z += -5
        elif c == "l":
            print("Down 20mm")
            pipette.robot._driver.move(z=-20,mode="relative")
            z += -20
        elif c == "x":
            print("Exit")
            break
        counter += 1
    pipette.calibrate_position((container,target.from_center(x=0, y=0, z=-1,reference=container)))

    if recalibrate:
        if counter > 1:
            print("Will recalibrate")
            redo = True
        else:
            print("Calibrated")
            redo = False
        return redo,z
    else:
        return z

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

def print_center(statement):
    columns = shutil.get_terminal_size().columns
    print('\n',statement.center(columns))

def request_info(statement,type='string',length=0):
    answer = input(statement)
    if answer == '':
        print("Please enter a value\n")
        return request_info(statement,type=type)
    elif type == 'int':
        try:
            int(answer)
            return int(answer)
        except:
            print("Invalid type\n")
            return request_info(statement,type=type)
    elif type == 'list':
        try:
            nums = [int(num) for num in answer.split(' ')]
            if len(nums) != length:
                print('Requires {} inputs'.format(length))
                return request_info(statement,type=type,length=length)
            return [int(num) for num in answer.split(' ')]
        except:
            print("Invalid type\n")
            return request_info(statement,type=type,length=length)
    else:
        return answer

def make_directory(path):
    dir_name = path.split("/")[-1]
    if os.path.exists(path):
        print("Directory {} already exists".format(dir_name))
    else:
        # Generates a new directory with the ID# as its name
        os.makedirs(path)
        print("Making directory for {}".format(dir_name))

def check_robot():
    try:
        robot_name = str(os.environ["ROBOT_DEV"][-5:])
    except:
        sys.exit("Not connected to a robot, run roboswitch <robot_name> to change the robot")
    robot_number = int(request_info("Run on this robot: {} ? 1-Yes, 2-No ".format(robot_name),type='int'))
    if robot_number == 1:
        print("Proceeding with run")
    else:
        sys.exit("Run `roboswitch <robot_name>` to change the robot")

def list_to_string(ls):
    string = ''
    for l in ls:
        string += "'{}',".format(l)
    return string[:-1]

def query_for_parts(status,enzyme,engine):
    query = "SELECT parts.part_id,parts.status,fragments.fragment_name,plates.plate_id,wells.address,wells.volume,plates.id FROM parts\
            INNER JOIN part_frag ON parts.id = part_frag.part_id\
            INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
            INNER JOIN wells ON fragments.id = wells.fragment_id\
            INNER JOIN plates on wells.plate_id = plates.id\
            WHERE parts.status IN ({})\
                AND parts.cloning_enzyme = '{}'".format(list_to_string(status),enzyme)

    return pd.read_sql_query(query, con=engine)

def query_for_plates(parts,engine):
    query = "SELECT parts.part_id,fragments.fragment_name,plates.plate_id,wells.address,wells.volume,plates.id,plates.plate_name FROM parts\
            INNER JOIN part_frag ON parts.id = part_frag.part_id\
            INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
            INNER JOIN wells ON fragments.id = wells.fragment_id\
            INNER JOIN plates on wells.plate_id = plates.id\
            WHERE parts.part_id IN ({})".format(list_to_string(parts))

    return pd.read_sql_query(query, con=engine)

def query_for_fragments(parts,engine):
    query = "SELECT parts.part_id,fragments.fragment_name FROM parts\
            INNER JOIN part_frag ON parts.id = part_frag.part_id\
            INNER JOIN fragments ON part_frag.fragment_id = fragments.id\
            WHERE parts.part_id IN ({})".format(list_to_string(parts))

    return pd.read_sql_query(query, con=engine)

def change_plates(locations,current_plates,source_slots):
    '''Allows the user to swap plates in the middle of a protocol'''
    source_plates = {}
    temp = locations
    plate_locations = list(zip(pd.unique(current_plates),source_slots[:len(current_plates)]))
    temp.update(plate_locations)
    print_layout(temp)
    for plate, slot in plate_locations:
        source_plates[plate] = containers.load('96-flat', slot)
    return source_plates

def make_gg_rxns(num_rxns,rxn_vol):
    '''Calculates the amount of each reagent to add to reach the desired master mix'''
    cutsmart = 1 * num_rxns
    atp = 1 * num_rxns
    vector = 0.25 * num_rxns
    ligase = 0.5 * num_rxns
    enzyme = 0.25 * num_rxns
    water = (rxn_vol - ((cutsmart + atp + vector + ligase + enzyme)/num_rxns)) * num_rxns
    master_mix = pd.DataFrame(
        {'Component':['H2O','Cutsmart','ATP','Vector','T4 Ligase','Restriction Enzyme','Total'],
        'Amount':[water,cutsmart,atp,vector,ligase,enzyme,rxn_vol*num_rxns]},
        columns=["Component","Amount"]
    )
    return master_mix



    #
