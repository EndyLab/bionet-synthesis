from opentrons import robot, containers, instruments
import numpy as np
import pandas as pd
import getch

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

def display_deck():
    df = pd.DataFrame(np.zeros((3,5)), columns=['A','B','C','D','E'], index=['3','2','1'])
    df.loc[:,:] = "---"

    for container, placeable in robot.containers().items():
        coord = list(placeable.get_parent().get_name())
        df.loc[coord[1],coord[0]] = placeable.get_name()

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











    #
