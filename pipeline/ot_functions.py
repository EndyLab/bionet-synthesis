from opentrons import robot, containers, instruments
import numpy as np
import pandas as pd

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
    # return robot













    #
