import pandas as pd
from opentrons import robot, containers, instruments

# Examples
locations = {
        "A3":"tiprack-200ul",
        "E3":"tiprack-10ul",
        "E2":"tiprack-10ul",
        "E1":"tiprack-10ul",
        "D1":"point",
        "C3":"PCR-strip-tall",
        "C2":"96-PCR-tall",
        "B1":"tube-rack-2ml"
        }

annotations = {
        "D1":"Trash"
        }

OT2_locations = {
        '1':'Output Plate',
        '2':'',
        '3':'',
        '4':'300uL Tips',
        '5':'DeepWell_Cells',
        '6':'H2O',
        '7':'300uL Tips',
        '8':'Dilutant',
        '9':'Lysis_buffer',
        '10':'',
        '11':''
        }

def load_robot_from_locations(locations, annotations, OT="OT1"):
    # Requires an annotated Trash
    if OT == "OT1":
        container_dict = {}
        p200_tipracks = []
        p10_tipracks = []
        # Setup p200
        for location,container in locations.items():
            if container == "tiprack-200ul":
                p200_tipracks.append(containers.load(container, location))
            if container == "tiprack-10ul":
                p10_tipracks.append(containers.load(container, location))
            else:
                if location in annotations:
                    container_dict[annotations[location]] = containers.load(container, location)
                else:
                    container_dict[location] = containers.load(container, location)
        p10 = instruments.Pipette(
            axis='a',
            max_volume=10,
            min_volume=0.5,
            tip_racks=p10_tipracks,
            trash_container=container_dict["Trash"],
            channels=8,
            name='p10-8',
            aspirate_speed=400,
            dispense_speed=800
        )
        p200 = instruments.Pipette(
            axis='b',
            max_volume=200,
            min_volume=20,
            tip_racks=p200_tipracks,
            trash_container=container_dict["Trash"],
            channels=1,
            name='p200-1',
            aspirate_speed=400,
            dispense_speed=800
        )

        output = { 
                'container_dict':container_dict,
                'p10':p10,
                'p200':p200,
                'p200_tipracks':p200_tipracks,
                'p10_tipracks':p10_tipracks                }
        return output
                


def show_opentrons_locations(locations, annotations, OT="OT1"):
    if OT == "OT1":
        locations = {**locations, **annotations}
        # Make the dataframe to represent the OT-1 deck
        deck = ['A1','B2','C3','D2','E1']
        slots = pd.Series(deck)
        columns = sorted(slots.str[0].unique())
        rows = sorted(slots.str[1].unique(), reverse=True)
        layout_table = pd.DataFrame(index=rows, columns=columns)
        layout_table.fillna("---", inplace=True)
        for key,obj in locations.items():
            layout_table.loc[key[1], key[0]] = obj

    if OT == "OT2":
        locations["12"] = "12"
        layout_dict = {
                "1":('Row_4','Col_1'),
                "2":('Row_4','Col_2'),
                "3":('Row_4','Col_3'),
                "4":('Row_3','Col_1'),
                "5":('Row_3','Col_2'),
                "6":('Row_3','Col_3'),
                "7":('Row_2','Col_1'),
                "8":('Row_2','Col_2'),
                "9":('Row_2','Col_3'),
                "10":('Row_1','Col_1'),
                "11":('Row_1','Col_2'),
                "12":('Row_1','Col_3')
                }
        layout_table = pd.DataFrame(index=['Row_1','Row_2','Row_3','Row_4'], columns=['Col_1','Col_2','Col_3']).fillna("---")
        for obj in locations:
            if str(locations[obj]):
                object_to_place = str(locations[obj])
            else:
                object_to_place = str(obj)
            layout_table.loc[layout_dict[obj][0], layout_dict[obj][1]] = object_to_place

    # Displays the required plate map and waits to proceed
    print("\n Please arrange the items in the following configuration: \n")
    print(layout_table,"\n")
    input("Press enter to continue")



