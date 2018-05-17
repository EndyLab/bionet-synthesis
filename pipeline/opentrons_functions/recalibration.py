import os
from opentrons import robot, containers, instruments
import getch

def change_height(container,target,location_data):
    p10 = location_data["p10"]
    x = 0
    counter = 0
    print("Change height - s-g:up h-l:down x:exit")
    while x == 0:
        c = getch.getch()
        if c == "s":
            print("Up 20mm")
            p10.robot._driver.move(z=20,mode="relative")
        elif c == "d":
            print("Up 5mm")
            p10.robot._driver.move(z=5,mode="relative")
        elif c == "f":
            print("Up 0.5mm")
            p10.robot._driver.move(z=0.5,mode="relative")
        elif c == "g":
            print("Up 0.1mm")
            p10.robot._driver.move(z=0.1,mode="relative")
        elif c == "h":
            print("Down 0.1mm")
            p10.robot._driver.move(z=-0.1,mode="relative")
        elif c == "j":
            print("Down 0.5mm")
            p10.robot._driver.move(z=-0.5,mode="relative")
        elif c == "k":
            print("Down 5mm")
            p10.robot._driver.move(z=-5,mode="relative")
        elif c == "l":
            print("Down 20mm")
            p10.robot._driver.move(z=-20,mode="relative")
        elif c == "x":
            print("Exit")
            x = 1
        counter += 1
    if counter > 1:
        print("Will recalibrate")
        redo = True
    else:
        print("Calibrated")
        redo = False
    p10.calibrate_position((container,target.from_center(x=0, y=0, z=-1,reference=container)))
    return redo


def simulation(simulate=True):
    if not simulate:
        port = os.environ["ROBOT_DEV"]
        print("Connecting robot to port {}".format(port))
        robot.connect(port)
    elif simulate:
        print("Simulating protcol run")
        robot.connect()

