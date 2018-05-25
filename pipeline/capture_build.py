import numpy as np
import cv2
import imutils
import glob
import math
import os
from db_config import *
session,engine = connect_db()
from config import *

def capture_plate(path,plate_name):
    '''
    Activates the webcam and streams its feed until you press 'q' on the keyboard
    while on the image. It will then save the frame to the desired path.
    '''
    cam = cv2.VideoCapture(0)
    counter = 0
    while(True):
        ret, frame = cam.read()

        img = np.copy(frame)
        img = imutils.resize(img, width=1000)
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

        cv2.imshow('img', gray)

        if cv2.waitKey(1) & 0xFF == ord('q'):
            name = input('enter build and plate number (e.g. 10p2): ')
            cv2.imwrite('{}/{}.jpg'.format(path,plate_name),frame)
            print("Photo saved, click on the image to move to next plate")
            cv2.waitKey(0)
            break

    cam.release()
    cv2.destroyAllWindows()

def make_directory(assembly):
    '''Checks to see if the directory already exists and if it doesn't it makes a new one'''

    build_num = assembly.builds.build_name
    path = '{}/builds/{}/{}_trans_pics'.format(BASE_PATH,build_num,build_num)
    if os.path.exists(path):
        print("Directory {}_trans_pics already exists".format(build_num))
    else:
        os.makedirs(path)
        print("Making directory for {}_trans_pics".format(build_num))
    return path

def capture_build():
    assemblies = []
    print("Choose which plate you would like to transform/plate:")
    for index,assembly in enumerate(session.query(Plate).join(Build,Plate.builds).filter(Plate.plated == 'not_plated').order_by(Build.build_name)):
        print("{}. {}".format(index,assembly.builds.build_name))
        assemblies.append(assembly)
    plate_num = int(input("Enter plate here: "))
    target_plate = assemblies[plate_num]
    photo_path = make_directory(target_plate)

    build_name = target_plate.builds.build_name
    num_reactions = len(target_plate.wells)
    num_plates = math.ceil(num_reactions/24)
    plate_names = [build_name + '_p' + str(num + 1) for num in range(num_plates)]
    print(plate_names)

    for plate in plate_names:
        print("Place the plate labelled {} inside the incubator".format(plate))
        capture_plate(photo_path,plate)
    print('All of the plates have been captured')

if __name__ == "__main__":
    capture_build()










#
