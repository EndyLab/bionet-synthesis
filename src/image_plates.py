import numpy as np
import cv2
import imutils
import glob
import math
import os
from db_config import *
session,engine = connect_db()
from config import *

import ot_functions as ot

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
            cv2.imwrite('{}/{}.jpg'.format(path,plate_name),frame)
            print("Photo saved, click on the image to move to next plate")
            cv2.waitKey(0)
            break

    cam.release()
    cv2.destroyAllWindows()

def capture_build():
    assemblies = []
    print("Choose which build you would like to photograph:")
    for index,assembly in enumerate(session.query(Plate).join(Build,Plate.builds).filter(Build.status == 'building').order_by(Build.build_name)):
        print("{}. {}".format(index,assembly.builds.build_name))
        assemblies.append(assembly)
    plate_num = int(input("Enter plate here: "))
    target_plate = assemblies[plate_num]

    build_num = assembly.builds.build_name
    photo_path = '{}/builds/{}/{}_trans_pics'.format(BASE_PATH,build_num,build_num)

    ot.make_directory(photo_path)

    build_name = target_plate.builds.build_name
    num_reactions = len(target_plate.wells)
    num_plates = math.ceil(num_reactions/24)
    plate_names = [build_name + '_p' + str(num + 1) for num in range(num_plates)]
    print(plate_names)

    for plate in plate_names:
        print("Place the plate labelled {} inside the image box".format(plate))
        capture_plate(photo_path,plate)
    print('All of the plates have been captured')

if __name__ == "__main__":
    capture_build()










#
