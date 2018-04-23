from matplotlib import pyplot as plt
from scipy.signal import butter, lfilter, freqz, argrelextrema

from imutils import contours
from skimage import measure
import numpy as np
import pandas as pd
import itertools
import argparse
import imutils
import cv2
import glob
import os
import getch
import math

from opentrons import robot, containers, instruments
import sys

def show_image(image):
	cv2.imshow("Image", image)
	# cv2.waitKey(1)
	input()
	return

def increase_contrast(image):
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(3,3))
    cl1 = clahe.apply(image)
    return cl1

def find_colonies(image):

    labels = measure.label(image, neighbors=8, background=0)
    mask = np.zeros(image.shape, dtype="uint8")

    cnts = cv2.findContours(image.copy(), cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)
    cnts = cnts[0] if imutils.is_cv2() else cnts[1]
    #
    color = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
    cv2.drawContours(color, cnts, -1, (255,0,0), 1)
    # show_image(color)

    counter = 0
    for cnt in cnts:
    # for label,cnt in zip(np.unique(labels),cnts):
    # for label in np.unique(labels):
        counter += 1

        # if label == 0:
        #     continue
        #
        # labelMask = np.zeros(image.shape, dtype="uint8")
        # labelMask[labels == label] = 255
        # numPixels = cv2.countNonZero(labelMask)

        M = cv2.moments(cnt)
        if M["m00"] == 0:
            continue
        cX = int(M["m10"] / M["m00"])
        cY = int(M["m01"] / M["m00"])

        perimeter = cv2.arcLength(cnt,True)
        approx = cv2.approxPolyDP(cnt,0.01*perimeter,True)
        area = cv2.contourArea(cnt)

        (rX,rY),(w,h),ang = cv2.minAreaRect(cnt)


        # if area > 10 and area < 300 and len(approx) < 15 and area > perimeter:
            # if w/h > 0.8 and w/h < 1.25:
        if area > 20 and area < 500 and w/h > 0.8 and w/h < 1.25:
            cv2.circle(color, (cX, cY), 1, (0, 0, 255), -1)
            cv2.putText(color, "{}".format(counter), (cX - 5, cY - 5),
            cv2.FONT_HERSHEY_SIMPLEX, 0.45, (0, 255, 255), 1)
            ((cX, cY), radius) = cv2.minEnclosingCircle(cnt)
            cv2.circle(cl1, (int(cX), int(cY)), int(radius+3),
                (255, 0, 0), 1)

            rect = cv2.minAreaRect(cnt)

            print(counter,".","Area",area,"Peri",perimeter,"comp",(perimeter**2)/(4*math.pi),"Approx",len(approx))
            print(w,h,w/h)
            show_image(color)

    # show_image(color)






if __name__ == "__main__":
    for file in sorted(glob.glob('./image_set/*_cropped*')):
        print(file)
        image = cv2.imread(file)
        # resized = imutils.resize(image,width=int(image.shape[1]/3))
        # show_image(resized)

        gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        cl1 = increase_contrast(gray)
        show_image(cl1)

        thresh = cv2.threshold(cl1, 180, 255, cv2.THRESH_BINARY)[1]
        show_image(thresh)

        find_colonies(thresh)







#
