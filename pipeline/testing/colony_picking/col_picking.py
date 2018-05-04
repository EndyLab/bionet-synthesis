from matplotlib import pyplot as plt
from scipy.signal import butter, lfilter, freqz, argrelextrema
# from scipy.optimize import curve_fit
# from pylab import *
import matplotlib.mlab as mlab
from scipy.stats import norm

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
from itertools import chain

from opentrons import robot, containers, instruments
import sys

from cam_calibrate import *

def order_points(pts):
	'''
	Takes in a list of (x,y) coordinates and returns them in clockwise order
	starting at top left
	'''
	# Allocate memory for points
	rect = np.zeros((4, 2), dtype = "float32")

	# The top-left point will have the smallest sum, whereas
	# the bottom-right point will have the largest sum
	s = np.sum(pts,axis=1)
	rect[0] = pts[np.argmin(s)]
	rect[2] = pts[np.argmax(s)]

	# The top-right point will have the smallest difference,
	# whereas the bottom-left will have the largest difference
	diff = np.diff(pts, axis = 1)
	rect[1] = pts[np.argmin(diff)]
	rect[3] = pts[np.argmax(diff)]

	# return the ordered coordinates
	return rect

def four_point_transform(image, pts):
	'''
	Uses four points in order to remap an object in the image to be a perfect rectangle
	'''
	# obtain a consistent order of the points and unpack them
	# individually
	rect = order_points(pts)
	(tl, tr, br, bl) = rect

	# compute the width of the new image, which will be the
	# maximum distance between bottom-right and bottom-left
	# x-coordiates or the top-right and top-left x-coordinates
	widthA = np.sqrt(((br[0] - bl[0]) ** 2) + ((br[1] - bl[1]) ** 2))
	widthB = np.sqrt(((tr[0] - tl[0]) ** 2) + ((tr[1] - tl[1]) ** 2))
	maxWidth = max(int(widthA), int(widthB))

	# compute the height of the new image, which will be the
	# maximum distance between the top-right and bottom-right
	# y-coordinates or the top-left and bottom-left y-coordinates
	heightA = np.sqrt(((tr[0] - br[0]) ** 2) + ((tr[1] - br[1]) ** 2))
	heightB = np.sqrt(((tl[0] - bl[0]) ** 2) + ((tl[1] - bl[1]) ** 2))
	maxHeight = max(int(heightA), int(heightB))

	# Set the desired dimensions based on the max dimensions of the original image
	dst = np.array([
		[0, 0],
		[maxWidth - 1, 0],
		[maxWidth - 1, maxHeight - 1],
		[0, maxHeight - 1]], dtype = "float32")

	# compute the perspective transform matrix and then apply it
	M = cv2.getPerspectiveTransform(rect, dst)
	warped = cv2.warpPerspective(image, M, (maxWidth, maxHeight))

	# return the warped image
	return warped

def find_corners(image):
	'''
	Searches the image to find the edges of the plate and calculates
	the intersection points
	'''
	image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
	blurred = cv2.GaussianBlur(image, (11, 11), 0)

	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(blurred)
	# show_image(cl1)
	thresh = cv2.threshold(cl1, 50, 255, cv2.THRESH_BINARY)[1]
	# show_image(thresh)

	labels = measure.label(thresh, neighbors=8, background=0)
	mask = np.zeros(thresh.shape, dtype="uint8")

	# loop over the unique components
	for label in np.unique(labels):
		# if this is the background label, ignore it
		if label == 0:
			continue

		# otherwise, construct the label mask and count the
		# number of pixels
		labelMask = np.zeros(thresh.shape, dtype="uint8")
		labelMask[labels == label] = 255
		numPixels = cv2.countNonZero(labelMask)

		# if the number of pixels in the component is sufficiently
		# large, then add it to our mask of "large blobs"
		if numPixels > 100000:
			mask = cv2.add(mask, labelMask)
	# show_image(mask)

	cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_NONE)
	cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	maskc = np.zeros(thresh.shape, dtype="uint8")
	cv2.drawContours(maskc, cnts, -1, (255,0,0), 1)

	epsilon = 0.1*cv2.arcLength(cnts[0],True)
	pts = cv2.approxPolyDP(cnts[0],epsilon,True)
	image = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
	cv2.polylines(image,[pts],True,(0,255,255),1)

	corners = [[x,y] for [[x,y]] in pts]
	corners = np.array(corners)
	print("Corners:",corners)

	return image, corners

def find_grid(img):
	'''
	Searches the image for rows and columns with no colonies to define a grid
	segmenting the plate based on part
	'''
	sums_y = []
	sums_x = []

	# create a CLAHE object (Arguments are optional).
	resized = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(resized)
	thresh = cv2.threshold(cl1, 170, 255, cv2.THRESH_BINARY)[1]
	show_image(thresh)
	height,width = thresh.shape[:2]
	for row in range(height):
		sum = np.sum(thresh[row,0:width])
		sums_y.append(sum)

	for col in range(width):
		sum = np.sum(thresh[0:height,col])
		sums_x.append(sum)

	# Generate a plot to visualize the sumations
	fig1 = plt.figure(0)
	ax1 = fig1.add_subplot(211)
	ax1.plot(sums_y)
	# Store only the points that have less than 20 pixels in the line
	bottoms = []
	for i,sum in enumerate(sums_y):
		if sum/255 < 40:
			ax1.plot(i,sum, 'r.')
			bottoms.append(i)

	# Finds the start and end points of streaks of low points
	edges = [bottoms[0]]
	for i,bot in enumerate(bottoms):
		if i == (len(bottoms) - 1):
			break
		if bot != (bottoms[i+1] - 1) or bot != (bottoms[i-1] + 1):
			if bottoms[i+1] - bot > 20 or bot - bottoms[i-1] > 20:
				edges.append(bot)
	edges.append(bottoms[-1])
	# Groups the endpoints together
	bounds = [edges[n:n+2] for n in range(0, len(edges), 2)]
	print("Bounds: ",bounds)
	temp = np.copy(cl1)

	# Calculates and plots the midpoints of the low regions
	mid_y = []
	for l,r in bounds:
		mid_point = int(((r-l) / 2)+l)
		mid_y.append(int(mid_point))#+top_bound/2))
		ax1.plot(mid_point,sums_y[mid_point], 'bo')
		# cv2.line(temp,(0,int(mid_point+top_bound/2)),(int(cl1.shape[1]),int(mid_point+top_bound/2)),(0,255,0),1)
		cv2.line(temp,(0,int(mid_point)),(int(height),int(mid_point)),(0,255,0),1)
	print("Number of rows: ",len(mid_y)-1)

	# Generate a plot to visualize the sumations
	ax2 = fig1.add_subplot(212)
	ax2.plot(sums_x)
	# ax2.ylabel('Row sum')

	# Store only the points that have less than 60 pixels in the line
	bottoms = []
	for i,sum in enumerate(sums_x):
		if sum/255 < 60:
			ax2.plot(i,sum, 'r.')
			bottoms.append(i)
	# Finds the start and end points of streaks of low points
	edges = [bottoms[0]]
	for i,bot in enumerate(bottoms):
		if i == (len(bottoms) - 1):
			break
		if bot != (bottoms[i+1] - 1) or bot != (bottoms[i-1] + 1):
			if bottoms[i+1] - bot > 20 or bot - bottoms[i-1] > 20:
				edges.append(bot)
	edges.append(bottoms[-1])

	# Groups the endpoints together
	bounds = [edges[n:n+2] for n in range(0, len(edges), 2)]
	print("Bounds: ",bounds)
	# Calculates and plots the midpoints of the low regions
	mid_x = []

	for t,b in bounds:
		mid_point = int(((b-t) / 2)+t)
		mid_x.append(mid_point)
		ax2.plot(mid_point,sums_x[mid_point], 'bo')
		cv2.line(temp,(int(mid_point),0),(int(mid_point),int(height)),(255,0,0),1)

	print("Number of columns: ",len(mid_x)-1)
	show_image(temp)

	fig1.show()
	input('Review grid')
	plt.close()

	return cl1, mid_x, mid_y

def group_sections(img,mid_x,mid_y,x_steps=1,y_steps=4):
	groups = []
	for y in range(0,len(mid_y)-1,y_steps):
		for x in range(0,len(mid_x)-1,x_steps):
			group = [[mid_x[x],mid_x[x+x_steps]],[mid_y[y],mid_y[y+y_steps]]]
			print(group)
			groups.append(group)
			cv2.rectangle(img,(mid_x[x],mid_y[y]),(mid_x[x+x_steps],mid_y[y+y_steps]),(0,0,255),3)
			cv2.imshow("Image", img)
			cv2.waitKey(10)
	return groups

def show_small(img):
	resized = imutils.resize(img,width=100)
	cv2.imshow("Image", resized)
	cv2.waitKey(1)
	input("next image")

def scale_linear_bycolumn(rawpoints, high=100.0, low=0.0):
    mins = np.min(rawpoints, axis=0)
    maxs = np.max(rawpoints, axis=0)
    rng = maxs - mins
    return high - (((high - low) * (maxs - rawpoints)) / rng)

def find_colonies(image):
	'''
	Searches the image for pickable colonies and returns the
	location of their centers
	'''
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(image)
	show_small(cl1)

	# nested = cl1.tolist()
	# unnest = list(chain(*nested))
	# acceptable = [x for x in unnest if x > 130]
	# counts = np.bincount(np.array(acceptable))
	# print(counts)
	# print(len(counts))
	# factor = 255/len(counts)
	# print(factor)
	# bound = np.argmax(counts[:170])
	# print(bound,counts[bound])
	# low = np.argmin(counts[bound:int(factor*200)])
	# print("pre",bound+low)
	# lowest = int((bound+low)*factor)
	# print("post",lowest)
	# print(counts[bound+low])
	# input('check counts')


	fig2 = plt.figure(1)
	ax2 = fig2.add_subplot(111)
	nums, bins, patches = ax2.hist(cl1.flatten(),bins=25)
	data = zip(bins[1:],nums)
	# print("nums",nums)
	# print("bins",bins)

	bound = np.argmax(nums)
	# print('Bound',bound)
	low = np.argmin(nums[bound:-3])
	# print("low",low+bound)
	threshold = bins[low+bound+2]
	if threshold >= 200:
		threshold = 190
	# print('thres',thr)

	# bottoms = []
	# for i,(b,n) in enumerate(data):
	# 	if n < 1000:
	# 		plt.plot(b,n, 'b.')
	# 		bottoms.append(i)
	# ends = []
	# print("bottoms",bottoms)
	# count = 0
	# for i in bottoms:
	# 	count += 1
	# 	if count == (len(bottoms)):
	# 		break
	# 	if i != bottoms[count]-1:
	# 		ends.append(i)
	# for e in ends:
	# 	ax2.plot(bins[:-1][e],nums[e], 'r.')
	ax2.plot(threshold*np.ones(5000), range(5000), marker='o', markersize=3, color="green")
	fig2.show()
	input('Check hist')
	plt.close()
	# threshold = bins[:-1][ends[-1]]
	# if threshold < 160:
	# 	threshold = 170
	# threshold = lowest
	print("thresh: ",threshold)
	thresh = cv2.threshold(cl1, threshold, 255, cv2.THRESH_BINARY)[1]
	show_small(thresh)
	cnts = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_NONE)
	cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	color = cv2.cvtColor(cl1,cv2.COLOR_GRAY2RGB)
	color_th = cv2.cvtColor(thresh,cv2.COLOR_GRAY2RGB)

	centers = []
	good_cols = []
	bad_cols = []
	good_center = []
	for counter,cnt in enumerate(cnts):
		print("counter: ",counter)
		M = cv2.moments(cnt)
		if M["m00"] == 0:
			continue
		cX = int(M["m10"] / M["m00"])
		cY = int(M["m01"] / M["m00"])

		perimeter = cv2.arcLength(cnt,True)
		approx = cv2.approxPolyDP(cnt,0.01*perimeter,True)

		perimeter = cv2.arcLength(approx,True)
		area = cv2.contourArea(approx)
		calc = area/((perimeter**2)/(4*math.pi))

		(rX,rY),(w,h),ang = cv2.minAreaRect(approx)
		ratio = w/h

		attributes = [area,perimeter,calc,ratio]

		# Set the acceptable ranges for the different attributes
		area_r = [10,375]
		peri_r = [12.5,75]
		calc_r = [0.6,1]
		ratio_r = [0.7,1.5]

		crit = [area_r,peri_r,calc_r,ratio_r]

		# Check each attribute against the ranges specified
		bad = False
		if perimeter > area:
			bad = True
		for att,(min,max) in zip(attributes,crit):
			print(att,min,max)
			if att < min or att > max:
				bad = True
		if bad:
			print('Bad')
			bad_cols.append(cnt)
			# cv2.drawContours(color, cnt, -1, (0,0,255), 1)
			# cv2.drawContours(color_th, cnt, -1, (0,0,255), 1)
			# show_small(color_th)
			# show_small(color)
		else:
			print('Good')
			good_cols.append(cnt)
			centers.append([[cX,cY]])
			# cv2.drawContours(color, cnt, -1, (0,255,0), 1)
			# cv2.drawContours(color_th, cnt, -1, (0,255,0), 1)
			# show_small(color_th)
			# show_small(color)
	if len(centers) == 0:
		print("No colonies found")
	else:
		good_center = centers[-1]
	cv2.drawContours(color, good_cols, -1, (0,255,0), 1)
	cv2.drawContours(color_th, good_cols, -1, (0,255,0), 1)
	cv2.drawContours(color, bad_cols, -1, (0,0,255), 1)
	cv2.drawContours(color_th, bad_cols, -1, (0,0,255), 1)
	show_small(color_th)
	show_small(color)

	# 	attributes = [area,perimeter,area/((perimeter**2)/(4*math.pi)),ratio]
	# 	data.append(attributes)
	# 	print("Approx",counter,".","Area",area,"Peri",perimeter,"comp",area/((perimeter**2)/(4*math.pi)),"Approx",len(approx),"Ratio",ratio)
	# list(data).append(targets)
	# normed = scale_linear_bycolumn(data)
	# normed_avg = normed[-1]
	# normed_scores = normed[:-1]
	# diffs = abs(normed_scores - normed_avg)
	# scores = enumerate(list(np.sum(diffs,axis=1)))
	# print(scores)
	# order = sorted(scores,key=lambda tup:tup[1])
	# print(order)
	#
	# for i,score in order:
	# 	cnts[i]
	# 	temp_c = np.copy(color)
	# 	temp_i = np.copy(image)
	# 	cv2.drawContours(temp_c, cnts[i], -1, (0,0,255), 1)
	# 	cv2.drawContours(temp_i, cnts[i], -1, (0,0,255), 1)
	# 	print(data[i])
	# 	print(score)
	# 	show_small(temp_c)
	# 	# show_small(temp_i)
	# 	input("next section")



		# temp_c = np.copy(color)
		# temp_i = np.copy(image)
		# cv2.circle(color_th, (cX, cY), 1, (0, 0, 255), -5)
		# cv2.circle(color, (cX, cY), 1, (0, 0, 255), -5)
		#
		# cv2.drawContours(color_th, [approx], -1, (0, 0, 255), 1)
		# cv2.drawContours(color, [approx], -1, (0, 0, 255), 1)
		# #
		# # show_small(color_th)
		# show_small(color)
		# col = input("g or b?: ")
		# if col == 'g':
		# 	print('good')
		# 	good_cols.append(attributes)
		# elif col == 'b':
		# 	print('bad')
		# 	bad_cols.append(attributes)



		# if area > 50 and area < 350 and w/h > 0.8 and w/h < 1.25 and area > perimeter:
		# if area > 40 and area < 200 and w/h > 0.8 and w/h < 1.25 and area > perimeter:
		# 	temp_c = np.copy(color)
		# 	temp_i = np.copy(image)
		# 	cv2.circle(temp_c, (cX, cY), 1, (0, 0, 255), -3)
		# 	cv2.circle(temp_i, (cX, cY), 1, (0, 0, 255), -3)
		# 	# cv2.putText(color, "{}".format(counter), (cX - 5, cY - 5),
		# 	# 	cv2.FONT_HERSHEY_SIMPLEX, 0.45, (0, 255, 255), 3)
		# 	((cX, cY), radius) = cv2.minEnclosingCircle(cnt)
		# 	cv2.circle(temp_c, (int(cX), int(cY)), int(radius+9),
		# 		(0, 255, 255), 2)
		# 	cv2.circle(temp_i, (int(cX), int(cY)), int(radius+9),
		# 		(0, 255, 255), 2)

			# print(counter,".","Area",area,"Peri",perimeter,"comp",(perimeter**2)/(4*math.pi),"Approx",len(approx))
			# print(w,h,w/h)
			# show_image(temp_c)
			# # show_image(temp_i)
			# print(cX,cY)
	# print("GOOD:",good_cols)
	# print("BAD:",bad_cols)
	return good_center

def change_loc(image,start):
	'''
	Allows for manually modification of points presented on images
	'''
	x = 0
	temp = np.copy(image)
	# cv2.circle(temp, (int(start[0]),int(start[1])), int(30),
	# 	(255,0,0), 9)
	cv2.line(temp,(int(start[0]),int(start[1])-40),(int(start[0]),int(start[1])+40),(0,0,255),2)
	cv2.line(temp,(int(start[0])-40,int(start[1])),(int(start[0])+40,int(start[1])),(0,0,255),2)
	temp_resized = imutils.resize(temp,width=400)
	cv2.imshow("Image", temp_resized)
	cv2.waitKey(5)
	print("please fix the location")
	while x == 0:
		# temp = image
		c = getch.getch()
		if c == "w":
			start[1] -= 1
			print("up")
		elif c == "s":
			start[1] += 1
			print("down")
		elif c == "a":
			start[0] -= 1
			print("left")
		elif c == "d":
			start[0] += 1
			print("right")
		elif c == "x":
			x = 1
			print("exit")
		temp = np.copy(image)
		# cv2.circle(temp, (int(start[0]),int(start[1])), int(30),
		# 	(255,0,0), 9)
		cv2.line(temp,(int(start[0]),int(start[1])-40),(int(start[0]),int(start[1])+40),(0,0,255),2)
		cv2.line(temp,(int(start[0])-40,int(start[1])),(int(start[0])+40,int(start[1])),(0,0,255),2)
		temp_resized = imutils.resize(temp,width=400)
		cv2.imshow("Image", temp_resized)
		cv2.waitKey(5)

	return start


def find_reference(image):
	'''
	Locate the corners of the agar to link the reference point on the robot
	to the image
	'''
	# image = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)

	all_points = [[19,27],[377,27],[380,612],[19, 614]]
	all_points = np.array(all_points)
	ratio = image.shape[1]/400
	all_points = np.multiply(all_points,ratio)
	rect = order_points(all_points)
	new_points = []
	for point in rect:
		new_points.append(change_loc(image,point))
		print("Started at: ",point)
		print("Ended at: ",new_points[-1],"\n")

	ref_point = new_points[3]
	print(new_points)
	print(ref_point)
	x_dim = new_points[1][0] - new_points[3][0]
	y_dim = new_points[3][1] - new_points[1][1]

	return ref_point,x_dim,y_dim,rect

def ot_coords(centers,image):
	'''
	Establish a relationship between pixels in the image to xy coords on the robot
	'''
	y_max = 121
	x_max = 78
	coords = []
	for cen in centers:
		x = float(((cen[0]/x_dim) * x_max))
		y = float(((image.shape[0]-cen[1])/y_dim) * y_max)
		if x < float(0) or y < float(0):
			print("invalid coordinate: ",x,y)
			continue
		coords.append([x,y])
	return coords

def move_motor():
	'''Allows for real-time calibration of the p10single channel xy location'''
	z = 0
	print("Change xy - w:back, d:forward, a:left, d:right, x:exit")
	x = 0
	y = 0
	while z == 0:
		c = getch.getch()
		if c == "w":
			p10s.robot._driver.move(y=0.5,mode="relative")
			print("up")
			y += 0.5
		elif c == "s":
			p10s.robot._driver.move(y=-0.5,mode="relative")
			print("down")
			y -= 0.5
		elif c == "a":
			p10s.robot._driver.move(x=-0.5,mode="relative")
			print("left")
			x -= 0.5
		elif c == "d":
			p10s.robot._driver.move(x=0.5,mode="relative")
			print("right")
			x += 0.5
		elif c == "x":
			z = 1
	print(x,y)
	return x,y

def calibrate_ot(image,coords,centers):
	first,_,last,_ = order_points(coords)
	print(last,first)
	input("check")
	_,last_cen,_,first_cen = order_points(centers)

	temp = np.copy(image)

	cv2.line(temp,(int(first_cen[0]),int(first_cen[1])-30),(int(first_cen[0]),int(first_cen[1])+30),(0,0,255),2)
	cv2.line(temp,(int(first_cen[0])-30,int(first_cen[1])),(int(first_cen[0])+30,int(first_cen[1])),(0,0,255),2)
	show_image(temp)
	p10s.move_to((trans_plate,[first[0],first[1],0]))
	print("Calibrate the first colony")
	off_x1,off_y1 = move_motor()

	cv2.line(temp,(int(last_cen[0]),int(last_cen[1])-30),(int(last_cen[0]),int(last_cen[1])+30),(0,0,255),2)
	cv2.line(temp,(int(last_cen[0])-30,int(last_cen[1])),(int(last_cen[0])+30,int(last_cen[1])),(0,0,255),2)
	show_image(temp)
	p10s.move_to((trans_plate,[last[0],last[1],0]))
	print("Calibrate the last colony")
	off_x2,off_y2 = move_motor()

	mX = (off_x2-off_x1)/(last[0]-first[0])
	bX = off_x1 - mX*first[0]
	mY = (off_y2-off_y1)/(last[1]-first[1])
	bY = off_y1 - mY*first[1]

	return mX,bX,mY,bY

def run_ot(image,coords,centers):
	'''
	Pass the coordinates to the robot to pick the colony
	'''
	mX,bX,mY,bY = calibrate_ot(image,coords,centers)
	print(mX,bX,mY,bY)

	for i,((x,y),cen) in enumerate(zip(coords,centers)):
		temp = image
		print(i,":",x,y)
		x_off = (mX * x) + bX
		y_off = (mY * y) + bY
		new_x = x + x_off
		new_y = y + y_off
		print(i,":",new_x,new_y)
		cv2.circle(temp, (int(cen[0]),int(cen[1])), int(20),
			(0,0,255), 2)
		show_image(temp)
		p10s.move_to((trans_plate,[new_x,new_y,0]))
		input("click to move to next colony")

#
def show_image(image,width=400):
	'''Scales and then presents the image'''
	resized = imutils.resize(image,width=width)
	cv2.imshow("Image", resized)
	cv2.waitKey(5)
	input()
	return


grids = []
shapes = []

parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
args = parser.parse_args()

if args.run:
	port = os.environ["ROBOT_DEV"]
	print("Connecting robot to port {}".format(port))
	robot.connect(port)
else:
	print("Simulating protcol run")
	robot.connect()

robot.home()

p10s_tipracks = [
    containers.load('tiprack-10ul', 'E1'),
    containers.load('tiprack-10ul', 'E2')
]
trash = containers.load('point', 'D1', 'holywastedplasticbatman')
well_plate = containers.load('96-deep-well', 'D2')
trans_plate = containers.load('point','B2')


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

# for file in sorted(glob.glob('./image_set/*.jpg')):
for file in sorted(glob.glob('./cam_photos/*.jpg')):
	print(file)
	img = cv2.imread(file)
	show_image(img)
	cali_file = './webcam_calibrations.npz'
	dst = undistort_img(img,cali_file)
	show_image(dst)
	height,width = img.shape[:2]
	rotated = imutils.rotate_bound(dst, 90)
	show_image(rotated)

	edges, intersections = find_corners(rotated)
	show_image(edges)
	warped = four_point_transform(edges,intersections)
	show_image(warped)

	ref, x_dim, y_dim,rect = find_reference(warped)
	agar = four_point_transform(warped,rect)
	show_image(agar)
	grid, mid_x, mid_y = find_grid(agar)
	show_image(grid)
	color = cv2.cvtColor(grid,cv2.COLOR_GRAY2RGB)
	groups = group_sections(grid,mid_x,mid_y)
	centers = []
	for group in groups:
		print("Next section:",group)
		partial = grid[group[1][0]:group[1][1],group[0][0]:group[0][1]]
		center = find_colonies(partial)
		if center == []:
			midg_y = int(((group[1][0]-group[0][0])/2)+group[0][0])
			midg_x = int(((group[1][1]-group[0][1])/2)+group[0][1])
			cv2.rectangle(color,(group[0][0],group[1][0]),(group[0][1],group[1][1]),(0,0,255),3)
			show_image(color)
			input("check image")
		else:
			cali_cenX = center[0][0]+group[0][0]
			cali_cenY = center[0][1]+group[1][0]
			cv2.circle(color, (cali_cenX,cali_cenY), int(20),
				(0,255,255), 3)
			cv2.circle(color, (cali_cenX,cali_cenY), int(2),
				(0,0,255), 2)
			show_image(color)
			centers.append([cali_cenX,cali_cenY])
			input("check image w/ col")
	show_image(color)
	input("check image w/ col")
	coords = ot_coords(centers,agar)
	run_ot(color,coords,centers)
	print("Complete")

# print(shapes)colonies, centers = find_colonies(grid)
# input("stop")
# res = np.hstack(grids)
# cv2.imshow("Image", res)
# cv2.imshow("Image", warped)
# cv2.waitKey(0)
cv2.destroyAllWindows()

## Finding the reference point
	# blurred = cv2.GaussianBlur(image, (15, 15), 0)
	# clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	# cl1 = clahe.apply(blurred)
	# show_image(cl1)
	# thresh = cv2.threshold(cl1, 140, 255, cv2.THRESH_BINARY)[1]
	# show_image(thresh)
	# labels = measure.label(thresh, neighbors=8, background=0)
	# mask = np.zeros(thresh.shape, dtype="uint8")
	#
	# # loop over the unique components
	# for label in np.unique(labels):
	# 	# if this is the background label, ignore it
	# 	if label == 0:
	# 		continue
	#
	# 	# otherwise, construct the label mask and count the
	# 	# number of pixels
	# 	labelMask = np.zeros(thresh.shape, dtype="uint8")
	# 	labelMask[labels == label] = 255
	# 	numPixels = cv2.countNonZero(labelMask)
	#
	# 	# if the number of pixels in the component is sufficiently
	# 	# large, then add it to our mask of "large blobs"
	# 	if numPixels > 5000:
	# 		mask = cv2.add(mask, labelMask)
	# cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
	# 	cv2.CHAIN_APPROX_SIMPLE)
	# cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	# edges = cv2.Canny(mask,50,150,apertureSize = 3)


## Finding hough lines and intersection points:
	# image = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
	# equ = []
	# def find_lines(min):
	# 	'''Finds the edges of the plate'''
	# 	lines = cv2.HoughLines(maskc,1,np.pi/180,min)
	# 	if min < 0:
	# 		print("reached the minimum threshold")
	# 		return lines
	# 	if type(lines) == type(None):
	# 		print("no lines were found")
	# 		return find_lines(min-5)
	# 	# Requires 4 edges and will continue to check until it finds them
	# 	if len(lines) < 4:
	# 		return find_lines(min-5)
	# 	elif len(lines) > 4:
	# 		return find_lines(min+5)
	# 	elif len(lines) == 4:
	# 		return lines
	#
	# lines = find_lines(300)
	# if len(lines) < 4:
	# 	input("Still lacking enough lines")
	#
	# for line in lines:
	# 	for rho,theta in line:
	# 		a = np.cos(theta)
	# 		b = np.sin(theta)
	# 		x0 = a*rho
	# 		y0 = b*rho
	# 		equ += [[a,b,rho]]
	# 		print((a/b),(rho/b))
	# 		x1 = int(x0 + 5000*(-b))
	# 		y1 = int(y0 + 5000*(a))
	# 		x2 = int(x0 - 5000*(-b))
	# 		y2 = int(y0 - 5000*(a))
	# 		cv2.line(image,(x1,y1),(x2,y2),(255,0,0),2)
	#
	# pts = []
	# for pair in itertools.combinations(equ,2):
	# 	print(pair)
	# 	pair = np.array(pair)
	# 	a = pair[0:2,0:2]
	# 	b = pair[0:2,2]
	# 	if round(pair[0,0],2)*round(pair[1,1],2) == round(pair[1,0],2)*round(pair[0,1],2):
	# 		print("lines are parallel")
	# 		continue
	# 	x = np.linalg.solve(a, b)
	# 	pts += [x]
	# 	cv2.circle(image, (int(x[0]), int(x[1])), int(10),
	# 			(255, 0, 0), 2)
	# print("number of intersections: ",len(pts))
	# pts = np.array(pts)
	# print("These are the intersection points:\n",pts)




## Finding rolling averages
	# plt.subplot(2, 1, 2)
	# print("Y dim",thresh.shape[0])
	# print("X dim",thresh.shape[1])
	# print("shape-bounds",thresh.shape[0] - (top_bound*2))
	# print("Len of sums",len(sums))
	# rolling = pd.Series(sums).rolling(window=10).mean()
	# print("Len of rolling",len(rolling))
	# # input("wait")
	# np_plot = np.array(rolling)
	# minima = argrelextrema(np_plot, np.less)
	# print("minima: ",*minima)
	# plt.plot(rolling)
	# plt.ylabel('Rolling average')
	# diffs = []
	# for i,min in enumerate(minima[0]):
	# 	if i == (len(minima[0])-2):
	# 		break
	# 	diffs.append(minima[0][i+1] - min)
	#
	# diff_avg = np.mean(diffs)
	# new_mins = [minima[0][0]]
	# for i,diff in enumerate(diffs):
	# 	if i == (len(diffs)-1):
	# 		break
	# 	print(i,diff)
	# 	print(diff, diffs[i+1],(diff_avg * 0.7))
	# 	# print(diffs[i+1] + diff,(2*diff_avg * 0.7))
	# 	if diff < diff_avg * 0.7 and diffs[i+1] < (diff_avg * 0.7):
	#
	#
	# 	# if (diffs[i+1] + diff) < (2*diff_avg * 0.7):
	# 		print("will exclude: ",diff, diffs[i+1],(diff_avg * 0.7))
	# 	else:
	# 		new_mins.append(minima[0][i+1])
	# 	# summs.append(diff[i+1] + diff)
	# new_mins.append(minima[0][-2])
	# new_mins.append(minima[0][-1])
	#
	# if new_mins[1] - new_mins[0] < diff_avg * 0.7:
	# 	print("Far left is too close")
	# 	new_mins = new_mins[1:]
	# rows = []
	# for row in marks:
	# 	cv2.line(thresh,(left_bound,row),(right_bound,row),(255,0,0),2)
	# 	cv2.putText(thresh, str(row), (left_bound, int(row)+10), cv2.FONT_HERSHEY_SIMPLEX,0.5, (255, 0, 0), 2)
	# for min in new_mins:
	# 	plt.plot(min,rolling[min], 'r.')
	# 	rows.append(min)
	# 	y = int(rolling[int(min)])
	# 	# print(y)
	# 	print("min+top",min+top_bound)
		# cv2.line(thresh,(0,min+top_bound),(int(thresh.shape[1]),min+top_bound),(255,0,0),1)


	# diffs = [diff_avg] + diffs
	# diffs.append(diff_avg)
	# print(minima)
	# print(diffs)
	# print(diff_avg)
	# rows = []
	# for min,diff in zip(new_mins,diffs):
	# 	print(min,diff,diff_avg/1.5)
	# 	if rolling[int(min)] > rolling[int(min)+10] and rolling[int(min)] > rolling[int(min)-10]:
	# 		print("not a minimum")
	# 		continue
	# 	# if diff < diff_avg/1.5:
	# 	# 	print("too close")
	# 	# 	continue
	# 	plt.plot(min,rolling[min], 'r.')
	# 	rows.append(min)
	# 	y = int(rolling[int(min)])
	# 	# print(y)
	# 	cv2.line(thresh,(0,min+top_bound),(int(thresh.shape[1]),min+top_bound),(255,0,0),1)



## Ordering of the intersection points
	# pts = pts[pts[:,0].argsort()]
	# print(pts)
	# if pts[0,1] < pts[1,1]:
	# 	inter.update({'topL' : {'x' : pts[0,0],'y' : pts[0,1]}})
	# 	inter.update({'botL' : {'x' : pts[1,0],'y' : pts[1,1]}})
	# else:
	# 	inter.update({'botL' : {'x' : pts[0,0],'y' : pts[0,1]}})
	# 	inter.update({'topL' : {'x' : pts[1,0],'y' : pts[1,1]}})
	# if pts[2,1] < pts[3,1]:
	# 	inter.update({'topR' : {'x' : pts[2,0],'y' : pts[2,1]}})
	# 	inter.update({'botR' : {'x' : pts[3,0],'y' : pts[3,1]}})
	# else:
	# 	inter.update({'botR' : {'x' : pts[2,0],'y' : pts[2,1]}})
	# 	inter.update({'topR' : {'x' : pts[3,0],'y' : pts[3,1]}})
	# print(inter)
	# print(inter['topR']['x'])

# def make_grid(image,p):
	# Find the change in pixels in both directions across all of the sides of the plate
	# top_lenX = abs(p['topR']['x'] - p['topL']['x'])
	# top_lenY = abs(p['topR']['y'] - p['topL']['y'])
	# offsetT = top_lenX / 10
	# top_lenX -= offsetT * 2

	# left_lenX = abs(p['botL']['x'] - p['topL']['x'])
	# left_lenY = abs(p['botL']['y'] - p['topL']['y'])
	# offsetL = left_lenY / 12
	# left_lenY -= offsetL * 2

	# right_lenX = abs(p['topR']['x'] - p['botR']['x'])
	# right_lenY = abs(p['topR']['y'] - p['botR']['y'])
	# offsetR = right_lenY / 12
	# right_lenY -= offsetR * 2

	# bot_lenX = abs(p['botL']['x'] - p['botR']['x'])
	# bot_lenY = abs(p['botL']['y'] - p['botR']['y'])
	# offsetB = bot_lenX / 10
	# bot_lenX -= offsetB * 2

	# y1 = int(offsetL+p['topL']['y'])
	# x1 = int(p['topL']['x'])
	# y2 = int(offsetR+p['topR']['y'])
	# x2 = int(p['topR']['x'])
	# cv2.line(image, (x1, y1),(x2, y2),(255,0,0),1)
	#
	# num_rows = 12
	# for row in range(num_rows):
	# 	y1 = (left_lenY - left_lenY*(row/num_rows)) + p['topL']['y'] + offsetL
	# 	x1 = (left_lenX - left_lenX*(row/num_rows)) + p['topL']['x']
	# 	y2 = (right_lenY - right_lenY*(row/num_rows)) + p['topR']['y'] + offsetR
	# 	x2 = (right_lenX - right_lenX*(row/num_rows)) + p['topR']['x']
	# 	print(x1,y1,x2,y2)
	#
	# 	cv2.line(image, (int(x1), int(y1)),(int(x2), int(y2)),(255,0,0),1)
	# num_cols = 8
	# cv2.line(image, (int(offsetT+p['topL']['x']), int(p['topL']['y'])),(int(offsetB+p['botL']['x']), int(p['botL']['y'])),(0,0,255),1)
	# for row in range(num_cols):
	# 	y1 = (bot_lenY - bot_lenY*(row/num_cols)) + p['topL']['y']
	# 	x1 = (bot_lenX - bot_lenX*(row/num_cols)) + p['topL']['x'] + offsetB
	# 	y2 = (top_lenY - top_lenY*(row/num_cols)) + p['botL']['y']
	# 	x2 = (top_lenX - top_lenX*(row/num_cols)) + p['botL']['x'] + offsetT
	# 	print(x1,y1,x2,y2)
	#
	# 	cv2.line(image, (int(x1), int(y1)),(int(x2), int(y2)),(0,0,255),1)
	# return image

# # create a CLAHE object (Arguments are optional).
# resized = cv2.cvtColor(resized, cv2.COLOR_BGR2GRAY)
# clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
# cl1 = clahe.apply(resized)
# thresh = cv2.threshold(cl1, 180, 255, cv2.THRESH_BINARY)[1]

# # perform a connected component analysis on the thresholded
# # image, then initialize a mask to store only the "large"
# # components
# labels = measure.label(thresh, neighbors=8, background=0)
# mask = np.zeros(thresh.shape, dtype="uint8")
#
# # loop over the unique components
# for label in np.unique(labels):
# 	# if this is the background label, ignore it
# 	if label == 0:
# 		continue
#
# 	# otherwise, construct the label mask and count the
# 	# number of pixels
# 	labelMask = np.zeros(thresh.shape, dtype="uint8")
# 	labelMask[labels == label] = 255
# 	numPixels = cv2.countNonZero(labelMask)
#
# 	# if the number of pixels in the component is sufficiently
# 	# large, then add it to our mask of "large blobs"
# 	# if numPixels > 50:
# 	if numPixels > 30 and numPixels < 100:
# 	# if numPixels < 50 and numPixels > 5:
# 		mask = cv2.add(mask, labelMask)
#
# cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
# 	cv2.CHAIN_APPROX_NONE)
# cnts = cnts[0] if imutils.is_cv2() else cnts[1]
#
# # edges = cv2.Canny(mask,50,150,apertureSize = 3)
# cl1 = cv2.cvtColor(cl1,cv2.COLOR_GRAY2RGB)
#
# # cv2.drawContours(cl1, cnts, -1, (255,0,0), 1)
#
# for c in cnts:
# 	# compute the center of the contour
# 	M = cv2.moments(c)
# 	if M["m00"] == 0:
# 		continue
# 	cX = int(M["m10"] / M["m00"])
# 	cY = int(M["m01"] / M["m00"])
# 	cv2.circle(cl1, (cX, cY), 1, (0, 0, 255), -1)
# 	((cX, cY), radius) = cv2.minEnclosingCircle(c)
# 	cv2.circle(cl1, (int(cX), int(cY)), int(radius+3),
# 		(255, 0, 0), 1)

# res = np.hstack((cl1,thresh,mask))
# cv2.imshow("Image", colonies)
# cv2.waitKey(0)
# cv2.destroyAllWindows()

# cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
# 	cv2.CHAIN_APPROX_SIMPLE)
# cnts = cnts[0] if imutils.is_cv2() else cnts[1]
# e = []
# mask = cv2.cvtColor(mask,cv2.COLOR_GRAY2RGB)
# cl1 = cv2.cvtColor(cl1,cv2.COLOR_GRAY2RGB)
# cv2.drawContours(cl1, cnts, -1, (0,255,0), 3)
#
#
# for c in cnts:
# 	M = cv2.moments(c)
# 	cX = int((M["m10"] / M["m00"]))
# 	cY = int((M["m01"] / M["m00"]))
# 	peri = cv2.arcLength(c, True)
# 	approx = cv2.approxPolyDP(c, 0.04 * peri, True)
# 	print(peri)
# 	print(approx)
# 	print(len(approx))
#
# 	c = c.astype("float")
# 	c *= 1.2
# 	c = c.astype("int")
# 	i = 0
# 	for set in approx:
# 		for (x,y) in set:
# 			print(x,y)
# 			e += [[x,y]]
# 			cv2.circle(mask, (int(x), int(y)), int(5),
# 					(255, 0, 0), 2)
# 			cv2.putText(mask, "#{}".format(i + 1), (x - 20, y),
# 				cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 0, 0), 2)
# 			i += 1
# 	print(e)
# 	e = np.array(e)
# 	leny1 = e[1,1] - e[0,1]
# 	lenx1 = e[1,0] - e[0,0]
# 	leny2 = e[2,1] - e[3,1]
# 	lenx2 = e[2,0] - e[3,0]
#
# 	leny3 = e[3,1] - e[0,1]
# 	lenx3 = e[3,0] - e[0,0]
# 	leny4 = e[2,1] - e[1,1]
# 	lenx4 = e[2,0] - e[1,0]
#
#
# 	for row in range(14):
# 		y1 = (leny1 - leny1*(row/14)) + e[0,1]
# 		x1 = (lenx1 - lenx1*(row/14)) + e[0,0]
# 		y2 = (leny2 - leny2*(row/14)) + e[3,1]
# 		x2 = (lenx2 - lenx2*(row/14)) + e[3,0]
# 		print(x1,y1,x2,y2)
#
# 		cv2.line(cl1, (int(x1), int(y1)),(int(x2), int(y2)),(0,0,255),1)
#
# 	for col in range(10):
# 		y3 = (leny3 - leny3*(col/10)) + e[0,1]
# 		x3 = (lenx3 - lenx3*(col/10)) + e[0,0]
# 		y4 = (leny4 - leny4*(col/10)) + e[1,1]
# 		x4 = (lenx4 - lenx4*(col/10)) + e[1,0]
# 		print(x3,y3,x4,y4)
# 		cv2.line(cl1, (int(x3), int(y3)),(int(x4), int(y4)),(0,0,255),1)


	# cv2.circle(mask, (int(midx1), int(midy1)), int(30),
	# 		(255, 0, 0), 2)
	# cv2.circle(mask, (int(midx2), int(midy2)), int(30),
	# 		(255, 0, 0), 2)
	# cv2.line(mask, (int(midx1)-50, int(midy1)),(int(midx2)+50, int(midy2)),(255,0,0),thickness=3)
	# cv2.line(cl1, (int(midx1)-50, int(midy1)),(int(midx2)+50, int(midy2)),(255,0,0),thickness=3)

	# cv2.drawContours(mask, [c], -1, (255, 0, 0), 2)
	# cv2.putText(image, shape, (cX, cY), cv2.FONT_HERSHEY_SIMPLEX,
	# 	0.5, (255, 255, 255), 2)

# find the contours in the mask, then sort them from left to
# right
# cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
# 	cv2.CHAIN_APPROX_SIMPLE)
# cnts = cnts[0] if imutils.is_cv2() else cnts[1]
# cnts = contours.sort_contours(cnts)[0]

# loop over the contours
# for (i, c) in enumerate(cnts):
# 	# draw the bright spot on the image
# 	(x, y, w, h) = cv2.boundingRect(c)
# 	((cX, cY), radius) = cv2.minEnclosingCircle(c)
# 	cv2.circle(mask, (int(cX), int(cY)), int(radius+3),
# 		(255, 0, 0), 2)
	# cv2.putText(mask, "#{}".format(i + 1), (x - 20, y),
	# 	cv2.FONT_HERSHEY_SIMPLEX, 0.45, (255, 0, 0), 2)

# show the output image
# res = np.hstack((cl1,mask))
# res = np.hstack((mask,maskc))
# cv2.imshow("Image", cl1)
# cv2.waitKey(0)
# cv2.destroyAllWindows()







#
