from matplotlib import pyplot as plt
import seaborn as sns
from scipy import interpolate

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
import re
import getch
import math
from itertools import chain
from datetime import datetime

from opentrons import robot, containers, instruments
from opentrons.helpers import helpers

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

def apply_threshold(image,size=201):
	'''Applies an adaptive threshold to the image and returns the result'''

	# The size and C were derived emperically through studying their effect on
	# the specific images used here. May need to be modified for other cameras
	thresh = cv2.adaptiveThreshold(image,255,cv2.ADAPTIVE_THRESH_MEAN_C,\
                cv2.THRESH_BINARY,size,-8)
	return thresh


def find_corners(image):
	'''
	Searches the image to find the edges of the agar plate
	'''
	# Disperse high frequency noise to find a cleaner edge
	image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
	blurred = cv2.GaussianBlur(image, (11, 11), 0)

	# Apply a tiling contrast enhancer to allow for more accurate thresholding
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(blurred)
	# show_image(cl1)

	# Apply a threshold that should only show the plate
	thresh = cv2.threshold(cl1, 50, 255, cv2.THRESH_BINARY)[1]
	# show_image(thresh)

	# Generate a mask to add only the plate back in
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

		# Only the complete plate should have enough pixels to clear this threshold
		if numPixels > 100000:
			mask = cv2.add(mask, labelMask)
	# show_image(mask)

	# Find the contour of the plate
	cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_NONE)
	cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	# Find the edges of the contour and draw them on the image
	epsilon = 0.1*cv2.arcLength(cnts[0],True)
	pts = cv2.approxPolyDP(cnts[0],epsilon,True)
	image = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
	cv2.polylines(image,[pts],True,(0,255,255),4)

	# Unnest the coordinates of the corners to create a numpy array
	corners = [[x,y] for [[x,y]] in pts]
	corners = np.array(corners)
	print("Corners:",corners)

	return image, corners

def find_minima(sums,number,threshold=40):
	# Sets a threshold to define regions considered 'empty' and stores their indexes
	bottoms = []
	for i,sum in enumerate(sums):
		if sum/255 < threshold:
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

	# Groups the endpoints together to define the bounds of each section
	bounds = [edges[n:n+2] for n in range(0, len(edges), 2)]
	print("Bounds: ",bounds)

	single = [bound for bound in bounds if len(bound) < 2]
	print(single,len(single))

	if len(single) > 0:
		print('single',threshold)
		return find_minima(sums,number,threshold=threshold+1)

	# Calculates and plots the midpoints of the low regions
	mid = []
	for l,r in bounds:
		mid_point = int(((r-l) / 2)+l)
		mid.append(int(mid_point))

	if len(mid)-1 > number:
		print('Too many minima',threshold)
		return find_minima(sums,number,threshold=threshold-5)
	if len(mid)-1 < number:
		print('Too few minima',threshold)
		return find_minima(sums,number,threshold=threshold+5)

	print("Number of minima: ",len(mid)-1)
	print("Threshold: ",threshold)

	return mid


def find_grid(img):
	'''
	Scans the image for the rows and columns of the image with the fewest colonies
	to define a grid segmenting the plate based on part and dilution
	'''
	sums_y = []
	sums_x = []

	# Gray-scale, enhance contrast and threshold to only see the colonies
	gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
	thresh = apply_threshold(gray)
	show_image(thresh)

	# Find the sum of the points in the rows and cols and store them in lists
	height,width = thresh.shape[:2]
	for row in range(height):
		sum = np.sum(thresh[row,0:width])
		sums_y.append(sum)
	for col in range(width):
		sum = np.sum(thresh[0:height,col])
		sums_x.append(sum)

	mid_y = find_minima(sums_y,12)
	mid_x = find_minima(sums_x,8)

	return gray, mid_x, mid_y

def group_sections(img,mid_x,mid_y,x_steps=1,y_steps=4):
	'''
	Groups the small section of the grid by the sections corresponding to
	a single gene
	'''
	groups = []
	for y in range(0,len(mid_y)-1,y_steps):
		for x in range(0,len(mid_x)-1,x_steps):
			group = [[mid_x[x],mid_x[x+x_steps]],[mid_y[y],mid_y[y+y_steps]]]
			print(group)
			groups.append(group)
			cv2.rectangle(img,(mid_x[x],mid_y[y]),(mid_x[x+x_steps],mid_y[y+y_steps]),(0,0,255),3)
			cv2.imshow("Image", img)
			cv2.waitKey(10)
		# cv2.imwrite('saved_pics/image{}.png'.format(datetime.now()),img)
	return groups

def show_small(img):
	'''Displays images at a smaller width'''
	resized = imutils.resize(img,width=100)
	cv2.imshow("Image", resized)
	# cv2.imwrite('saved_pics/image{}.png'.format(datetime.now()),img)
	cv2.waitKey(1)
	input("next image")

def find_colonies(image,check,size=201):
	'''
	Searches the image for pickable colonies and returns the
	location of their centers
	'''

	thresh = apply_threshold(image,size=size)
	show_small(thresh)
	show_small(image)

	# Find the contours of each colony
	cnts = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_NONE)
	cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	# Convert the color to be able to draw on them in color
	color = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
	color_th = cv2.cvtColor(thresh,cv2.COLOR_GRAY2RGB)

	## Iterate over each colony, extract several different parameters, then
	## Determine if they are acceptable colonies

	centers = []
	good_cols = []
	bad_cols = []
	good_center = []
	for counter,cnt in enumerate(cnts):
		print("counter: ",counter)

		# Find their moments and use it to find the center of mass
		M = cv2.moments(cnt)
		if M["m00"] == 0:
			continue
		cX = int(M["m10"] / M["m00"])
		cY = int(M["m01"] / M["m00"])

		# Find the perimeter and use it to generate an approximation of the shape
		perimeter = cv2.arcLength(cnt,True)
		approx = cv2.approxPolyDP(cnt,0.01*perimeter,True)

		# Use the approximate shape to more reliably find colonies
		perimeter = cv2.arcLength(approx,True)
		area = cv2.contourArea(approx)

		# Calc represents the proportion of area to calculated area based on
		# perimeter
		calc = area/((perimeter**2)/(4*math.pi))

		# Generates a bound rectangle that can rotate so the ratio of its width to
		# height is a metric for circularity
		(rX,rY),(w,h),ang = cv2.minAreaRect(approx)
		ratio = w/h

		attributes = [area,perimeter,calc,ratio]

		# Set the acceptable ranges for the different attributes. Values were
		# determined by generating a test set and extracting the acceptable ranges
		area_r = [10,375]
		peri_r = [12.5,80]
		calc_r = [0.55,1]
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

		# Assess and draw a green contour for good and red for bad
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
		# Currently picks the good colony closest to the top of the plate
		# TODO: Modify the selection process to allow for sorting by proximity
		good_center = centers[-1]
	cv2.drawContours(color, good_cols, -1, (0,255,0), 1)
	cv2.drawContours(color_th, good_cols, -1, (0,255,0), 1)
	cv2.drawContours(color, bad_cols, -1, (0,0,255), 1)
	cv2.drawContours(color_th, bad_cols, -1, (0,0,255), 1)
	# show_small(color_th)
	# show_small(color)

	if size < 100:
		print("No colonies found")
	elif len(centers) == 0:
		print("No colonies found, changing threshold")
		return (find_colonies(image,check,size=size-50))
	else:
		# Currently picks the good colony closest to the top of the plate
		def choose_colony(img,centers):
			print("Enter 'g'->good, 'b'->bad")
			show_small(color)
			x = len(centers)
			counter = 0
			temp_o = np.copy(img)
			while x != 0:
				temp = np.copy(temp_o)
				counter -= 1
				good_center = centers[counter]
				print(good_center)
				cv2.circle(temp, (int(good_center[0][0]),int(good_center[0][1])), int(20),
					(0,255,255), 2)
				resized = imutils.resize(temp,width=100)
				cv2.imshow("Image", resized)
				cv2.waitKey(10)
				c = getch.getch()
				if c == "g":
					print("Good")
					return good_center
				elif c == "b":
					print("Going to next option")
					x -= 1
				else:
					print("Going to next option")
					x -= 1
			print('Iterated through all colonies')
			ans = str(input('r-retry, q-quit'))
			if ans == 'r':
				return choose_colony(img,centers)
			elif ans == 'q':
				return []
		if check:
			good_center = choose_colony(color,centers)
		else:
			good_center = centers[-1]

	return good_center

def change_loc(image,start):
	'''
	Allows for manually modification of points presented on images
	'''
	# Shows the initial starting position
	x = 0
	temp = np.copy(image)
	cv2.line(temp,(int(start[0]),int(start[1])-40),(int(start[0]),int(start[1])+40),(0,0,255),2)
	cv2.line(temp,(int(start[0])-40,int(start[1])),(int(start[0])+40,int(start[1])),(0,0,255),2)
	temp_resized = imutils.resize(temp,width=400)
	cv2.imshow("Image", temp_resized)
	cv2.waitKey(5)

	# Allows the user to move the position of the cursor to ensure that the location
	# is correct
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
	# The agar outline
	all_points = [[19,27],[377,27],[376,605],[19, 610]]
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

def ot_coords(centers,image,x_dim,y_dim):
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

def move_motor(pipette):
	'''Allows for real-time calibration of the p10single channel xy location'''
	z = 0
	print("Change xy - w:back, d:forward, a:left, d:right, x:exit")
	x = 0
	y = 0
	while z == 0:
		c = getch.getch()
		if c == "w":
			pipette.robot._driver.move(y=0.5,mode="relative")
			print("up")
			y += 0.5
		elif c == "s":
			pipette.robot._driver.move(y=-0.5,mode="relative")
			print("down")
			y -= 0.5
		elif c == "a":
			pipette.robot._driver.move(x=-0.5,mode="relative")
			print("left")
			x -= 0.5
		elif c == "d":
			pipette.robot._driver.move(x=0.5,mode="relative")
			print("right")
			x += 0.5
		elif c == "x":
			z = 1
	print(x,y)
	return x,y

def find_offset(x,y,z):
	f = interpolate.interp2d(x, y, z, kind='linear')

	xnew = np.arange(min(x), max(x))
	ynew = np.arange(min(y), max(y))

	xx, yy = np.meshgrid(xnew, ynew)

	znew = f(xnew, ynew)
	z_flat = [z for sublist in znew.tolist() for z in sublist]

	xx, yy = np.meshgrid(xnew, ynew)
	x_flat = [x for sublist in xx.tolist() for x in sublist]
	y_flat = [y for sublist in yy.tolist() for y in sublist]

	all_data = pd.DataFrame({
		'X':x_flat,
		'Y':y_flat,
		'Offset':z_flat
	})
	result = all_data.pivot(index='Y', columns='X', values='Offset')
	fig, ax = plt.subplots()

	ax = sns.heatmap(result,ax=ax,xticklabels=10,yticklabels=10)
	xt = np.arange(int(min(x)), int(max(x)), step=10)
	xl = [str(num) for num in xt]
	plt.xticks(xt,xl)

	fig.show()
	print('Y-axis is reversed')
	input("Check heatmap")
	plt.close()

	return f

def calibrate_ot(pipette,agar_plate,image,coords,centers):
	no_well = [[x,y] for [x,y,w] in centers]
	bl_O, tl_O, tr_O, br_O = order_points(coords)
	tl_I, tr_I, br_I, bl_I = order_points(no_well)

	paired_coords = [[bl_O,bl_I],[tl_O,tl_I],[tr_O,tr_I],[br_O,br_I]]

	temp = np.copy(image)
	Ox_l = []
	Oy_l = []
	Ozx = []
	Ozy = []
	pipette.pick_up_tip()
	print('picking up tip')
	for [[Ox,Oy],[Ix,Iy]] in paired_coords:
		cv2.line(temp,(int(Ix),int(Iy)-30),(int(Ix),int(Iy)+30),(0,0,255),2)
		cv2.line(temp,(int(Ix)-30,int(Iy)),(int(Ix)+30,int(Iy)),(0,0,255),2)
		show_image(temp)
		pipette.move_to((agar_plate,[Ox,Oy,0]))
		print("Calibrate the position for the colony")
		off_x,off_y = move_motor(pipette)
		Ox_l.append(Ox)
		Oy_l.append(Oy)
		Ozx.append(off_x)
		Ozy.append(off_y)
	print(paired_coords)
	fx = find_offset(Ox_l,Oy_l,Ozx)
	fy = find_offset(Ox_l,Oy_l,Ozy)


	return fx,fy

def mix_in_well(pipette,depth=-0.75,location=None,radius=0.7):
    well_edges = [
        location.from_center(x=radius, y=0, z=depth),       # right edge
        location.from_center(x=0, y=radius, z=depth),       # back edge
        location.from_center(x=radius * -1, y=0, z=depth),  # left edge
        location.from_center(x=0, y=radius * -1, z=depth)   # front edge
    ]
    [pipette.move_to((location, e), strategy='direct') for e in well_edges]

    return

def inoculate(pipette,agar_plate,colony,well,depth=-0.75,radius=0.7,mix=3):

	_description = 'Inoculating'
	pipette.robot.add_command(_description)

	x,y,z = colony

	if not pipette.current_tip():
		print('no tip')
		pipette.pick_up_tip()
		print('picking up tip')
	else:
		print('Already has a tip')
		print('Calibrate Z')
		print('Use "r"->up, "f"->down to pick colony')
		pipette.move_to((agar_plate,[x,y,z]), strategy='arc')
		while True:
			c = getch.getch()
			if c == 'r':
				z += 0.5
				print('up')
			elif c == 'f':
				z -= 0.5
				print('down')
			elif c == 'x':
				break
			pipette.move_to((agar_plate,[x,y,z]), strategy='direct')
		print('Final z:', z)

	print(x,y,z)
	pipette.move_to((agar_plate,[x,y,z]), strategy='arc')
	input('Picked colony?')

	pipette.move_to(well, strategy='arc')

	for num in range(mix):
		mix_in_well(pipette,location=well,depth=depth,radius=radius)

	pipette.move_to((well,well.from_center(x=0, y=0, z=depth)),strategy='direct')

	pipette.motor.move(pipette._get_plunger_position('drop_tip'))
	pipette.motor.move(pipette._get_plunger_position('bottom'))

	pipette.current_volume = 0
	pipette.current_tip(None)

	return z


def run_ot(pipette,agar_plate,deep_well,image,coords,centers,pick):
	'''
	Pass the coordinates to the robot to pick the colony
	'''

	fx,fy = calibrate_ot(pipette,agar_plate,image,coords,centers)
	z = 0
	for i,((Ox,Oy),(Ix,Iy,well)) in enumerate(zip(coords,centers)):
		temp = image
		print(i,":",Ox,Oy)
		x_off = fx(Ox,Oy)[0]
		y_off = fy(Ox,Oy)[0]
		new_x = Ox + x_off
		new_y = Oy + y_off
		print(i,":",new_x,new_y)
		cv2.circle(temp, (int(Ix),int(Iy)), int(20),
			(0,0,255), 2)
		show_image(temp)
		print(well)
		colony = [new_x,new_y,z]
		if pick:
			print('inoculating')
			z = inoculate(pipette,agar_plate,colony,deep_well.wells(well))
		else:
			pipette.move_to((agar_plate,[new_x,new_y,z]))
		input("click to move to next colony")

def show_image(image,width=400):
	'''Scales and then presents the image'''
	resized = imutils.resize(image,width=width)
	cv2.imshow("Image", resized)
	# cv2.imwrite('saved_pics/image{}.png'.format(datetime.now()),image)
	cv2.waitKey(5)
	input('Enter to move to next image')
	return

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

def find_colony_coordinates(file,check):
	build_num,plate_num = re.match(r'.*\/.+.([0-9]{3})p([0-9]).jpg',file).groups()
	print('Build: {}, Plate: {}'.format(build_num,plate_num))
	print((int(plate_num)-1)*24,int(plate_num)*24)
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
	groups = groups[16:24] + groups[8:16] + groups[:8]
	centers = []
	missing = []
	wells = well_addresses()[(int(plate_num)-1)*24:int(plate_num)*24]
	print(wells)
	input("check wells")
	agar = cv2.cvtColor(agar, cv2.COLOR_BGR2GRAY)
	for well,group in zip(wells,groups):
		print("Next section:",group)
		partial = agar[group[1][0]:group[1][1],group[0][0]:group[0][1]]
		center = find_colonies(partial,check)
		if center == []:
			midg_y = int(((group[1][0]-group[0][0])/2)+group[0][0])
			midg_x = int(((group[1][1]-group[0][1])/2)+group[0][1])
			cv2.rectangle(color,(group[0][0],group[1][0]),(group[0][1],group[1][1]),(0,0,255),3)
			missing.append(well)
		else:
			cali_cenX = center[0][0]+group[0][0]
			cali_cenY = center[0][1]+group[1][0]
			cv2.circle(color, (cali_cenX,cali_cenY), int(20),
				(0,255,255), 3)
			cv2.circle(color, (cali_cenX,cali_cenY), int(2),
				(0,0,255), 2)
			centers.append([cali_cenX,cali_cenY,well])
	show_image(color)
	input("check image")
	agar = cv2.cvtColor(agar,cv2.COLOR_GRAY2RGB)
	return centers,agar,x_dim,y_dim,missing



def pick_colonies():

	parser = argparse.ArgumentParser(description="Resuspend a plate of DNA on an Opentrons OT-1 robot.")
	parser.add_argument('-r', '--run', required=False, action="store_true", help="Send commands to the robot and print command output.")
	parser.add_argument('-i', '--inoculate', required=False, action="store_true", help="Picks colonies and inoculates a deep well plate.")
	parser.add_argument('-c', '--check', required=False, action="store_true", help="Allows the user to choose the colonies.")
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
	deep_well = containers.load('96-deep-well', 'C2')
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

	for file in sorted(glob.glob('./new_photos/*.jpg'),reverse=True):
		print(file)
		skip = input('Skip? y/n: ')
		if skip == 'y':
			continue
		print("Inoculate: ",args.inoculate)
		input()
		centers,agar,x_dim,y_dim,missing = find_colony_coordinates(file,args.check)
		coords = ot_coords(centers,agar,x_dim,y_dim)
		run_ot(p10s,trans_plate,deep_well,agar,coords,centers,args.inoculate)
		print('The following wells were skipped:\n',missing)
		print("Complete")
		robot.home()
	cv2.destroyAllWindows()

if __name__ == "__main__":
    pick_colonies()






#
