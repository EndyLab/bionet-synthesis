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
	# show_image(thresh)
	# Set bounds on all of the edges to crop out the edge of the plate
	left_bound = int(thresh.shape[1]/10)
	right_bound = int((9*thresh.shape[1])/10)
	top_bound = int(thresh.shape[0]/10)
	bot_bound = int((9*thresh.shape[0])/10)

	# Sum all of the rows in he thresholded image
	for row in range(int(top_bound/2),int(bot_bound+top_bound/2)):
		sum = np.sum(thresh[row,left_bound:right_bound])
		sums_y.append(sum)

	# Sum all of the columns in he thresholded image
	for col in range(int(left_bound/2),int(right_bound+left_bound/2)):
		sum = np.sum(thresh[top_bound:bot_bound,col])
		sums_x.append(sum)

	# cl1 = cv2.cvtColor(cl1,cv2.COLOR_GRAY2RGB)

	# Generate a plot to visualize the sumations
	plt.subplot(2, 1, 1)
	plt.plot(sums_y)
	plt.ylabel('Row sum')

	# Store only the points that have less than 20 pixels in the line
	bottoms = []
	for i,sum in enumerate(sums_y):
		if sum/255 < 30:
			plt.plot(i,sum, 'r.')
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

	temp = np.copy(cl1)

	# Calculates and plots the midpoints of the low regions
	mid_y = []
	for l,r in bounds:
		mid_point = int(((r-l) / 2)+l)
		mid_y.append(int(mid_point+top_bound/2))
		plt.plot(mid_point,sums_y[mid_point], 'bo')
		cv2.line(temp,(0,int(mid_point+top_bound/2)),(int(cl1.shape[1]),int(mid_point+top_bound/2)),(0,255,0),1)
	print("Number of rows: ",len(mid_y)-1)

	# Generate a plot to visualize the sumations
	plt.subplot(2, 1, 2)
	plt.plot(sums_x)
	plt.ylabel('Row sum')

	# Store only the points that have less than 60 pixels in the line
	bottoms = []
	for i,sum in enumerate(sums_x):
		if sum/255 < 60:
			plt.plot(i,sum, 'r.')
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
		mid_x.append(mid_point+top_bound)
		plt.plot(mid_point,sums_x[mid_point], 'bo')
		cv2.line(temp,(int(mid_point+left_bound/2),0),(int(mid_point+left_bound/2),int(cl1.shape[0])),(255,0,0),1)

	print("Number of columns: ",len(mid_x)-1)
	show_image(temp)
	# # Test rectangle
	# cv2.rectangle(cl1,(mid_x[2],mid_y[4]),(mid_x[3],mid_y[8]),(0,0,255),3)

	plt.show()
	return cl1, mid_x, mid_y

def find_colonies(image):
	'''
	Searches the image for pickable colonies and returns the
	location of their centers
	'''
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(image)
	thresh = cv2.threshold(cl1, 180, 255, cv2.THRESH_BINARY)[1]


	cnts = cv2.findContours(thresh.copy(), cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_NONE)
	cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	# edges = cv2.Canny(mask,50,150,apertureSize = 3)
	color = cv2.cvtColor(cl1,cv2.COLOR_GRAY2RGB)
	cv2.drawContours(color, cnts, -1, (255,0,0), 1)

	# cv2.drawContours(cl1, cnts, -1, (255,0,0), 1)
	centers = []
	for counter,cnt in enumerate(cnts):
		# compute the center of the contour
		# print("counter: ",counter)
		M = cv2.moments(cnt)
		if M["m00"] == 0:
			continue
		cX = int(M["m10"] / M["m00"])
		cY = int(M["m01"] / M["m00"])

		perimeter = cv2.arcLength(cnt,True)
		approx = cv2.approxPolyDP(cnt,0.01*perimeter,True)
		area = cv2.contourArea(cnt)

		(rX,rY),(w,h),ang = cv2.minAreaRect(cnt)

		if area > 50 and area < 350 and w/h > 0.8 and w/h < 1.25 and area > perimeter:
			cv2.circle(color, (cX, cY), 1, (0, 0, 255), -3)
			# cv2.putText(color, "{}".format(counter), (cX - 5, cY - 5),
			# 	cv2.FONT_HERSHEY_SIMPLEX, 0.45, (0, 255, 255), 3)
			((cX, cY), radius) = cv2.minEnclosingCircle(cnt)
			cv2.circle(color, (int(cX), int(cY)), int(radius+9),
				(0, 255, 255), 3)

			rect = cv2.minAreaRect(cnt)

			print(counter,".","Area",area,"Peri",perimeter,"comp",(perimeter**2)/(4*math.pi),"Approx",len(approx))
			print(w,h,w/h)
			show_image(color)
			centers += [[cX,cY]]
	return color,centers

def change_loc(image,start):
	'''
	Allows for manually modification of points presented on images
	'''
	x = 0
	temp = np.copy(image)
	cv2.circle(temp, (int(start[0]),int(start[1])), int(3),
		(0,0,255), 9)
	cv2.circle(temp, (int(start[0]),int(start[1])), int(30),
		(255,0,0), 9)
	temp_resized = imutils.resize(temp,width=400)
	cv2.imshow("Image", temp_resized)
	cv2.waitKey(5)
	print("please fix the location")
	while x == 0:
		# temp = image
		c = getch.getch()
		if c == "w":
			start[1] -= 3
			print("up")
		elif c == "s":
			start[1] += 3
			print("down")
		elif c == "a":
			start[0] -= 3
			print("left")
		elif c == "d":
			start[0] += 3
			print("right")
		elif c == "x":
			x = 1
			print("exit")
		temp = np.copy(image)
		cv2.circle(temp, (int(start[0]),int(start[1])), int(3),
			(0,0,255), 9)
		cv2.circle(temp, (int(start[0]),int(start[1])), int(30),
			(255,0,0), 9)
		temp_resized = imutils.resize(temp,width=400)
		cv2.imshow("Image", temp_resized)
		cv2.waitKey(5)

	return start


def find_reference(image):
	'''
	Locate the corners of the agar to link the reference point on the robot
	to the image
	'''
	image = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)

	# all_points = [[22,19],[475,19],[475,305],[22, 305]]
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
	x_dim = new_points[2][0] - new_points[0][0]
	y_dim = new_points[2][1] - new_points[0][1]

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
		if x < 0 or y < 0:
			print("invalid coordinate")
			continue
		# print(cen[0],ref[0],x_dim,x_max,"=",x)
		# print(cen[1],ref[1],y_dim,y_max,"=",y)
		# input()
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
			y += 0.5
		elif c == "s":
			p10s.robot._driver.move(y=-0.5,mode="relative")
			y -= 0.5
		elif c == "a":
			p10s.robot._driver.move(x=-0.5,mode="relative")
			x -= 0.5
		elif c == "d":
			p10s.robot._driver.move(x=0.5,mode="relative")
			x += 0.5
		elif c == "x":
			z = 1
	print(x,y)
	# p10s.calibrate_position((container,target.from_center(x=0, y=0, z=-1,reference=container)))


def run_ot(image,coords,centers):
	'''
	Pass the coordinates to the robot to pick the colony
	'''
	for i,(coord,cen) in enumerate(zip(coords,centers)):
		temp = image
		print(i,":",coord[0],coord[1])
		cv2.circle(temp, (int(cen[0]),int(cen[1])), int(20),
			(0,0,255), 3)
		show_image(temp)
		p10s.move_to((trans_plate,[coord[0]-0.5,coord[1]-2,0]))
		move_motor()


		# input("click to move to next colony")

#
def show_image(image,width=400):
	'''Scales and then presents the image'''
	resized = imutils.resize(image,width=width)
	cv2.imshow("Image", resized)
	cv2.waitKey(1)
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
	# dst = undistort_img(img,cali_file)
	# show_image(dst)
	# resized = imutils.resize(new,width=10000)
	# show_image(resized)
	edges, intersections = find_corners(rotated)
	show_image(edges)
	warped = four_point_transform(edges,intersections)
	show_image(warped)
	grid, mid_x, mid_y = find_grid(warped)
	show_image(grid)
	ref, x_dim, y_dim,rect = find_reference(grid)
	agar = four_point_transform(grid,rect)
	show_image(agar)

	colonies, centers = find_colonies(agar)
	# coords = ot_coords(centers,ref, x_dim, y_dim)
	coords = ot_coords(centers,agar)
	run_ot(colonies,coords,centers)
	# show_image(colonies)
	print()

# print(shapes)
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
