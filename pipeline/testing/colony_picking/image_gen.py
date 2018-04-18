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

def order_points(pts):
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
	image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
	blurred = cv2.GaussianBlur(image, (11, 11), 0)

	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(blurred)
	thresh = cv2.threshold(cl1, 100, 255, cv2.THRESH_BINARY)[1]

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
		# if numPixels > 50:
		if numPixels > 100000:
		# if numPixels < 50 and numPixels > 5:
			mask = cv2.add(mask, labelMask)

	cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_NONE)
	cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	maskc = np.zeros(thresh.shape, dtype="uint8")
	cv2.drawContours(maskc, cnts, -1, (255,0,0), 1)
	image = cv2.cvtColor(image,cv2.COLOR_GRAY2RGB)
	equ = []
	def find_lines(min):
		lines = cv2.HoughLines(maskc,1,np.pi/180,min)
		if min < 0:
			print("reached the minimum threshold")
			return lines
		if len(lines) < 4:
			return find_lines(min-5)
		elif len(lines) > 4:
			return find_lines(min+5)
		elif len(lines) == 4:
			return lines

	lines = find_lines(300)
	if len(lines) < 4:
		input("Still lacking enough lines")

	for line in lines:
		for rho,theta in line:
			a = np.cos(theta)
			b = np.sin(theta)
			x0 = a*rho
			y0 = b*rho
			equ += [[a,b,rho]]
			print((a/b),(rho/b))
			x1 = int(x0 + 5000*(-b))
			y1 = int(y0 + 5000*(a))
			x2 = int(x0 - 5000*(-b))
			y2 = int(y0 - 5000*(a))
			cv2.line(image,(x1,y1),(x2,y2),(255,0,0),2)

	pts = []
	for pair in itertools.combinations(equ,2):
		print(pair)
		pair = np.array(pair)
		a = pair[0:2,0:2]
		b = pair[0:2,2]
		if round(pair[0,0],2)*round(pair[1,1],2) == round(pair[1,0],2)*round(pair[0,1],2):
			print("lines are parallel")
			continue
		x = np.linalg.solve(a, b)
		pts += [x]
		cv2.circle(image, (int(x[0]), int(x[1])), int(10),
				(255, 0, 0), 2)
	pts = np.array(pts)
	print("These are the intersection points:\n",pts)
	return image, pts

def find_colonies(image):
	# create a CLAHE object (Arguments are optional).
	# resized = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(image)
	thresh = cv2.threshold(cl1, 180, 255, cv2.THRESH_BINARY)[1]

	# perform a connected component analysis on the thresholded
	# image, then initialize a mask to store only the "large"
	# components
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
		# if numPixels > 50:
		if numPixels > 30 and numPixels < 150:
		# if numPixels < 50 and numPixels > 5:
			mask = cv2.add(mask, labelMask)

	cnts = cv2.findContours(mask.copy(), cv2.RETR_EXTERNAL,
		cv2.CHAIN_APPROX_NONE)
	cnts = cnts[0] if imutils.is_cv2() else cnts[1]

	# edges = cv2.Canny(mask,50,150,apertureSize = 3)
	cl1 = cv2.cvtColor(cl1,cv2.COLOR_GRAY2RGB)

	# cv2.drawContours(cl1, cnts, -1, (255,0,0), 1)
	centers = []
	for c in cnts:
		# compute the center of the contour
		M = cv2.moments(c)
		if M["m00"] == 0:
			continue
		cX = int(M["m10"] / M["m00"])
		cY = int(M["m01"] / M["m00"])
		centers += [[cX,cY]]
		cv2.circle(cl1, (cX, cY), 1, (0, 0, 255), -1)
		# Rough approx if the colony is circular
		(x, y, w, h) = cv2.boundingRect(c)
		if w/h > 1.2 or h/w > 1.2:
			# print("not circular")
			continue
		((cX, cY), radius) = cv2.minEnclosingCircle(c)
		if int(radius) > 7:
			# print('too big')
			continue
		cv2.circle(cl1, (int(cX), int(cY)), int(radius+3),
			(255, 0, 0), 1)
	return cl1,centers

def find_grid(img):
	sums_y = []
	sums_x = []

	# create a CLAHE object (Arguments are optional).
	resized = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(resized)
	thresh = cv2.threshold(cl1, 180, 255, cv2.THRESH_BINARY)[1]

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
		if sum/255 < 20:
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

	# Calculates and plots the midpoints of the low regions
	mid_y = []
	for l,r in bounds:
		mid_point = int(((r-l) / 2)+l)
		mid_y.append(int(mid_point+top_bound/2))
		plt.plot(mid_point,sums_y[mid_point], 'bo')
		cv2.line(cl1,(0,int(mid_point+top_bound/2)),(int(cl1.shape[1]),int(mid_point+top_bound/2)),(0,255,0),1)
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

	# Calculates and plots the midpoints of the low regions
	mid_x = []
	for t,b in bounds:
		mid_point = int(((b-t) / 2)+t)
		mid_x.append(mid_point+top_bound)
		plt.plot(mid_point,sums_x[mid_point], 'bo')
		cv2.line(cl1,(int(mid_point+left_bound/2),0),(int(mid_point+left_bound/2),int(cl1.shape[0])),(255,0,0),1)
	print("Number of columns: ",len(mid_x)-1)

	# # Test rectangle
	cv2.rectangle(cl1,(mid_x[2],mid_y[4]),(mid_x[3],mid_y[8]),(0,0,255),3)

	plt.show()
	return cl1, mid_x, mid_y

def show_image(image):
	cv2.imshow("Image", image)
	cv2.waitKey(0)
	return

grids = []
shapes = []
for file in sorted(glob.glob('./image_set/*.jpg')):
	print(file)
	img = cv2.imread(file)
	resized = imutils.resize(img,width=1000)
	show_image(resized)
	edges, intersections = find_corners(resized)
	show_image(edges)
	warped = four_point_transform(edges,intersections)
	show_image(warped)
	grid, mid_x, mid_y = find_grid(warped)
	show_image(grid)
	colonies, centers = find_colonies(grid)
	show_image(colonies)
	print()

# print(shapes)
# input("stop")
# res = np.hstack(grids)
# cv2.imshow("Image", res)
# cv2.imshow("Image", warped)
# cv2.waitKey(0)
cv2.destroyAllWindows()

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
