from matplotlib import pyplot as plt

from imutils import contours
from skimage import measure
import numpy as np
import itertools
import argparse
import imutils
import cv2

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
	# lines = cv2.HoughLines(maskc,1,np.pi/180,300)
	def find_lines(min):
		lines = cv2.HoughLines(maskc,1,np.pi/180,min)
		# print("Number of lines:",len(lines))
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
			# print((a/b),(rho/b))
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
		# print(a,b)
		if round(pair[0,0],2)*round(pair[1,1],2) == round(pair[1,0],2)*round(pair[0,1],2):
			print("lines are parallel")
			continue
		x = np.linalg.solve(a, b)
		# print(x)
		pts += [x]
		cv2.circle(image, (int(x[0]), int(x[1])), int(10),
				(255, 0, 0), 2)
	inter = {}
	print(pts)
	pts = np.array(pts)
	print(pts)
	pts = pts[pts[:,0].argsort()]
	print(pts)
	if pts[0,1] < pts[1,1]:
		inter.update({'topL' : {'x' : pts[0,0],'y' : pts[0,1]}})
		inter.update({'botL' : {'x' : pts[1,0],'y' : pts[1,1]}})
	else:
		inter.update({'botL' : {'x' : pts[0,0],'y' : pts[0,1]}})
		inter.update({'topL' : {'x' : pts[1,0],'y' : pts[1,1]}})
	if pts[2,1] < pts[3,1]:
		inter.update({'topR' : {'x' : pts[2,0],'y' : pts[2,1]}})
		inter.update({'botR' : {'x' : pts[3,0],'y' : pts[3,1]}})
	else:
		inter.update({'botR' : {'x' : pts[2,0],'y' : pts[2,1]}})
		inter.update({'topR' : {'x' : pts[3,0],'y' : pts[3,1]}})
	print(inter)
	print(inter['topR']['x'])
	return image,inter

def find_colonies(image):
	# create a CLAHE object (Arguments are optional).
	resized = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
	clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
	cl1 = clahe.apply(resized)
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
			print("not circular")
			continue
		((cX, cY), radius) = cv2.minEnclosingCircle(c)
		if int(radius) > 7:
			print('too big')
			continue
		cv2.circle(cl1, (int(cX), int(cY)), int(radius+3),
			(255, 0, 0), 1)
	return cl1,centers

def make_grid(image,p):
	# Find the change in pixels in both directions across all of the sides of the plate
	top_lenX = abs(p['topR']['x'] - p['topL']['x'])
	top_lenY = abs(p['topR']['y'] - p['topL']['y'])
	offsetT = top_lenX / 10
	top_lenX -= offsetT * 2

	left_lenX = abs(p['botL']['x'] - p['topL']['x'])
	left_lenY = abs(p['botL']['y'] - p['topL']['y'])
	offsetL = left_lenY / 12
	left_lenY -= offsetL * 2

	right_lenX = abs(p['topR']['x'] - p['botR']['x'])
	right_lenY = abs(p['topR']['y'] - p['botR']['y'])
	offsetR = right_lenY / 12
	right_lenY -= offsetR * 2

	bot_lenX = abs(p['botL']['x'] - p['botR']['x'])
	bot_lenY = abs(p['botL']['y'] - p['botR']['y'])
	offsetB = bot_lenX / 10
	bot_lenX -= offsetB * 2

	y1 = int(offsetL+p['topL']['y'])
	x1 = int(p['topL']['x'])
	y2 = int(offsetR+p['topR']['y'])
	x2 = int(p['topR']['x'])
	cv2.line(image, (x1, y1),(x2, y2),(255,0,0),1)

	num_rows = 12
	for row in range(num_rows):
		y1 = (left_lenY - left_lenY*(row/num_rows)) + p['topL']['y'] + offsetL
		x1 = (left_lenX - left_lenX*(row/num_rows)) + p['topL']['x']
		y2 = (right_lenY - right_lenY*(row/num_rows)) + p['topR']['y'] + offsetR
		x2 = (right_lenX - right_lenX*(row/num_rows)) + p['topR']['x']
		print(x1,y1,x2,y2)

		cv2.line(image, (int(x1), int(y1)),(int(x2), int(y2)),(255,0,0),1)
	num_cols = 8
	cv2.line(image, (int(offsetT+p['topL']['x']), int(p['topL']['y'])),(int(offsetB+p['botL']['x']), int(p['botL']['y'])),(0,0,255),1)
	for row in range(num_cols):
		y1 = (bot_lenY - bot_lenY*(row/num_cols)) + p['topL']['y']
		x1 = (bot_lenX - bot_lenX*(row/num_cols)) + p['topL']['x'] + offsetB
		y2 = (top_lenY - top_lenY*(row/num_cols)) + p['botL']['y']
		x2 = (top_lenX - top_lenX*(row/num_cols)) + p['botL']['x'] + offsetT
		print(x1,y1,x2,y2)

		cv2.line(image, (int(x1), int(y1)),(int(x2), int(y2)),(0,0,255),1)
	return image


img = cv2.imread('plate3.jpg')
resized = imutils.resize(img,width=1000)

edges, intersections = find_corners(resized)
colonies, centers = find_colonies(resized)
grid = make_grid(edges,intersections)

res = np.hstack((grid,colonies))
cv2.imshow("Image", res)
cv2.waitKey(0)
cv2.destroyAllWindows()


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
