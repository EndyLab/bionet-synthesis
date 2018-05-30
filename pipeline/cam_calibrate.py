import numpy as np
import cv2
import glob
import imutils


def take_calibration_photos():
    cam = cv2.VideoCapture(0)

    # termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)
    # prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
    objp = np.zeros((6*7,3), np.float32)
    objp[:,:2] = np.mgrid[0:7,0:6].T.reshape(-1,2)
    # Arrays to store object points and image points from all the images.
    objpoints = [] # 3d point in real world space
    imgpoints = [] # 2d points in image plane.

    counter = 12

    while(True):
        if counter == 24:
            break

        ret, frame = cam.read()

        frame = imutils.resize(frame, width=500)
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        img = np.copy(frame)

        cv2.imshow('img', img)

        # Find the chess board corners
        ret, corners = cv2.findChessboardCorners(gray, (7,6), None)
        # If found, add object points, image points (after refining them)
        if ret == True:
            # counter += 1
            objpoints.append(objp)
            corners2=cv2.cornerSubPix(gray,corners, (11,11), (-1,-1), criteria)
            imgpoints.append(corners)

            # Draw and display the corners
            cv2.drawChessboardCorners(img, (7,6), corners2, ret)
            cv2.imshow('img', img)
            cv2.imwrite('./cali_photos/photo{}.jpg'.format(str(counter)),frame)
            cv2.waitKey(0)

        if cv2.waitKey(1) & 0xFF == ord('q'):
            break

def calibrate_cam(images,output):
    # termination criteria
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

    # prepare object points, like (0,0,0), (1,0,0), (2,0,0) ....,(6,5,0)
    objp = np.zeros((6*7,3), np.float32)
    objp[:,:2] = np.mgrid[0:7,0:6].T.reshape(-1,2)

    # Arrays to store object points and image points from all the images.
    objpoints = [] # 3d point in real world space
    imgpoints = [] # 2d points in image plane.

    for fname in images:
        print(fname)

        img = cv2.imread(fname)
        resized = imutils.resize(img, width=1920)
        gray = cv2.cvtColor(resized, cv2.COLOR_BGR2GRAY)

        # Find the chess board corners
        ret, corners = cv2.findChessboardCorners(gray, (7,6), None)

        # If found, add object points, image points (after refining them)
        if ret == True:
            print("Found chess board")
            objpoints.append(objp)
            corners2=cv2.cornerSubPix(gray,corners, (11,11), (-1,-1), criteria)
            imgpoints.append(corners)

            # Draw and display the corners
            cv2.drawChessboardCorners(resized, (7,6), corners2, ret)
            cv2.imshow('img', resized)
            cv2.waitKey(500)
        else:
            print("not found")

    cv2.destroyAllWindows()

    # Calibrate
    ret, mtx, dist, rvecs, tvecs = cv2.calibrateCamera(objpoints, imgpoints, gray.shape[::-1], None, None)
    np.savez(output,ret=ret,mtx=mtx,dist=dist,rvecs=rvecs,tvecs=tvecs)

    return # ret, mtx, dist, rvecs, tvecs

def undistort_img(img,cali_file):# ret, mtx, dist, rvecs, tvecs):
    # for fname in images:
    h,  w = img.shape[:2]
    calibrations = np.load(cali_file)
    newcameramtx, roi=cv2.getOptimalNewCameraMatrix(calibrations['mtx'], calibrations['dist'], (w,h), 1, (w,h))

    # undistort
    mapx, mapy = cv2.initUndistortRectifyMap(calibrations['mtx'], calibrations['dist'], None, newcameramtx, (w,h), 5)
    dst = cv2.remap(img, mapx, mapy, cv2.INTER_LINEAR)

    # crop the image
    x, y, w, h = roi
    dst = dst[y:y+h, x:x+w]

    return dst

if __name__ == "__main__":
    cali_images = sorted(glob.glob('./cali_photos/*.jpg'))
    cali_file = './webcam_calibrations.npz'
    calibrate_cam(cali_images,cali_file)
    images = sorted(glob.glob('./cam_photos/*.jpg'))
    for fname in images:
        img = cv2.imread(fname)
        print(img.shape[0],img.shape[1])
        # img = imutils.resize(img, width=500)
        dst = undistort_img(img,cali_file)
        cv2.imshow('original',img)
        cv2.waitKey(0)
        cv2.imshow('undistorted',dst)
        cv2.waitKey(0)
        cv2.destroyAllWindows()






#
