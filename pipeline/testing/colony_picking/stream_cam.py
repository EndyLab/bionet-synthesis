import numpy as np
import cv2
import imutils

cam = cv2.VideoCapture(0)

def stream_cam():
    counter = 0
    while(True):
        ret, frame = cam.read()

        img = np.copy(frame)
        img = imutils.resize(img, width=1000)
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

        cv2.imshow('img', gray)

        if cv2.waitKey(1) & 0xFF == ord('q'):
            counter += 1
            cv2.imwrite('./cam_photos/photo{}.jpg'.format(str(counter)),frame)
            cv2.waitKey(0)
            # break

if __name__ == "__main__":
    stream_cam()


# When everything done, release the capture
cam.release()
cv2.destroyAllWindows()
