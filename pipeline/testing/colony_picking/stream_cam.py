import numpy as np
import cv2
import imutils

def stream_cam():
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
            cv2.imwrite('./new_photos/new_build_{}.jpg'.format(name),frame)
            cv2.waitKey(10)
            counter += 1
            if counter == 5:
                break

    cam.release()
    cv2.destroyAllWindows()

if __name__ == "__main__":
    stream_cam()
