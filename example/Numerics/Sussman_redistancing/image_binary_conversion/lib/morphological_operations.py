import cv2
import numpy as np

def closure(img, kernel, iterations):
    lz, ly, lx = img.shape[0], img.shape[1], img.shape[2]
    erosion = np.zeros((lz, ly, lx))
    dilation = np.zeros((lz, ly, lx))
    
    # Closure
    # -> connects separated regions, increase connectivity
    i = 0
    while i < 1:

        for z_slice in range(lz):
            dilation[z_slice] = cv2.dilate(img[z_slice], kernel, iterations = iterations)
        img = dilation
        for z_slice in range(lz):
            erosion[z_slice] = cv2.erode(img[z_slice], kernel, iterations = iterations)
        img = erosion
        i += 1;
    return img

def opening(img, kernel, iterations):
    lz, ly, lx = img.shape[0], img.shape[1], img.shape[2]
    erosion = np.zeros((lz, ly, lx))
    dilation = np.zeros((lz, ly, lx))
    
    # Opening
    # -> removes islands / small regions
    i = 0
    while i < 1:

        for z_slice in range(lz):
            erosion[z_slice] = cv2.erode(img[z_slice], kernel, iterations = iterations)
        img = erosion
        for z_slice in range(lz):
            dilation[z_slice] = cv2.dilate(img[z_slice], kernel, iterations = iterations)
        img = dilation
        i += 1;
    return img