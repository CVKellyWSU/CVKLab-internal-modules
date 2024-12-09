# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 17:59:20 2022

@author: cvkelly
"""

import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.optimize import curve_fit

def import_npz(npz_file,allow_pickle=False):
    # This doesn't work as part of a Python module
    Data = np.load(npz_file,allow_pickle=allow_pickle)
    for varName in Data:
        globals()[varName] = Data[varName]
        
def get_time(file):
    t1 = file.find('time')
    t2 = file.find('_z')
    t = int(file[(t1+4):t2])
    return t

def calc_droplet_volume(x,z):
    dz = z[1:]-z[:-1]
    meanx = (x[1:]+x[:-1])/2
    # unitvol = pi*r^2 * dz
    vol = np.pi*meanx**2*dz
    return np.sum(vol)

def analyze_ilastik(ilastik):
    contours, hierarchy = cv2.findContours(ilastik, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    if len(contours)>0:
        areas = np.zeros(len(contours))
        for ci in range(len(contours)):
            areas[ci] = cv2.contourArea(contours[ci])
        bestcont = np.argmax(areas)
        area = areas[bestcont]
        contour = contours[bestcont].reshape(-1,2)
    else:
        area = 0
        contour = 0
    return contour, area

def y_circ_plus(x,x0,y0,r):
    r = np.abs(r)
    y = np.zeros(x.shape)
    keep = np.all([x<(x0+r),x>x0-r],axis=0)
    y[keep] = np.sqrt(r**2-(x[keep]-x0)**2)+y0
    return y

def y_circ_neg(x,x0,y0,r):
    r = np.abs(r)
    y = np.zeros(x.shape)
    keep = np.all([x<(x0+r),x>x0-r],axis=0)
    y[keep] = -np.sqrt(r**2-(x[keep]-x0)**2)+y0
    return y

def fit_R0(x,y,p0=(1,0,1),bounds = ((-0.05,0,0),(0.05,5,10))):
    x2 = -y.copy()
    y2 = x.copy()
    x2 = x2-(np.max(x2)+np.min(x2))/2
    y2 = y2-np.min(y2)

    fres = curve_fit(y_circ_neg,x2,y2,p0=p0,bounds=bounds)
    return fres[0],x2,y2


