# -*- coding: utf-8 -*-
"""
Standard FCCS functions

Created on Mon Apr 16 14:32:31 2018
@author: CVKelly


Updates
Apr 21, 2021:   to load TIF files from Solis
Aug 2022:       for calibration shifting from time-average data
Sept 22, 2022:  to make plotting the fit spectra easier
May 4, 2023:    Added get zeros for file name funct
Jan 15, 2024:   Added four-color FCCS for proteins and Cy5, fit and plot many G in one function
- copied to IX83 on April 24, 2024
April 29, 2024: Added subfolder inspectioan definitions
July 1, 2024: Added analysis functions, like subfolder finding
July 22, 2024: Made calc_xcorr analyze both +/- tau and average the results, and added more oscilation options
Oct 22, 2024: Added spectrum for CFP and Cy5
Dec 9, 2024: Updated CFF and Cy5 spectrum

"""
import numpy as np
from scipy.optimize import curve_fit as fit
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from os import listdir, remove
import os
from os.path import exists, isdir
#from PIL import Image
from skimage import io
import matplotlib
import sif_parser
from scipy.signal import savgol_filter


# %% functions to find files
class dumb():
    pass

def get_subfolders(fold):
	consid = listdir(fold)
	subfolds = []
	for i,con in enumerate(consid):
		if isdir(fold+con):
			subfolds.append(con)
	return subfolds

def get_new_subfolds(fold,oldsubfolds):
	subfolds = get_subfolders(fold)
	found = np.zeros(len(subfolds))
	newsubfolds = []
	for i,sf in enumerate(subfolds):
		for j,sff in enumerate(oldsubfolds):
			if sf == sff:
				found[i] = 1
				break
		if found[i] == 0:
			newsubfolds.append(sf)
	return newsubfolds

def assess_subfolds(fold_new_data,ff_status):
    subfolds = get_subfolders(fold_new_data)   
    subfolds_done = {}
    f_status = open(ff_status,'r')
    l_status = f_status.readlines()
    f_status.close()
    newsubfolds = []
    oldsubfolds = []
    for i,subfold in enumerate(subfolds):
        found = 0
        for j,line in enumerate(l_status):
            if line.find(subfold)>-1:
                oldsubfolds.append(subfold)
                found = 1
                break
        if found == 0:
            if subfold[-1] == '\\':
                newsubfolds.append(subfold)
            else:
                newsubfolds.append(subfold+'\\')

        subfolds_done[subfold] = found
    return subfolds,subfolds_done,oldsubfolds,newsubfolds
    
	


# %% functions for display
def get_colors(num_cols,cmap='jet'):
# num_cols = 4
    var = np.linspace(0,1,num_cols)
    cols = np.zeros((num_cols,4))
    for i in range(num_cols):
        if cmap == 'jet':
            cols[i,:] = matplotlib.cm.jet(var[i])
        elif cmap == 'hsv':
            cols[i,:] = matplotlib.cm.hsv(var[i])
        elif cmap == 'gray':
            cols[i,:] = matplotlib.cm.gray(var[i])
        elif cmap == 'rainbow':
            cols[i,:] = matplotlib.cm.rainbow(var[i])
        elif cmap == 'tab20b':
            cols[i,:] = matplotlib.cm.tab20b(var[i])
        elif cmap == 'nipy_spectral':
            cols[i,:] = matplotlib.cm.nipy_spectral(var[i])
        elif cmap == 'viridis':
            cols[i,:] = matplotlib.cm.viridis(var[i])
        else:
            cols[i,:] = matplotlib.cm.jet(var[i])
            print('   ', cmap,' not found. Jet used.')
    cols = cols.astype(np.float16)
    return(cols)
    
def fname_num(num,length):
    # makes num a length
    # fname_num(55,5) returns '00055'
    zs = '0'*length
    num = str(int(num))
    nl = len(num)
    if nl < length:
        zs = zs[nl:]
    else:
        zs = ''
    return zs+num
        
def load_meta_funct(file)  :
    file = file.replace('-RawData','')
    file = file.replace('.dat','')
    file = file.replace('-Corr','')
    file = file.replace('Fit','')
    file = file.replace('Para','')
    file = file.replace('-Meta','')
    file = file+'-Meta.dat'
    a = np.loadtxt(file)
    meta = dumb()
    meta.numframes = a[0].astype(np.uint32)
    meta.realexp = a[1]
    meta.realHz = a[2]
    meta.setexp = a[3]
    meta.sensorw = a[4].astype(np.uint16)
    meta.sensorh = a[5].astype(np.uint16)
    meta.binw = a[6].astype(np.uint16)
    meta.binh = a[7].astype(np.uint16)
    meta.width1 = a[8].astype(np.uint16)
    meta.height1 = a[9].astype(np.uint16) # width for line extraction
    meta.top1 = a[10].astype(np.uint16)
    meta.bot1 = a[11].astype(np.uint16) # center line for extraction
    meta.width2 = a[12].astype(np.uint16)
    meta.height2 = a[13].astype(np.uint16)
    meta.top2 = a[14].astype(np.uint16)
    meta.bot2 = a[15].astype(np.uint16)
    meta.numTau = a[16]
    try:
        meta.temp = a[18]
        meta.maxTau = a[17]
    except:
        meta.maxTau = 1
        meta.temp = a[17]

    return meta,a

def load_metaTIF_funct(file)  :
    file = file.replace('-RawData','')
    file = file.replace('.dat','')
    file = file.replace('-Corr','')
    file = file.replace('Fit','')
    file = file.replace('Para','')
    file = file.replace('-Meta','')
    file = file+'-Meta.dat'
    a = np.loadtxt(file)
    meta = dumb()
    meta.numframes = a[0].astype(np.uint32)
    meta.realexp = a[1]
    meta.realHz = a[2]
    meta.setexp = a[3]
    meta.sensorw = a[4].astype(np.uint16)
    meta.sensorh = a[5].astype(np.uint16)
    meta.binw = a[6].astype(np.uint16)
    meta.binh = a[7].astype(np.uint16)
    meta.width1 = a[8].astype(np.uint16)
    meta.height1 = a[9].astype(np.uint16) # width for line extraction
    meta.top1 = a[10].astype(np.uint16)
    meta.bot1 = a[11].astype(np.uint16) # center line for extraction
    meta.width2 = a[12].astype(np.uint16)
    meta.height2 = a[13].astype(np.uint16)
    meta.top2 = a[14].astype(np.uint16)
    meta.bot2 = a[15].astype(np.uint16)
    meta.numTau = a[16]
    try:
        meta.temp = a[18]
        meta.maxTau = a[17]
    except:
        meta.maxTau = 1
        meta.temp = a[17]

    return meta,a


def load_FRET_meta_funct(file)  :
    file = file.replace('-RawData','')
    file = file.replace('.dat','')
    file = file.replace('-Corr','')
    file = file.replace('Fit','')
    file = file.replace('Para','')
    file = file.replace('-Meta','')
    file = file+'-Meta.dat'
    a = np.loadtxt(file)
    meta = dumb()
    meta.numframes = a[0].astype(np.uint32)
    meta.realexp = a[1]
    meta.setexp = a[2]
    meta.sensorw = a[3].astype(np.uint16)
    meta.sensorh = a[4].astype(np.uint16)
    meta.binw = a[5].astype(np.uint16)
    meta.binh = a[6].astype(np.uint16)
    meta.bot1 = a[7].astype(np.uint16) # center line for extraction
    meta.height1 = a[8].astype(np.uint16) # width for line extraction
    meta.fret_im_height = a[9].astype(np.uint16)
    meta.fret_im_width = a[10].astype(np.uint16)
    meta.fret_im_step = a[11].astype(np.uint16)
    meta.sensor_temp = a[12].astype(np.uint16)

    return meta,a


def load_data_funct(file,meta):
    file = file.replace('-RawData','')
    file = file.replace('.dat','')
    file = file.replace('-Corr','')
    file = file.replace('Fit','')
    file = file.replace('Para','')
    file = file.replace('-Meta','')
    file = file+'-RawData.dat'

    f = open(file, "r")
    a = np.fromfile(f, dtype=np.uint16,count=-1,sep='')
    data = a.reshape((meta.numframes,meta.sensorw,meta.sensorh))
    data = data.swapaxes(0, 2)
    data = data.swapaxes(0, 1)

    f.close()

    return data


def load_TIF_funct(file,isolate_pos = False):
    #file must be the full path to the 3D TIF file
    # isolate_pos == extract one locatoin from a micromanager MDA multi-position array
    
#    im = Image.open(file)
    im = io.imread(file).astype(np.float64)

    if len(im.shape)==4 and isolate_pos:
        # find pos
        pos = -1
        trynow = 0
        while pos == -1:
            if file.find('Pos'+str(trynow)+'.ome.tif')>-1:
                pos = trynow
            trynow += 1
        im = im[pos,:,:,:]


    n_frames = im.shape[0]
    hei = im.shape[1]
    wid = im.shape[2]

    data = im.transpose(1,2,0)

    frame_time = -1
    return data,n_frames,wid,hei,frame_time

def load_SIF_funct(file):
    data,info = sif_parser.np_open(file)
    n_frames = data.shape[0]
    hei = data.shape[1]
    wid = data.shape[2]
    data = data.transpose(1,2,0)

    frame_time = info['AccumulatedCycleTime']
    
    return data,n_frames,wid,hei,frame_time 



def get_lines_from_array(data,linenum=-1,width=4):

#    a = avea

    # Collapse from 2D image to 1D line at each time
    # linenum = meta.bot1 = center of data to extract
    # width = meta.height1 = width of data to extract
    
    # consider data = np.moveaxis(data,[0],[2])
    
    s0 = data.shape[0]  # data height
    s1 = data.shape[1]  # data width
    s2 = data.shape[2]  # number of frames, i.e. time
#    alldatalines = np.zeros((s0,s2)) ## for vertical lines with bad prism orientation
    alldatalines = np.zeros((s1,s2)) ## for horizontal lines with proper prism orientation
    
    
    data = data[1:,:,:]  # always remove top line that has noise from cropped sensor mode

    if s0 > 4:
        if linenum == -1:
            d1 = np.sum(data,axis=2)
            d1 = np.sum(d1,axis=1)
            if len(d1)>10:
                d1[:3] = 0
                d1[-3:] = 0
            args = np.argsort(d1)
            linenum = args[-2]
            
        if width/2 == np.round(width/2):  # even
            down = width/2-1
            up = width/2
        else:  # odd
            down = (width-1)/2
            up = (width-1)/2
        down = np.max(np.array([1,linenum-down])) ## changed 0 to 1 on 240627
        up = np.min(np.array([s0-1,linenum+up+1])) ## for vertical lines with bad prism orientation
        down=down.astype(np.uint16)
        up=up.astype(np.uint16)
        data = data[down:up,:,:]

    alldatalines[:] = np.sum(data,axis=0)

    return alldatalines

# %% get calibration peak if only 1 fluor is present
def get_calib_from_alldata_ave(alldatalines):
    # get standard fit funct from the time average line
    # used when only one color diffuser is present
    calib = np.mean(alldatalines,axis=1)
    calib = calib/np.sum(calib)
    return calib

def convert_and_save_3Ddata_as_calib(data,file_to_save,linenum=-1,width=4):
    # convert_data_to_calib(data,file_to_save)
    # data = array of (width, height, frame_num) if data.dim == 3
    # data = array of (width, height) if data.dim = 2
    #   use  data = load_TIF_funct(file) if necessary
    # linenum = center of good data
    # width = how many rows to ave togeter

    alldatalines = get_lines_from_array(data,linenum,width)
    spec = get_calib_from_alldata_ave(alldatalines)
    spec = spec  - np.min(spec)
    spec = spec/np.sum(spec)
    np.savetxt(file_to_save,spec)

    return spec

def TIF_to_calib_file(file,file_to_save='',linenum=-1,width=4):
    # file = full path to 3D .TIF file from Solis
    # file_to_save (optional) = new .npz file towrite the 1D matrix
    data,a,b,c = load_TIF_funct(file)
    if len(file_to_save) == 0:
        file_to_save = file[:file.find('.tif')]+'.dat'
    if exists(file_to_save):
        print('ERROR: FILE TO SAVE ALREADY EXISTS, line 232 in sff')
        print('   No file being saved. Process aborted.')
        return 'Error in sff.TIF_to_calib_file'

    spec = convert_and_save_3Ddata_as_calib(data,file_to_save,linenum,width)
    return spec

def resave_TIF_as_smaller_TIF(file,new_width=32):
    f = open('D:/Default FCS Data Folder/00000-FileWrite.txt')
    file_write = f.readline()
    f.close()

    data,a,b,c = load_TIF_funct(file)

    old_width = data.shape[1]
    cut = int((old_width-new_width)/2)

    if (cut < old_width/2) and (cut > 0):
        data_smaller = data[:,cut:-cut,:].astype(np.uint16)
    else:
        data_smaller = data.astype(np.uint16)

#    np.save(file[:file.find('.tif')]+'.npy',data_smaller)
#         remove(file)

    if data_smaller.shape != data.shape:
        io.call_plugin('imsave',plugin='tifffile',file=file_write,data=data_smaller.transpose(2,0,1),imagej=True, metadata={'axes': 'XYT', 'fps': 10.0})

    return data_smaller # 'successfully saved TIF as NPY'


def resave_DATA_as_smaller_TIF(data,file_write,new_width=32,linenum=-1,new_height=4):
#    f = open('D:/Default FCS Data Folder/00000-FileWrite.txt')
#    file_write = f.readline()
#    f.close()

# crop width to be new_width=32
    old_width = data.shape[1]
    cut = int((old_width-new_width)/2)
    if (cut < old_width/2) and (cut > 0):
        data_smaller = data[:,cut:-cut,:].astype(np.uint16)
    else:
        data_smaller = data.astype(np.uint16)

# crop height to be new_height=4 around brightest line
    old_height = data.shape[0]
    if new_height < old_height:
        if linenum == -1:
            d1 = np.sum(data,axis=2)
            d1 = np.sum(d1,axis=1)
            args = np.argsort(d1)
            linenum = args[-1]

        if new_height/2 == np.round(new_height/2):  # even
            down = new_height/2-1
            up = new_height/2
        else:  # odd
            down = (new_height-1)/2
            up = (new_height-1)/2
        down = np.max(np.array([1,linenum-down]))
        up = np.min(np.array([old_height,linenum+up+1]))
        down=down.astype(np.uint16)
        up=up.astype(np.uint16)

        data_smaller2 = data_smaller[down:up,:,:]

    else:
        data_smaller2 = data_smaller

    if data_smaller2.shape != data.shape:
        print('   Data has been resized and resaved. New size:',data_smaller2.shape)
#    else:
#        print('Data was not sized.')

    if data_smaller.shape != data.shape:
        io.call_plugin('imsave',plugin='tifffile',file=file_write,data=data_smaller2.transpose(2,0,1),imagej=True, metadata={'axes': 'XYT', 'fps': 10.0})

    return data_smaller #, txt


# %% get simple intensity vs time for whole image minus background
def simp_inten_vs_time(alldatalines,subave = 1,numbackpix=6):
    # subave = 0 -- super noisy at long tau
    # subave = 1 -- best
    # subave = 2 -- noisier at short tau
    # numbackpix doesn't matter much until >10
    if subave == 0: # subtract nothing
        pass
    elif subave == 1: # subtract mean background
        tavespec = np.mean(alldatalines,axis=1)
        alldatalines = alldatalines-np.mean(tavespec[-numbackpix:])
    elif subave == 2: # subtract background from each frame
        for t in range(alldatalines.shape[1]):
            alldatalines[:,t] = alldatalines[:,t]-np.mean(alldatalines[-numbackpix:,t])
    inten = np.zeros(alldatalines.shape[1])
    for t in range(alldatalines.shape[1]):
        inten[t] = np.sum(alldatalines[:,t])
    return inten


# %% find the inverse matrices to use for fitting
def get_inv_mat1(calib):
    o = np.ones(calib.shape[0])
    o = o/np.sum(o)
    m = np.zeros(4)
    m[0] = np.sum(calib**2)
    m[1] = np.sum(o*calib)
    m[2] = m[1]
    m[3]= np.sum(o**2)
    M = np.matrix(((m[0],m[1]),(m[2],m[3])))
    Mi = np.linalg.inv(M)
    return Mi

def get_inv_mat2(calib1,calib2):
    calib1 = calib1/np.sum(calib1)
    calib2 = calib2/np.sum(calib2)
    o = np.ones(calib1.shape[0])
    o = o/np.sum(o)
    M = np.matrix(((np.sum(calib1**2)    ,np.sum(calib1*calib2),np.sum(calib1*o)),
                   (np.sum(calib1*calib2),np.sum(calib2**2)    ,np.sum(calib2*o)),
                   (np.sum(calib1*o)     ,np.sum(o*calib2)     ,np.sum(o**2) )))
    Mi = np.linalg.inv(M)
    return Mi

def get_inv_mat2_noBack(calib1,calib2):
    M = np.matrix(((np.sum(calib1**2)    ,np.sum(calib1*calib2)),
                   (np.sum(calib1*calib2),np.sum(calib2**2))))
    Mi = np.linalg.inv(M)
    return Mi

def get_inv_mat3(calib1,calib2,calib3):
    o = np.ones(calib1.shape[0])
    o = o/np.sum(o)
    M = np.matrix(((np.sum(calib1**2)    ,np.sum(calib1*calib2),np.sum(calib1*calib3),np.sum(calib1*o)),
                   (np.sum(calib1*calib2),np.sum(calib2**2)    ,np.sum(calib2*calib3),np.sum(calib2*o)),
                   (np.sum(calib1*calib3),np.sum(calib2*calib3),np.sum(calib3**2),    np.sum(calib3*o)),
                   (np.sum(calib1*o)     ,np.sum(o*calib2)     ,np.sum(o*calib3),     np.sum(o**2) )))
    Mi = np.linalg.inv(M)
    return Mi

def get_inv_mat3_noBack(calib1,calib2,calib3):
    M = np.matrix(((np.sum(calib1**2)    ,np.sum(calib1*calib2),np.sum(calib1*calib3)),
                   (np.sum(calib1*calib2),np.sum(calib2**2)    ,np.sum(calib2*calib3)),
                   (np.sum(calib1*calib3),np.sum(calib2*calib3),np.sum(calib3**2)) ))
    Mi = np.linalg.inv(M)
    return Mi

def get_inv_mat4(calib1,calib2,calib3,calib4):
    o = np.ones(calib1.shape[0])
    o = o/np.sum(o)
    M = np.matrix(((np.sum(calib1**2)    ,np.sum(calib1*calib2),np.sum(calib1*calib3),np.sum(calib1*calib4),np.sum(calib1*o)),
                   (np.sum(calib1*calib2),np.sum(calib2**2)    ,np.sum(calib2*calib3),np.sum(calib2*calib4),np.sum(calib2*o)),
                   (np.sum(calib1*calib3),np.sum(calib2*calib3),np.sum(calib3**2),    np.sum(calib3*calib4),np.sum(calib3*o)),
                   (np.sum(calib1*calib4),np.sum(calib2*calib4),np.sum(calib3*calib4),np.sum(calib4**2),    np.sum(calib4*o)),
                   (np.sum(calib1*o)     ,np.sum(o*calib2)     ,np.sum(o*calib3),     np.sum(o**calib4),    np.sum(o**2) )))
    Mi = np.linalg.inv(M)
    return Mi


def get_inv_mat4_noBack(calib1,calib2,calib3,calib4):
    M = np.matrix(((np.sum(calib1**2)    ,np.sum(calib2*calib1),np.sum(calib1*calib3),np.sum(calib1*calib4)),
                   (np.sum(calib1*calib2),np.sum(calib2**2)    ,np.sum(calib2*calib3),np.sum(calib2*calib4)),
                   (np.sum(calib1*calib3),np.sum(calib2*calib3),np.sum(calib3**2),    np.sum(calib3*calib4)),
                   (np.sum(calib1*calib4),np.sum(calib2*calib4),np.sum(calib4*calib3),np.sum(calib4**2)   )))
    Mi = np.linalg.inv(M)
    return Mi

def get_inv_matN_noBack(calib_list):
    if len(calib_list) == 2:
        Mi = get_inv_mat2_noBack(calib_list[0],calib_list[1])
    elif len(calib_list) == 3:
        Mi = get_inv_mat3_noBack(calib_list[0],calib_list[1],calib_list[2])
    elif len(calib_list) == 4:
        Mi = get_inv_mat4_noBack(calib_list[0],calib_list[1],calib_list[2],calib_list[3])
    else:
        Mi = -1
    return Mi
# %% fit a single I vs X to find the contribution of each peak
def fit_one_time1(calib,Mi,data1):
    o = np.ones(calib.shape[0])
    o = o/np.sum(o)
    v1 = np.sum(calib*data1)
    v2 = np.sum(o*data1)
    V = np.matrix((v1,v2))
    V = V.transpose()
    inten = Mi*V
    return inten[0],inten[1]

def fit_one_time1_noback(data1):   # don't need this in a loop, just do the mean
    inten = np.mean(data1,axis=0)
    return inten[0],inten[1]

def fit_one_time2(calib1,calib2,Mi,data1):
    calib1 = calib1/np.sum(calib1)
    calib2 = calib2/np.sum(calib2)
    o = np.ones(calib1.shape[0])
    o = o/np.sum(o)
    v1 = np.sum(calib1*data1)
    v2 = np.sum(calib2*data1)
    v3 = np.sum(o*data1)
    V = np.matrix((v1,v2,v3))
    V = V.transpose()
    inten = Mi*V
    return inten[0][0,0],inten[1][0,0],inten[2][0,0]

def fit_one_time2_noBack(calib1,calib2,Mi,data1):
    v1 = np.sum(calib1*data1)
    v2 = np.sum(calib2*data1)
    V = np.matrix((v1,v2))
    V = V.transpose()
    inten = Mi*V
    return inten[0][0,0],inten[1][0,0]

def fit_one_time3(calib1,calib2,calib3,Mi,data1):
    o = np.ones(calib1.shape[0])
    o = o/np.sum(o)
    v1 = np.sum(calib1*data1)
    v2 = np.sum(calib2*data1)
    v3 = np.sum(calib3*data1)
    v4 = np.sum(o*data1)
    V = np.matrix((v1,v2,v3,v4))
    V = V.transpose()
    inten = Mi*V
    return inten[0][0,0],inten[1][0,0],inten[2][0,0],inten[3][0,0]

def fit_one_time3_noBack(calib1,calib2,calib3,Mi,data1):
    v1 = np.sum(calib1*data1)
    v2 = np.sum(calib2*data1)
    v3 = np.sum(calib3*data1)
    V = np.matrix((v1,v2,v3))
    V = V.transpose()
    inten = Mi*V
    return inten[0][0,0],inten[1][0,0],inten[2][0,0]

def fit_one_time4(calib1,calib2,calib3,calib4,Mi,data1):
    o = np.ones(calib1.shape[0])
    o = o/np.sum(o)
    v1 = np.sum(calib1*data1)
    v2 = np.sum(calib2*data1)
    v3 = np.sum(calib3*data1)
    v4 = np.sum(calib4*data1)
    v5 = np.sum(o*data1)
    V = np.matrix((v1,v2,v3,v4,v5))
    V = V.transpose()
    inten = Mi*V
    return inten[0][0,0],inten[1][0,0],inten[2][0,0],inten[3][0,0],inten[4][0,0]

def fit_one_time4_noBack(calib1,calib2,calib3,calib4,Mi,data1):
    o = np.ones(calib1.shape[0])
    o = o/np.sum(o)
    v1 = np.sum(calib1*data1)
    v2 = np.sum(calib2*data1)
    v3 = np.sum(calib3*data1)
    v4 = np.sum(calib4*data1)
    V = np.matrix((v1,v2,v3,v4))
    V = V.transpose()
    inten = Mi*V
    return inten[0][0,0],inten[1][0,0],inten[2][0,0],inten[3][0,0]

# %% repeat the fitting for all times in
def fit_many_times1(calib,Mi,alldatalines):
    s1 = alldatalines.shape[1] # number of times
    inten0 = np.zeros(s1)
    inten1 = np.zeros(s1)
    for i in range(s1):
        data1 = alldatalines[:,i]
        inten0[i],inten1[i] = fit_one_time1(calib,Mi,data1)
    return inten0,inten1

def fit_many_times1_noback(alldatalines):
    inten0 = np.sum(alldatalines,axis=0)
    return inten0

def fit_many_times2(calib1,calib2,Mi,alldatalines):
    if len(alldatalines.shape)==2:
        s1 = alldatalines.shape[1] # number of times
        inten0 = np.zeros(s1)
        inten1 = np.zeros(s1)
        inten2 = np.zeros(s1)
        for i in range(s1):
            data1 = alldatalines[:,i]
            inten0[i],inten1[i],inten2[i] = fit_one_time2(calib1,calib2,Mi,data1)
    else:
        inten0,inten1,inten2 = fit_one_time2(calib1,calib2,Mi,alldatalines)
    return inten0,inten1,inten2

def fit_many_times2_noBack(calib1,calib2,Mi,alldatalines):
    if len(alldatalines.shape)==2:
        s1 = alldatalines.shape[1] # number of times
        inten0 = np.zeros(s1)
        inten1 = np.zeros(s1)
        for i in range(s1):
            data1 = alldatalines[:,i]
            inten0[i],inten1[i] = fit_one_time2_noBack(calib1,calib2,Mi,data1)
    else:
        inten0,inten1 = fit_one_time2_noBack(calib1,calib2,Mi,alldatalines)
    return inten0,inten1

def fit_many_times3(calib1,calib2,calib3,Mi,alldatalines):
    s1 = alldatalines.shape[1] # number of times
    inten0 = np.zeros(s1)
    inten1 = np.zeros(s1)
    inten2 = np.zeros(s1)
    inten3 = np.zeros(s1)
    for i in range(s1):
        data1 = alldatalines[:,i]
        inten0[i],inten1[i],inten2[i],inten3[i] = fit_one_time3(calib1,calib2,calib3,Mi,data1)
    return inten0,inten1,inten2,inten3

def fit_many_times3_noBack(calib1,calib2,calib3,Mi,alldatalines):
    if len(alldatalines.shape)>1:
        s1 = alldatalines.shape[1] # number of times
    else:
        s1 = 1
    inten0 = np.zeros(s1)
    inten1 = np.zeros(s1)
    inten2 = np.zeros(s1)
    for i in range(s1):
        if len(alldatalines.shape)>1:
            data1 = alldatalines[:,i]
        else:
            data1 = alldatalines.copy()
        inten0[i],inten1[i],inten2[i] = fit_one_time3_noBack(calib1,calib2,calib3,Mi,data1)
    return inten0,inten1,inten2

def fit_many_times4(calib1,calib2,calib3,calib4,Mi,alldatalines):
    if len(alldatalines.shape)>1:
        s1 = alldatalines.shape[1] # number of times
    else:
        s1 = 1
    inten0 = np.zeros(s1)
    inten1 = np.zeros(s1)
    inten2 = np.zeros(s1)
    inten3 = np.zeros(s1)
    inten4 = np.zeros(s1)
    for i in range(s1):
        if len(alldatalines.shape)>1:
            data1 = alldatalines[:,i]
        else:
            data1 = alldatalines.copy()
        inten0[i],inten1[i],inten2[i],inten3[i],inten4[i] = fit_one_time4(calib1,calib2,calib3,calib4,Mi,data1)
    return inten0,inten1,inten2,inten3,inten4

def fit_many_times4_noBack(calib1,calib2,calib3,calib4,Mi,alldatalines):
    if len(alldatalines.shape)>1:
        s1 = alldatalines.shape[1] # number of times
        inten0 = np.zeros(s1)
        inten1 = np.zeros(s1)
        inten2 = np.zeros(s1)
        inten3 = np.zeros(s1)
        for i in range(s1):
            data1 = alldatalines[:,i]
            inten0[i],inten1[i],inten2[i],inten3[i] = fit_one_time4_noBack(calib1,calib2,calib3,calib4,Mi,data1)
    else:
        inten0,inten1,inten2,inten3 = fit_one_time4_noBack(calib1,calib2,calib3,calib4,Mi,alldatalines)

    return inten0,inten1,inten2,inten3

def fit_many_timesN_noBack(calib_list,Mi,alldatalines):
    if len(calib_list) == 1:
        i0 = fit_many_times1_noback(alldatalines)
        i1,i2,i3 = [0],[0],[0]
    elif len(calib_list) == 2:
        i0,i1 = fit_many_times2_noBack(calib_list[0],calib_list[1],Mi,alldatalines)
        i2,i3 = [0],[0]
    elif len(calib_list) == 3:
        i0,i1,i2 = fit_many_times3_noBack(calib_list[0],calib_list[1],calib_list[2],Mi,alldatalines)
        i3 = [0]
    elif len(calib_list) == 4:
        i0,i1,i2,i3 = fit_many_times4_noBack(calib_list[0],calib_list[1],calib_list[2],calib_list[3],Mi,alldatalines)
    else:
        inten_list = -1
    
        
    i0 = i0[:-1]   
    if len(i1)>1:
        i1 = i1[:-1]
    if len(i2)>1:
        i2 = i2[:-1]
    if len(i3)>1:
        i3 = i3[:-1]
        
    return i0,i1,i2,i3

# %% do correlation analysis
# def calc_xcorr(inten0,inten1,time_per_frame,maxt,numdt):
#     dt = np.logspace(0,np.log10(maxt/time_per_frame),numdt)
#     dt = dt.astype(np.uint32)
#     dt = np.unique(dt)
#     dt = dt[dt<(inten0.shape[0]/3)]
#     G = np.zeros(dt.shape[0])
#     for i in range(dt.shape[0]):
#         G[i] = np.mean(inten0[:(len(inten0)-dt[i])]*inten1[dt[i]:])
#     G = G/np.mean(inten0)/np.mean(inten1)-1
#     return dt,G

def calc_xcorr(inten0,inten1,time_per_frame,maxt,numdt):
    dt = np.logspace(0,np.log10(maxt/time_per_frame),numdt)
    dt = dt.astype(np.uint32)
    dt = np.unique(dt)
    dt = dt[dt<(inten0.shape[0]/3)]
    G1 = np.zeros(dt.shape[0])
    G2 = np.zeros(dt.shape[0])
    for i in range(dt.shape[0]):
        G1[i] = np.mean(inten0[:(len(inten0)-dt[i])]*inten1[dt[i]:])
        G2[i] = np.mean(inten1[:(len(inten1)-dt[i])]*inten0[dt[i]:])
    G1 = G1/np.mean(inten0)/np.mean(inten1)-1
    G2 = G2/np.mean(inten0)/np.mean(inten1)-1
    G = (G1+G2)/2
    return dt,G

def calc_xcorr_G0(inten0,inten1):
    G = np.zeros(2)
    for dt in [1,2]:
        G[dt-1] = np.mean(inten0[:(inten0.shape[0]-dt)]*inten1[dt:])
    G = G/np.mean(inten0)/np.mean(inten1)-1
    G0 = np.mean(G)
    return G0

def Corr2D(i1,i2,i3,time_per_frame,max_tau,num_tau):
    tau_list  = np.logspace(np.log10(10),np.log10(max_tau/time_per_frame),num_tau)
    tau_list = tau_list.astype(np.uint64)
    tau_list = np.unique(tau_list)

    G = np.zeros((len(tau_list),len(tau_list)))
    for i,tau1 in enumerate(tau_list):
        for j,tau2 in enumerate(tau_list):
            numtimes = int(len(i1)-np.max([tau1,tau2]))-1
            G[i,j] = np.mean(i1[:numtimes]*i2[tau1:int(numtimes+tau1)]*i3[tau2:int(numtimes+tau2)])
    G = G/np.mean(i1)/np.mean(i2)/np.mean(i3)-1 #normalize G
    return G,tau_list

def Corr3D(i1,i2,i3,i4,time_per_frame,max_tau,num_tau):
    tau_list  = np.logspace(np.log10(10),np.log10(max_tau/time_per_frame),num_tau)
    tau_list = tau_list.astype(np.uint64)
    tau_list = np.unique(tau_list)

    G = np.zeros((len(tau_list),len(tau_list),len(tau_list)))
    for i,tau1 in enumerate(tau_list):
        for j,tau2 in enumerate(tau_list):
            for k,tau3 in enumerate(tau_list):
                numtimes = int(len(i1)-np.max([tau1,tau2,tau3]))-1
                G[i,j,k] = np.mean(i1[:numtimes]*i2[tau1:int(numtimes+tau1)]*i3[tau2:int(numtimes+tau2)]*i4[tau3:int(numtimes+tau3)])
    G = G/np.mean(i1)/np.mean(i2)/np.mean(i3)/np.mean(i4)-1 #normalize G
    return G,tau_list

# def fit_xcorr(G,dtlistsec,model):
#     ## Model:
#     #    '2d'     = Brownian 2D model
#     #    '2danom' = Anomalous 2D model
#     #    '3d'     = Brownian 3D model

#     if (model.find('2danom') > -1 ):
#         f = lambda x,td,a,alpha,c: a/(1+(x/td)**alpha)+c
#         numparam = 3
#         start = (0.005,1,0.9,0.0)
#         bound = ((0.0,0.0,0.0,0.0),(1.0,np.inf,np.inf,np.inf))
#     elif (model.find('2d') > -1 ):
#         f = lambda x,td,a,c: a/(1+(x/td))+c
#         numparam = 3
#         start = (0.005,1,0)
#         bound = ((0.0,-np.inf,-np.inf),(1.0,np.inf,np.inf))##sonali changed lower bound from 0.0 to -np.inf on 11/22/21

#     elif (model.find('3dx') > -1):
#         f = lambda x,td,td2,a,a2,r: a/(1+(x/td))/(1+r**-2*np.sqrt(x/td)) + a2/(1+(x/td2))/(1+r**-2*np.sqrt(x/td2))
#         numparam = 3
#         start = (0.005,0.01,0.1,0.5,0.5)
#         bound = ((0.0,0.0,0.0,0.0,0.0),(1.0,np.inf,np.inf,np.inf,np.inf))
#     #elif (model.find('3d') > -1):
#         #f = lambda x,td,a,r: a/(1+(x/td))/(1+r**-2*np.sqrt(x/td))
#         #numparam = 3
#         #start = (0.005,1,1)
#         #bound = ((0.0,-np.inf,0.01),(1.0,np.inf,np.inf))#sonali changed lower bound from 0.0 to -np.inf on 11/22/21
#     elif (model.find('3d') > -1):
#         f = lambda x,td,a,r: a/(1+(x/td))/np.sqrt(1+(r**-2)*(x/td))
#         numparam = 3
#         start = (0.005,1,1)
#         bound = ((0.0,-np.inf,0.01),(1.0,np.inf,np.inf))

#     else:
#         print('no fitting function used')
#         f = lambda x,a,b: a/x+b
#         numparam = 2
#         start = (1.0,1.0)
#         bound = ((0.0,0.0),(np.inf,np.inf))

#     fres,covf = fit(f,dtlistsec,G,p0=start,bounds=bound)
# #    print(fres)
# #    param = np.zeros(numparam)
# #    for i in range(numparam):
# #        param[i] = fres[i][0]

#     return f, fres, covf, numparam  #, param

def fit_xcorr(G,dtlistsec,model,
              guess_taud = 0.01, 
              guess_G0 = 0.05,
              guess_lam = 1.5, 
              alpha=1,
              guess_oscamp = 0.05,
              guess_r = 1):
    
    
    ## Model:
    #    '2d'     = Brownian 2D model
    #    '2dosc'  = Brownian 2D model with oscilations of period "guess_lam" and fixed "alpha"
    #    '2da2'   = Brownian 2D model with fixed alpha = "alpha"
    #    '2danom' = Anomalous 2D model
    #    '3d'     = Brownian 3D model

    if guess_G0 == -1:
        guess_G0 = np.mean(G[:3])

    if (model.find('2danom') > -1 ):
        f = lambda x,td,a,alpha,c: a/(1+(x/td)**alpha)+c
        numparam = 3
        start = (guess_taud,guess_G0,0.9,0.0)
        bound = ((0.0,0.0,0.0,0.0),(1.0,np.inf,np.inf,np.inf))
    elif (model.find('2dOsc') > -1):
        f = lambda x,td,a,b,lam,c: a/(1+(x/td)**alpha)+b*np.cos(x*2*np.pi/lam)+c
        numparam = 5
        start = (guess_taud,guess_G0,guess_oscamp,guess_lam,0)
        bound = ((0.0,0.0,0.0,0.5*guess_lam,-np.inf),(1.0,np.inf,np.inf,2*guess_lam,np.inf))##sonali changed lower bound from 0.0 to -np.inf on 11/22/21
    elif (model.find('2dFixedOsc') > -1):
        f = lambda x,td,a,c: a/(1+(x/td)**alpha)+guess_oscamp*np.cos(x*2*np.pi/guess_lam)+c
        numparam = 3
        start = (guess_taud,guess_G0,0)
        bound = ((0.0,0.0,-np.inf),(1.0,np.inf,np.inf))##sonali changed lower bound from 0.0 to -np.inf on 11/22/21
    elif (model.find('2da2') > -1 ):
        f = lambda x,td,a,c: a/(1+(x/td)**alpha)+c
        numparam = 3
        start = (guess_taud,guess_G0,0)
        bound = ((0.0,-np.inf,-np.inf),(1.0,np.inf,np.inf))##sonali changed lower bound from 0.0 to -np.inf on 11/22/21
    elif (model.find('2dnoc') > -1 ):
        f = lambda x,td,a: a/(1+(x/td))
        numparam = 2
        start = (guess_taud,1)
        bound = ((0.0,-np.inf),(1.0,np.inf))##sonali changed lower bound from 0.0 to -np.inf on 11/22/21
    elif (model.find('2d') > -1 ):
        f = lambda x,td,a,c: a/(1+(x/td))+c
        numparam = 3
        start = (guess_taud,guess_G0,0)
        bound = ((0.0,-np.inf,-np.inf),(1.0,np.inf,np.inf))##sonali changed lower bound from 0.0 to -np.inf on 11/22/21
    elif (model.find('3dx') > -1):
        f = lambda x,td,td2,a,a2,r: a/(1+(x/td))/(1+r**-2*np.sqrt(x/td)) + a2/(1+(x/td2))/(1+r**-2*np.sqrt(x/td2))
        numparam = 3
        start = (guess_taud*2,guess_taud/5,guess_G0,guess_G0,0.5)
        bound = ((0.0,0.0,0.0,0.0,0.0),(1.0,np.inf,np.inf,np.inf,np.inf))
    #elif (model.find('3d') > -1):
        #f = lambda x,td,a,r: a/(1+(x/td))/(1+r**-2*np.sqrt(x/td))
        #numparam = 3
        #start = (0.005,1,1)
        #bound = ((0.0,-np.inf,0.01),(1.0,np.inf,np.inf))#sonali changed lower bound from 0.0 to -np.inf on 11/22/21
    elif (model.find('3d') > -1):
        f = lambda x,td,a,r: a/(1+(x/td))/np.sqrt(1+(r**-2)*(x/td))
        numparam = 3
        start = (guess_taud,guess_G0,guess_r)
        bound = ((0.0,-np.inf,0.01),(10.0,np.inf,np.inf))

    else:
        print('no fitting function used')
        f = lambda x,a,b: a/x+b
        numparam = 2
        start = (1.0,1.0)
        bound = ((0.0,0.0),(np.inf,np.inf))

    fres,covf = fit(f,dtlistsec,G,p0=start,bounds=bound)


    return f, fres, covf, numparam  #, param



def plot_inten0(inten0,time_per_frame,actuallyplot):
    t = np.linspace(0,inten0.shape[0]*time_per_frame,inten0.shape[0])
    N = 20  # make a round number
#    smoothedI = np.convolve(inten0, np.ones((N,))/N, mode='valid')
#    offset = np.array([N/2],dtype=np.int16)
#    smoothedt = t[offset[0]:-(offset[0]-1)]
    if actuallyplot:
        plt.plot(t,inten0,label='intensity')
#        plt.plot(smoothedt,smoothedI,label='smoothed')
        plt.xlabel('Time (sec)')
        plt.ylabel('Intensity')
        plt.title('Single I vs. t')
        plt.xlim([0,np.max(t)])
        plt.legend()
        plt.show()

    return t, inten0  # , smoothedt, smoothedI

def plot_xcorr(G,dtlistsec,actuallyplot):
    xfit = np.logspace(np.log10(np.min(dtlistsec)),np.log10(np.max(dtlistsec)),1000)
    yfit = np.zeros(xfit.shape[0])
    if actuallyplot:
        plt.semilogx(dtlistsec,G,'o',label='data')
        plt.xlabel('tau (sec)')
        plt.ylabel('G(tau)')
#        plt.title('Autocorrelation')
        plt.show()

    return dtlistsec, G, xfit, yfit

def plot_xcorr_and_fit(G,dtlistsec,f,fres,actuallyplot):
    xfit = np.logspace(np.log10(np.min(dtlistsec)),np.log10(np.max(dtlistsec)),1000)
    yfit = xfit*0
    for i in range(xfit.shape[0]):
        yfit[i] = f(xfit[i],*fres)
    if actuallyplot:
        plt.semilogx(dtlistsec,G,'o',label='data')
        plt.semilogx(xfit,yfit,'-',label='fit')
        plt.xlabel('tau (sec)')
        plt.ylabel('G(tau)')
#        plt.title('Correlation')
        plt.show()

    return dtlistsec, G, xfit, yfit

def plot_xcorr_and_fit_all(Glist,dtlistsec,model='2d',actuallyplot=True,title = ''):
    col = get_colors(len(Glist))
    for autoi,G in enumerate(Glist):
        xfit = np.logspace(np.log10(np.min(dtlistsec)),np.log10(np.max(dtlistsec)),1000)
        try:
            f, fres, covf, numparam = fit_xcorr(G,dtlistsec,model)
            yfit = f(xfit,*fres)
            tit = 'fit '+str(autoi)
        except:
            tit = 'fit '+str(autoi) + ' failed'
        if actuallyplot:
            plt.semilogx(dtlistsec,G,'o',color = col[autoi],label='data '+str(autoi))
            plt.semilogx(xfit,yfit,'-',color = col[autoi],label=tit)

    plt.xlabel('tau (sec)')
    plt.ylabel('G(tau)')
    plt.legend()
    plt.title(title)
    plt.show()

    return dtlistsec, G, xfit, yfit



# %% assorted functions for analysis
# def make_file_string(num):
    # z = '00000'
    # if type('sadf') == type(num):
        # num = int(num)
    # if not (type('sadf') == type(num)):
        # num = np.array([num])
    # out = []
    # for i in range(num.shape[0]):
        # n = str(int(num[i]))
        # zn = z[:(5-len(n))]
        # out.append(zn+n)
    # if len(out) == 1:
        # out = out[0]
    # return out

def make_file_string(num):
    z = '00000'
    if (type(num)==type(3)) | (type(num)==type(3.3)):
        num = np.array([num])
    out = []
    for i in range(len(num)):
        n = str(int(num[i]))
        zn = z[:(5-len(n))]
        out.append('\\'+zn+n)
    return out


def get_filename_from_num(fold,num):
    fnumstr = num
    if type(num) != type('asdf'):
        fnumstr = make_file_string(num)
    files_in_fold = listdir(fold)
    for file in files_in_fold:
        if file.startswith(fnumstr) and file.endswith('RawData.dat'):
            file = file[:-12]
            break
    return file


# added 2022-08-11
def shift_calib(shifts,calib1,calib2,aveavelinedata):

    shift = shifts.copy()
    maxx = len(aveavelinedata)
    calib1 = calib1[:maxx]
    calib2 = calib2[:maxx]
    calib1 = calib1/np.sum(calib1)
    calib2 = calib2/np.sum(calib2)
    x = np.arange(maxx)

    f = lambda x,a1,a2,x1,x2,s1,s2,c : a1*np.exp(-(x-x1)**2/2/s1**2)+ \
                                                a2*np.exp(-(x-x2)**2/2/s2**2)+ \
                                                +c
    p0 = (0.3,0.3,9,9,1,1,0)
    fres1 = curve_fit(f,x,calib1,p0=p0,bounds=((0,0,0,0,0,0,0),(1000,1000,100,100,10,10,100000)))
    fres2 = curve_fit(f,x,calib2,p0=p0,bounds=((0,0,0,0,0,0,0),(1000,1000,100,100,10,10,100000)))

    a11,a21,x11,x21,s11,s21,c1 = fres1[0] # be careful of copy vs memory blah
    f1 = lambda x,a,x0 : a*(a11*np.exp(-(x-x11-x0)**2/2/s11**2)+ \
                            a21*np.exp(-(x-x21-x0)**2/2/s21**2))+ \
                            +c1

    a12,a22,x12,x22,s12,s22,c2 = fres2[0] # be careful of copy vs memory blah
    f2 = lambda x,a,x0 : a*(a12*np.exp(-(x-x12-x0)**2/2/s12**2)+ \
                            a22*np.exp(-(x-x22-x0)**2/2/s22**2))+ \
                            +c2

    # % get best shift fitting
    f = lambda x,a1,a2,x1,x2,c : f1(x,a1,x1) + f2(x,a2,x2) + c
    p0=(20,100,0,0,np.mean(aveavelinedata[-3:]))
    fres = curve_fit(f,x,aveavelinedata,p0=p0)
    a13,a23,x13,x23,c3 = fres[0] #  be careful of copy vs memory blah

    # % apply best AND custom shifting for non-linear LS fitting:

    calib1 = f1(x,1,x13+shift[0])
    calib2 = f2(x,1,x23+shift[1])
    calib1 = calib1/np.sum(calib1)
    calib2 = calib2/np.sum(calib2)
    calib1 = calib1-np.mean(calib1[-2:])
    calib2 = calib2-np.mean(calib2[-2:])

    Mi = get_inv_mat2(calib1,calib2)
    aveavelinedata = aveavelinedata-np.mean(aveavelinedata[-2:])
    i1,i2,ic = fit_many_times2(calib1,calib2,Mi,aveavelinedata)

    return calib1,calib2,i1,i2


# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # #
# %% functions for getting the proper calibration shifts
# This is a new way of doing things with fixed calibration spectra.

def get_cal_functs_SimpleAuto(x00 = 0):
    # analyze the sum intensity vs time for the whole image without fitting
    fcal1 = lambda x,a0,x0 : a0+(x+x0+x00)*0
    cal_functs = [fcal1]
    return cal_functs

def get_cal_functs_4beads(x00 = 0.5):
    # The parameters in this function are found from minimizing the crosscorrelations to the four-color beads experiments that should have had no cross-correlations. The individual continuous calibration functions were fit to the isolated samples of just one color, then the continuous functions were shifted to fit the four-bead samples. The shifts between the four calibrations were fixed and the four were only allowed to shift together for the future samples of unknown cross-correlations with the function "get_calib_sum_funct"

    # output: cal_functs is a list of the four continous functions for each of the four beads. Inputs to these functions are just amplitude and a shared shift.

    cal_fres = np.array([[ 1.01755012e+00, -5.95791167e-03,  4.77725493e+00, 2.37676233e+01,  9.93809631e-01,  7.11386183e+00, 5.06077394e-03],
                       [ 6.07318378e-01,  4.36963355e-01,  8.10367394e+00, 7.15765339e+00,  9.52707668e-01,  2.16991470e+00, 3.90770722e-04],
                       [ 3.87106582e-01,  9.64674964e-01,  8.18338601e+00, 1.19257260e+01,  1.64262927e+00,  1.12007239e+00, 3.76959928e-03],
                       [ 1.05083661e-01,  9.88523163e-01,  1.23750776e+01, 1.51297056e+01,  1.67860097e+00,  1.12225143e+00, 9.99210771e-04]])
    cal_funct = []
    f2 = lambda x,a,x0,a1,a2,x1,x2,s1,s2,c : a*(a1*np.exp(-(x-x1-x0)**2/2/s1**2)+ \
                                                 a2*np.exp(-(x-x2-x0)**2/2/s2**2)+ \
                                                 +c)
    fcal1 = lambda x,a,x0 : f2(x,a,x00+x0-0.35746571,*cal_fres[0])
    fcal2 = lambda x,a,x0 : f2(x,a,x00+x0-0.31373466,*cal_fres[1])
    fcal3 = lambda x,a,x0 : f2(x,a,x00+x0-0.54969673,*cal_fres[2])
    fcal4 = lambda x,a,x0 : f2(x,a,x00+x0-0.488479,*cal_fres[3])
    cal_functs = [fcal1,fcal2,fcal3,fcal4]
    return cal_functs

def get_cal_functs_CFPYFPmCherryCy5(x00 =0.5):
    cal_fres = [np.array([ 3.52614527e+01, -3.50497503e+01,  1.05468370e+00,  1.09445042e+01, 1.09351876e+01,  1.55977067e+01,  1.36243273e+00,  1.35352638e+00, 1.21936334e+00, -4.67242434e-05]),
                 np.array([ 5.75069912e-01, -2.88585262e-01,  6.84275447e-01,  1.00779767e+01, 9.59890996e+00,  1.18144094e+01,  1.86662362e+00,  9.35915820e-01, 9.99607785e-01,  8.82043169e-05]),
                 np.array([2.43643187e-01, 9.11137379e-01, 1.32761033e-01, 4.83673726e+00, 7.63223497e+00, 7.85430786e+00, 1.00616406e+00, 1.13885975e+00, 1.73108170e+00, 2.12601846e-04]),
                 np.array([4.95870095e-01, 6.34787342e-01, 1.00885136e-11, 5.88566435e+00, 4.88221889e+00, 4.96833983e-02, 1.17322429e+00, 9.54524032e-01, 0.80000000e+00, 5.00695133e-03])]
    cal_funct = []
    g3 = lambda x,a,x0,a1,a2,a3,x1,x2,x3,s1,s2,s3,c : a*(a1*np.exp(-(x-x1-x0)**2/2/s1**2)+ \
                                                         a2*np.exp(-(x-x2-x0)**2/2/s2**2)+ \
                                                         a3*np.exp(-(x-x3-x0)**2/2/s3**2)+ \
                                                         +c)
    fcal1 = lambda x,a,x0 : g3(x,a,x0+x00-0.206,*cal_fres[0])
    fcal2 = lambda x,a,x0 : g3(x,a,x0+x00+0.264,*cal_fres[1])
    fcal3 = lambda x,a,x0 : g3(x,a,x0+x00+0.683,*cal_fres[2])
    fcal4 = lambda x,a,x0 : g3(x,a,x0+x00,*cal_fres[3])
    cal_functs = [fcal1,fcal2,fcal3,fcal4]
    return cal_functs
    
def get_cal_functs_CFPCy5(x00 =0):
    cal_fres = [np.array([ 3.53814527e+01, -3.50497503e+01,  1.05468370e+00,  1.09445042e+01, 1.09351876e+01,  1.55977067e+01,  1.36243273e+00,  1.35352638e+00, 1.21936334e+00, -4.67242434e-05]),
                 np.array([4.95870095e-01, 6.34787342e-01, 1.00885136e-11, 5.88566435e+00, 4.88221889e+00, 4.96833983e-02, 1.17322429e+00, 9.54524032e-01, 0.80000000e+00, 5.00695133e-03])]
    cal_funct = []
    g3 = lambda x,a,x0,a1,a2,a3,x1,x2,x3,s1,s2,s3,c : a*(a1*np.exp(-(x-x1-x0)**2/2/s1**2)+ \
                                                         a2*np.exp(-(x-x2-x0)**2/2/s2**2)+ \
                                                         a3*np.exp(-(x-x3-x0)**2/2/s3**2)+ \
                                                         +c)
    fcal1 = lambda x,a,x0 : g3(x,a,x0+x00-3.16,*cal_fres[0])
    fcal2 = lambda x,a,x0 : g3(x,a,x0+x00-1.68,*cal_fres[1])
    print('6')
    cal_functs = [fcal1,fcal2]
    return cal_functs
    

def get_cal_functs_CFPYFPmCherry(x00 =0.5):
    cal_fres = [np.array([ 3.52614527e+01, -3.50497503e+01,  1.05468370e+00,  1.09445042e+01, 1.09351876e+01,  1.55977067e+01,  1.36243273e+00,  1.35352638e+00, 1.21936334e+00, -4.67242434e-05]),
                 np.array([ 5.75069912e-01, -2.88585262e-01,  6.84275447e-01,  1.00779767e+01, 9.59890996e+00,  1.18144094e+01,  1.86662362e+00,  9.35915820e-01, 9.99607785e-01,  8.82043169e-05]),
                 np.array([2.43643187e-01, 9.11137379e-01, 1.32761033e-01, 4.83673726e+00, 7.63223497e+00, 7.85430786e+00, 1.00616406e+00, 1.13885975e+00, 1.73108170e+00, 2.12601846e-04])]
    cal_funct = []
    g3 = lambda x,a,x0,a1,a2,a3,x1,x2,x3,s1,s2,s3,c : a*(a1*np.exp(-(x-x1-x0)**2/2/s1**2)+ \
                                                         a2*np.exp(-(x-x2-x0)**2/2/s2**2)+ \
                                                         a3*np.exp(-(x-x3-x0)**2/2/s3**2)+ \
                                                         +c)
    fcal1 = lambda x,a,x0 : g3(x,a,x00+x0-0.206,*cal_fres[0])
    fcal2 = lambda x,a,x0 : g3(x,a,x00+x0+0.264,*cal_fres[1])
    fcal3 = lambda x,a,x0 : g3(x,a,x00+x0+0.683,*cal_fres[2])
    cal_functs = [fcal1,fcal2,fcal3]
    return cal_functs

def get_cal_functs_CFPYFP(x00 =0.5):
    cal_fres = [np.array([ 3.52614527e+01, -3.50497503e+01,  1.05468370e+00,  1.09445042e+01, 1.09351876e+01,  1.55977067e+01,  1.36243273e+00,  1.35352638e+00, 1.21936334e+00, -4.67242434e-05]),
                 np.array([ 5.75069912e-01, -2.88585262e-01,  6.84275447e-01,  1.00779767e+01, 9.59890996e+00,  1.18144094e+01,  1.86662362e+00,  9.35915820e-01, 9.99607785e-01,  8.82043169e-05])]
    g3 = lambda x,a,x0,a1,a2,a3,x1,x2,x3,s1,s2,s3,c : a*(a1*np.exp(-(x-x1-x0)**2/2/s1**2)+ \
                                                         a2*np.exp(-(x-x2-x0)**2/2/s2**2)+ \
                                                         a3*np.exp(-(x-x3-x0)**2/2/s3**2)+ \
                                                         +c)
    fcal1 = lambda x,a,x0 : g3(x,a,x0+x00-0.206,*cal_fres[0])
    fcal2 = lambda x,a,x0 : g3(x,a,x0+x00+0.264,*cal_fres[1])
    cal_functs = [fcal1,fcal2]
    return cal_functs

def get_cal_functs_CFPmCherry(x00 =0.5):
    cal_fres = [np.array([ 3.52614527e+01, -3.50497503e+01,  1.05468370e+00,  1.09445042e+01, 1.09351876e+01,  1.55977067e+01,  1.36243273e+00,  1.35352638e+00, 1.21936334e+00, -4.67242434e-05]),
                np.array([2.43643187e-01, 9.11137379e-01, 1.32761033e-01, 4.83673726e+00, 7.63223497e+00, 7.85430786e+00, 1.00616406e+00, 1.13885975e+00, 1.73108170e+00, 2.12601846e-04])]
    g3 = lambda x,a,x0,a1,a2,a3,x1,x2,x3,s1,s2,s3,c : a*(a1*np.exp(-(x-x1-x0)**2/2/s1**2)+ \
                                                         a2*np.exp(-(x-x2-x0)**2/2/s2**2)+ \
                                                         a3*np.exp(-(x-x3-x0)**2/2/s3**2)+ \
                                                         +c)
    fcal1 = lambda x,a,x0 : g3(x,a,x0-0.206,*cal_fres[0])
    fcal2 = lambda x,a,x0 : g3(x,a,x0+0.683,*cal_fres[1])
    cal_functs = [fcal1,fcal2]
    return cal_functs

def get_cal_functs_YFPmCherry(x00 =0.5):
    cal_fres = [np.array([ 5.75069912e-01, -2.88585262e-01,  6.84275447e-01,  1.00779767e+01, 9.59890996e+00,  1.18144094e+01,  1.86662362e+00,  9.35915820e-01, 9.99607785e-01,  8.82043169e-05]),
                np.array([2.43643187e-01, 9.11137379e-01, 1.32761033e-01, 4.83673726e+00, 7.63223497e+00, 7.85430786e+00, 1.00616406e+00, 1.13885975e+00, 1.73108170e+00, 2.12601846e-04])]
    g3 = lambda x,a,x0,a1,a2,a3,x1,x2,x3,s1,s2,s3,c : a*(a1*np.exp(-(x-x1-x0)**2/2/s1**2)+ \
                                                         a2*np.exp(-(x-x2-x0)**2/2/s2**2)+ \
                                                         a3*np.exp(-(x-x3-x0)**2/2/s3**2)+ \
                                                         +c)
    fcal1 = lambda x,a,x0 : g3(x,a,x00+x0+6.764,*cal_fres[0])
    fcal2 = lambda x,a,x0 : g3(x,a,x00+x0+7.183,*cal_fres[1])
    cal_functs = [fcal1,fcal2]
    return cal_functs

def get_calsum_funct(cal_functs):
    # output = calsum_funct = the continuous function with the four continuous calibration files. inputs to the function are four amplitudes and a single shift. All four calibration functions must shift together, but they have independent amplitude.
    if len(cal_functs) == 1:
        calsum_funct = lambda x,a0,x0 : cal_functs[0](x,a0,x0) # SimpleAuto
    elif len(cal_functs)==2:
        calsum_funct = lambda x,a0,a1,x0 : \
                            cal_functs[0](x,a0,x0) + \
                            cal_functs[1](x,a1,x0)
    elif len(cal_functs)==3:
        calsum_funct = lambda x,a0,a1,a2,x0 : \
                            cal_functs[0](x,a0,x0) + \
                            cal_functs[1](x,a1,x0) + \
                            cal_functs[2](x,a2,x0)
    elif len(cal_functs)==4:
        calsum_funct = lambda x,a0,a1,a2,a3,x0 : \
                            cal_functs[0](x,a0,x0) + \
                            cal_functs[1](x,a1,x0) + \
                            cal_functs[2](x,a2,x0) + \
                            cal_functs[3](x,a3,x0)
    else:
        calsum_funct = -1
    return calsum_funct

def get_filespecific_shift(tave_data,calsum_funct):
    # finds the file-specific shift of the continous calsum_funct from the time-averaged data. This assumes that alignment may change between files but the relative alignment of the calibration functions never changes and the alignment is constant during a single acquisition
    if calsum_funct.__code__.co_argcount == 4:
        p0 = (100,100,0)
        bounds = ((0,0,-3),(10**6,10**6,3))
    elif calsum_funct.__code__.co_argcount == 5:
            p0 = (100,100,100,0)
            bounds = ((0,0,0,-3),(10**6,10**6,10**6,3))
    elif calsum_funct.__code__.co_argcount == 6:
        p0 = (100,100,100,100,0)
        bounds = ((0,0,0,0,-3),(10**6,10**6,10**6,10**6,3))
    else:
        p0 = -1
        bounds = -1

    x = np.arange(len(tave_data))
    fres = curve_fit(calsum_funct,x,tave_data,
                     p0 = p0,
                     bounds = bounds)
    shift = fres[0][(calsum_funct.__code__.co_argcount)-2]

    return shift

def get_singleshift_calibs(cal_functs,shift=0,numpoints=32):
    # takes the continuous functions in the list of cal_functs and makes the descrete calib arrays to be used for linear fitting
    calibs = []
    x = np.arange(numpoints)
    for i in range(len(cal_functs)):
        cal = cal_functs[i](x,1,shift)
        cal = cal/np.sum(cal)
        calibs.append(cal)
    calibs = np.array(calibs)
    return calibs

def get_filespecific_calibs(tave_data, 
                            which_calibration_set = '4beads', 
                            shifttechnique=1,
                            x00 = 0):
    # combines together the above codes. Tell this function (1) the time-ave data and (2) the calib files to use and it gives you out the shifted goodness ready for linear fitting
    calfound = 1
    if which_calibration_set == 'SimpleAuto':
        cal_functs = get_cal_functs_SimpleAuto(x00=0) # four color beads
    elif which_calibration_set == '4beads':
        cal_functs = get_cal_functs_4beads(x00=x00) # four color beads
    elif which_calibration_set == 'CFPYFP':
        cal_functs = get_cal_functs_CFPYFP(x00=x00) # four color beads
    elif which_calibration_set == 'CFPmCherry':
        cal_functs = get_cal_functs_CFPmCherry(x00=x00) # four color beads
    elif which_calibration_set == 'YFPmCherry':
        cal_functs = get_cal_functs_YFPmCherry(x00=x00) # four color beads
    elif which_calibration_set == 'CFPYFPmCherry':
        cal_functs = get_cal_functs_CFPYFPmCherry(x00=x00) # CFP, YFP, and mCherry
    elif which_calibration_set == 'CFPYFPmCherryCy5':
        cal_functs = get_cal_functs_CFPYFPmCherryCy5(x00=x00) # CFP, YFP, mCherry, and Cy5
    elif which_calibration_set == 'CFPCy5':
        cal_functs = get_cal_functs_CFPCy5(x00=x00) # CFP and Cy5
    else:
        cal_functs = -1
        calibs = -1
        calfound = 0
        print('NO CALIBRATION FUNCTS FOUND - line 1049')

    if calfound:            
        calsum_funct = get_calsum_funct(cal_functs)
        
        if shifttechnique == 1 and which_calibration_set != 'SimpleAuto':
            shift = get_filespecific_shift(tave_data,calsum_funct)
        else:
            shift = 0
        calibs = get_singleshift_calibs(cal_functs,shift,len(tave_data))

    return calibs,shift

# %%
def smooth_IvT(inten,window=1,time_per_step = 0.00116):
    # window in seconds
    pnts = int(window/time_per_step)
    in_new = savgol_filter(inten, window_length=pnts, polyorder=2)
    return in_new


# %%
# New code for FCCS in iBio
# do full analysis from name of TIF/SIF file


def get_corrs_from_TIF(tiffile,
                       time_per_frame, maxt, numdt,
                       which_calibration_set='4beads',
                       shifttechnique=1,
                       changeWidth = 32,
                       SubtractBackground = False,
                       isolate_pos = False,
                       x00 = 0.5,
                       SmoothIvT = False):
    if tiffile[-4:] == '.tif':
        data,n_frames,wid,hei,frame_time_temp = load_TIF_funct(tiffile,isolate_pos = isolate_pos)
    elif tiffile[-4:] == '.sif':
        data,n_frames,wid,hei,frame_time_temp = load_SIF_funct(tiffile)
    
    if frame_time_temp > 0:
        time_per_frame = frame_time_temp
        
    if changeWidth > 0 and wid > changeWidth:
        start = int(wid/2 - changeWidth/2)
        end = int(start+changeWidth)
        wid = changeWidth
        data = data[:,start:end,:]
        
    data_lines = get_lines_from_array(data)
    tdata = np.mean(data_lines,axis=1) # time average line
    data_lines = data_lines - np.mean(tdata[-5:])
    tdata = tdata - np.mean(tdata[-5:])

    # data_lines2 = np.mean(data[:3,:,:],axis=0)
    # print('cn',data_lines)

    # print(which_calibration_set)
    calibs_now,shift = get_filespecific_calibs(tdata,
                                               which_calibration_set = which_calibration_set,
                                               shifttechnique = shifttechnique,
                                               x00 = x00)
    
    # print('cn',data_lines)simp_inten_vs_time

    if which_calibration_set == 'SimpleAuto':
        SubtractBackground = False
        
    if len(calibs_now)<4 and SubtractBackground == True:
        cn = np.zeros((3,32))
        cn[:2,:] = calibs_now
        cn[2,:] = (1/32)
        calibs_now = cn.copy()
        
    Mi = get_inv_matN_noBack(calibs_now)
    i0,i1,i2,i3 = fit_many_timesN_noBack(calibs_now,Mi,data_lines)
    # meanIs = np.array([np.mean(i0),np.mean(i1),np.mean(i2),np.mean(i3)])
    
    if SmoothIvT:
        i0s = smooth_IvT(i0,window=maxt,time_per_step = time_per_frame)
        i0 = i0-i0s+ np.mean(i0s)
        if len(i1)>1000:    
            i1s = smooth_IvT(i1,window=maxt,time_per_step = time_per_frame)
            i1 = i1-i1s+ np.mean(i1s)
        if len(i2)>1000:
            i2s = smooth_IvT(i2,window=maxt,time_per_step = time_per_frame)
            i2 = i2-i2s+ np.mean(i2s)
        if len(i3)>1000:
            i3s = smooth_IvT(i3,window=maxt,time_per_step = time_per_frame)
            i3 = i3-i3s + np.mean(i3s)

   # if len(calibs_now)<4 and SubtractBackground == True:
   #     calibs_now = calibs_now[:-1,:]
        
    
    if len(i0)>99:
        if len(calibs_now) > 0:
            dt,G0 = calc_xcorr(i0,i0,time_per_frame,maxt,numdt)
        if len(calibs_now) > 1:
            dt,G1 = calc_xcorr(i1,i1,time_per_frame,maxt,numdt)
            dt,G01 = calc_xcorr(i0,i1,time_per_frame,maxt,numdt)
        if len(calibs_now) > 2:
            dt,G2 = calc_xcorr(i2,i2,time_per_frame,maxt,numdt)
            dt,G02 = calc_xcorr(i0,i2,time_per_frame,maxt,numdt)
            dt,G12 = calc_xcorr(i1,i2,time_per_frame,maxt,numdt)
        if len(calibs_now) > 3:
            dt,G3 = calc_xcorr(i3,i3,time_per_frame,maxt,numdt)
            dt,G03 = calc_xcorr(i0,i3,time_per_frame,maxt,numdt)
            dt,G13 = calc_xcorr(i1,i3,time_per_frame,maxt,numdt)
            dt,G23 = calc_xcorr(i2,i3,time_per_frame,maxt,numdt)
        if len(calibs_now) > 4 or len(calibs_now) < 1:
            print('Wrong number of calibrations?? Line 1088')
            dt = -1
            corrs = -1
            # meanIs = -1
            
            print('Need updated which_calibration_set')
            
        if len(calibs_now) == 1:
            corrs = [G0]
            intens = [i0]
        elif len(calibs_now) == 2:
            corrs = [G0,G1,G01]
            intens = [i0,i1]
        elif len(calibs_now) == 3:
            corrs = [G0,G1,G2,G01,G02,G12]
            intens = [i0,i1,i2]
        elif len(calibs_now) == 4:
            corrs = [G0,G1,G2,G3,G01,G02,G03,G12,G13,G23]
            intens = [i0,i1,i2,i3]
    
        dt = dt*time_per_frame

    else:
        dt,corrs,intens,shift,tdata,calibs_now = [],[],[],[],[],[]
        print('   Movie too short to analyze. Line 1210')
        
    return dt,corrs,intens,shift,tdata,calibs_now,time_per_frame

def plot_all_calib_inten_G(dt,corrs,intens,tdata,calibs_now,title='',normalizeG = True,savename=''):
    f, (a0, a1, a2) = plt.subplots(3, 1, 
                                   dpi=100,
                                   figsize = (6,7),
                                   gridspec_kw={'height_ratios': [1, 1, 2]})
    
    if len(calibs_now) == 1:
        cols = 'g'
    elif len(calibs_now) == 2:
        cols = 'gr'
    elif len(calibs_now) == 3:
        cols = 'bgr'
    elif len(calibs_now) == 4:
        cols = 'bgrm'
    else:
        cols = get_colors(len(calibs_now))
        
    a0.plot(tdata,'ko',label='Time Ave Spectra')        
    calsum = calibs_now[0]*0
    for i,cal in enumerate(calibs_now):
        calsum = calsum + cal*np.mean(intens[i])
        a0.plot(cal*np.mean(intens[i]),'-',color=cols[i],label='Calib'+str(i+1))
    a0.plot(calsum,'-k',label='Sum fit')
    a0.set_xlabel('pixel number')
    a0.set_ylabel('intensity')
    a0.set_title(title)
    a0.legend(fontsize = 6, loc='upper right')
    
    x = np.arange(len(intens[0]))*dt[0]
    stdmax = 0
    normI = False
    for i,inten in enumerate(intens):
        if normI:
            a1.plot(x,inten-np.mean(inten),'-',color=cols[i],label='<Fluor'+str(i+1)+'>='+str(int(np.mean(inten))))
        else:
            a1.plot(x,inten,'-',color=cols[i],label='<Fluor'+str(i+1)+'>='+str(int(np.mean(inten))),alpha = 0.7)
        stdmax = np.max([stdmax,np.std(inten-np.mean(inten))])
    a1.set_xlabel('Time (sec)')
    if normI:
        a1.set_ylabel('Intensity-<Intensity>')
        a1.set_ylim([-stdmax*5,stdmax*5])
    else:
        a1.set_ylabel('Intensity')
    a1.set_xlim([0,x[-1]])
    a1.legend(loc='lower right')  
        
    nn = ''
    for i,G in enumerate(corrs[:len(intens)]):
        norm = 1
        if normalizeG:
            norm = np.mean(G[:2])
            nn = 'Normalized '
        a2.semilogx(dt,G/norm,'o',color=cols[i],label='Acorr'+str(i+1))
    sha = '<>^vsop<>v^'
    if len(corrs)<7:
        lab = ['12','13','23']
    else:
        lab = ['12','13','14','23','24','34']
    for i,G in enumerate(corrs[len(intens):]):
        norm = 1
        if normalizeG:
            norm = np.mean(G[:2])
        a2.semilogx(dt,G/norm,sha[i],color='k',label='Xcorr'+lab[i])
    a2.semilogx(dt,dt*0,'-',color='k')
    a2.set_xlabel('Tau (sec)')
    a2.set_ylabel(nn+'G(tau)')
    a2.legend(loc='upper right')  
    if normalizeG:
        a2.set_ylim([-0.3,1.3])
    a2.set_xlim([dt[0],dt[-1]])
    plt.tight_layout()
    if len(savename)>0:
        if type(savename) == type('a'):
            plt.savefig(savename)
        else:
            for sn in savename:
                plt.savefig(sn)
    plt.show()
    
    return f, (a0, a1, a2)


def plot_calibs_w_tave(tdata,calibs_now,which_calibration_set):
    # make pretty three-panel plot with
    # intensity vs pixel
    # intensity vs time
    # correlations
    
    Mi = get_inv_matN_noBack(calibs_now)
    ilist = fit_many_timesN_noBack(calibs_now,Mi,tdata)

    cols = get_colors(len(calibs_now),'rainbow')
    plt.figure(dpi=200)
    plt.plot(tdata,'ko',label='Time Ave Data')
    sum_calibs = calibs_now[0]*0
    for i in range(len(calibs_now)):
        plt.plot(ilist[i]*calibs_now[i],'-',color=cols[i],label='Calib '+str(i))
        sum_calibs += ilist[i]*calibs_now[i]
    plt.plot(sum_calibs,'k-',label='Fit')
    plt.legend()
    plt.xlabel('pixel')
    plt.ylabel('intensity')
    plt.title('Calibration set: '+which_calibration_set)
    plt.show()
    return ilist


# %% look through the folder and try to identify which files need to be analyzed
def new_tif(fold_new_data,analyzed_tifs):
    files = listdir(fold_new_data)
    newtifs = []
    for fi,f in enumerate(files):
        found = 0
        for of in analyzed_tifs:
            if f == of:
                found = 1
                break
        if found == 0:
            if f[-4:] == ".tif":
                newtifs.append(f)
    return newtifs    

def png_done(fold_new_data,analyzed_tifs):
    files = listdir(fold_new_data)
    for fi,f in enumerate(files):
        if f[-4:] == ".png":
            analyzed_tifs.append(f[:-4]+'.tif')
    return analyzed_tifs            

def get_subfolds(fold):
    files = os.listdir(fold)
    subfolds = []
    for fi,file in enumerate(files):
        if os.path.isdir(fold+file):
            if file[-1] != '\\':
                file = file+'\\'
            subfolds.append(file)
    return subfolds
                
def init_list(length):
    a = [[] for x in range(length)]
    return a

def get_subfold_tree(fold):
    # goes down 5 levels
    
    # subfoldtree[i] where i = folder level
    # subfoldtree[i][j] where j = folder number in first level
    # subfoldtree[i][j][k] where k = folder number in second level. Only valid when i >= 2, etc.
    # subfoldtree[i][j][k][L] where L = folder number in third level. Only valid when i >= 3, etc.

    levels = 5
    subfoldtree = init_list(levels)
    
    subfolds1 = get_subfolds(fold)
    subfoldtree[0] = subfolds1
    for i in range(1,levels):
        subfoldtree[i] = init_list(len(subfolds1))
      
    for i1,sf1 in enumerate(subfolds1):
        try:
            subfolds2 = get_subfolds(fold+sf1)
        except:
            subfolds2 = []
            
        subfoldtree[1][i1] = subfolds2
        for i in range(2,levels):
            subfoldtree[i][i1] = init_list(len(subfolds2))
      
            
        for i2,sf2 in enumerate(subfolds2):
            try:
                subfolds3 = get_subfolds(fold+sf1+sf2)
            except:
                subfolds3 = []
            subfoldtree[2][i1][i2] = subfolds3
            for i in range(3,levels):
                subfoldtree[i][i1][i2] = init_list(len(subfolds3))
    
        
            for i3,sf3 in enumerate(subfolds3):
                try:
                    subfolds4 = get_subfolds(fold+sf1+sf2+sf3)
                except:
                    subfolds4 = []
                subfoldtree[3][i1][i2][i3] = subfolds4
                for i in range(4,levels):
                    subfoldtree[i][i1][i2][i3] = init_list(len(subfolds4))
                    
                    for i4,sf4 in enumerate(subfolds4):
                        try:
                            subfolds5 = get_subfolds(fold+sf1+sf2+sf3+sf4)
                        except:
                            subfolds5 = []
                        subfoldtree[4][i1][i2][i3][i4] = subfolds5
                        
    return subfoldtree

def get_subfold_list(fold):
    # 5 levels
    subfoldtree = get_subfold_tree(fold)
    foldlist = []
    foldlist.append(fold)
    
    for i in range(len(subfoldtree[0])):
        foldlist.append(fold+subfoldtree[0][i])
        
        for j in range(len(subfoldtree[1][i])):
            foldlist.append(fold+subfoldtree[0][i]+subfoldtree[1][i][j])
                    
            for k in range(len(subfoldtree[2][i][j])):
                foldlist.append(fold+subfoldtree[0][i]+subfoldtree[1][i][j]+subfoldtree[2][i][j][k])
    
                for q in range(len(subfoldtree[3][i][j][k])):
                    foldlist.append(fold+subfoldtree[0][i]+subfoldtree[1][i][j]+subfoldtree[2][i][j][k]+subfoldtree[3][i][j][k][q])
    
                    for r in range(len(subfoldtree[4][i][j][k][q])):
                        foldlist.append(fold+subfoldtree[0][i]+subfoldtree[1][i][j]+subfoldtree[2][i][j][k]+subfoldtree[3][i][j][k][q]+subfoldtree[4][i][j][k][q][r])
                        
    return foldlist

                   
                        
