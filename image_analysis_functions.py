# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 19:11:32 2021

@author: cvkelly

last updated:
2022 02 18

"""
import cv2
import matplotlib
import numpy as np
from skimage import io
import tifffile
from os import rename, getcwd
from scipy.fft import fftn, ifftn, fftshift
from time import gmtime, strftime
import h5py
from scipy.ndimage import gaussian_filter
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

def get_drive():
    drive = getcwd()
    return drive[0]

def get_image_num(num,length=5):
    z = '00000'
    num = str(int(num))
    z = z[:-len(num)]
    return z+num
    
def import_npz(npz_file,allow_pickle=False):
    # This doesn't work as part of a Python module
    Data = np.load(npz_file,allow_pickle=allow_pickle)
    imported = 0
    for varName in Data:
        imported += 1
        globals()[varName] = Data[varName]
    return imported
        
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

def create_cmap(col):
    N = 256
    vals = np.ones((N, 4))
    for j in range(3):
        vals[:,j] = np.linspace(0,col[j], N)
    cmap = ListedColormap(vals)
    return cmap

def save_tiff(data,file_write,fmt='TXYCZ'):
    succ = 'save success'
    if len(data.shape) == 2:
        io.call_plugin('imsave',
                    plugin='tifffile',
                    file=file_write,
                    data=data,
                    imagej=True)
    elif len(data.shape)==3:    # data should be in format (x,y,t)
        io.call_plugin('imsave',
                        plugin='tifffile',
                        file=file_write,
                        data=data.transpose(2,0,1),
                        imagej=True,
                        metadata={'axes': 'TYX', 'fps': 10.0})
    elif len(data.shape)==4:    # data should be in format (x,y,c,t)
        io.imsave(file_write, data,
                    imagej=True,
                    metadata={'axes': 'TCYX', 'fps': 10.0})
    elif len(data.shape)==5:    # data should be in format (t,x,y,c,z)
        io.imsave(file_write, data,
                    imagej=True,
                    metadata={'axes': 'TCYXZ', 'fps': 10.0})
    else: 
        succ = 'Save TIF failed'
    return succ

def load_tif_stack(file):
    data = io.imread(file)
    # data = im.transpose(1,2,0)
    return data

def load_tif_stack_as_float(file):
    # I don't get why this code isn't working...

    data = io.imread(file).astype(np.float64)
    # if len(data.shape) == 3:
    #    data = data.transpose(1,2,0)
    return data


def get_OME_deet(file,descr='ExposureTime'):
    # other good options include
    # tmeta.imagej_metadata
    # tmeta.micromanager_metadata
    tmeta = tifffile.TiffFile(file)
    s = tmeta.ome_metadata
    w1 = s.find(descr)
    s = s[w1:]
    w1 = s.find('"')
    w2 = s[(w1+1):].find('"')
    val = s[(w1+1):(w1+w2+1)]
    try:
        val = float(val)
    except:
        pass
    return val
    
def load_MMtif_time(file,numtime=1,numcol=1):
    tmeta = tifffile.TiffFile(file)

    s = tmeta.ome_metadata
    d1 = 'DeltaT='
    d2 = 'DeltaTUnit'
    w1 = s.find(d1)
    w2 = s.find(d2)
    ttt = []
    cnt = 0
    while (w1 > 0 and cnt < 10**6):
        cnt += 1
        ttt.append(float(s[(w1+8):(w2-2)])) # originally in milliseconds
        s = s[(w2+5):]
        w1 = s.find(d1)
        w2 = s.find(d2)
    if not (numtime == 1 and numcol == 1):
        ttt = np.array(ttt).reshape((numttt,numcol))
    ttt = np.array(ttt)/1000/60
    return ttt # in minutes

def get_MM_sec_per_frame(file):
    ttt = load_MMtif_time(file)*60
    dt = ttt[1:]-ttt[:-1]
    dt = dt[dt<np.percentile(dt,98)]
    dt = dt[dt>np.percentile(dt,2)]
    sec_per_frame = np.mean(dt)
    return sec_per_frame

def stretch_image(data,frac_min=0.01,frac_max=0.01):
    # make the color bar useful
    # don't let the extreme pixels dominate the colorbar
    # data = data.astype(np.float64)
    # dl = data.flatten()
    # dl.sort()
    # argmin = int(dl.shape[0]*frac_min)
    # argmax = int(dl.shape[0]*frac_max)
    # old_min = dl[argmin]
    # old_max = dl[-argmax]-old_min

    # data = data-old_min
    # data = data/old_max
    # data[data<0]=0
    # data[data>1]=1
    # data = data*old_max
    # data = data+old_min
    top = np.percentile(data,frac_min*100)
    bot = np.percentile(data,100-frac_max*100)
    data[data<bot] = bot
    data[data>top] = top

    return data #,old_max,old_min


def smooth_line(a,sigma=3):
    b = gaussian_filter(a, sigma)
    return b

def smooth_image(a,size=3):
    kernel = np.ones((size,size),np.float32)
    kernel = kernel/np.sum(kernel)
    b = cv2.filter2D(a,-1,kernel)
    return b

#%% assorted image quantifications
def radial_average(c,maxr,numr):
    # assume the center of the image is the origin
    mid = int(c.shape[0]/2)
    a = 0
    if mid < c.shape[0]/2:
        a = 1
    x = np.arange(-mid,mid+a,1)

    X,Y = np.meshgrid(x,x)
    R = np.sqrt(X**2+Y**2)

    stepr = int(maxr/numr)
    stepr = np.max([stepr,1])

    r_edge = np.arange(0,maxr+stepr,stepr)
    g = np.zeros(len(r_edge)-1)
    for i in range(len(r_edge)-1):
        keep = np.all([R>r_edge[i],R<=r_edge[i+1]],axis=0)
        g[i] = np.mean(c[keep])

    r = (r_edge[1:]+r_edge[:-1])/2
    return r,g

def radial_average_orig(c,orig,numr=np.inf):

    # specify the origin
    maxr = np.min([c.shape[0]-orig[0],
                    c.shape[1]-orig[1],
                    orig[0],
                    orig[1]])

    x = np.arange(0,c.shape[0],1)-orig[1]
    y = np.arange(0,c.shape[1],1)-orig[0]

    X,Y = np.meshgrid(y,x)
    R = np.sqrt(X**2+Y**2)

    stepr = int(maxr/numr)
    stepr = np.max([stepr,1])

    r_edge = np.arange(0,maxr+stepr,stepr)
    g = np.zeros(len(r_edge)-1)
    for i in range(len(r_edge)-1):
        keep = np.all([R>r_edge[i],R<=r_edge[i+1]],axis=0)
        g[i] = np.mean(c[keep])

    r = (r_edge[1:]+r_edge[:-1])/2
    return r,g

def image_com(data):
    # data should be a 2D array
    background = np.mean([np.mean(data[:,:2]),
                          np.mean(data[:,-2:]),
                          np.mean(data[-2:,:]),
                          np.mean(data[:2,:])])
    data = data-background
    data[data<0] = 0
    xlist = np.arange(0,data.shape[0],1)
    ylist = np.arange(0,data.shape[1],1)
    com_y = np.sum(xlist*np.sum(data,axis=1))/np.sum(data)
    com_x = np.sum(ylist*np.sum(data,axis=0))/np.sum(data)
    return [com_x,com_y]


# %% do spatial correlations -- not verified
def spatial_autocorr_pre(image):
    b = fftn(image)
    b = np.abs(b)**2
    b = ifftn(b)
    b = np.real(fftshift(b))
    return b

def spatial_crosscorr_pre(image1,image2):
    b1 = fftn(image1)
    b2 = fftn(image2)
    b = b1 * np.conj(b2)
    b = ifftn(b)
    b = np.real(fftshift(b))
    return b

def spatial_corr_FFT(image1, image2, mask=0, maxr=100, numr=100):
    # uses FFT - unverified normalizations
    if np.all(mask==0):
        mask = np.ones(image1.shape)
    mask = mask/np.sum(mask)
    a = spatial_autocorr_pre(mask)
    b = spatial_crosscorr_pre(image1,image2)


    a[a<=0] = np.min(a[a>0])
    b[b<=0] = np.min(b[b>0])

    r1,g1 = radial_average(b,maxr,numr)
    r2,g2 = radial_average(a,maxr,numr)
    r = r1
    g = g1/g2

    return g,r

def spatial_corr_pre_manual(image1,image2,maxr = 50):
    dxlist = np.arange(-maxr,maxr+1,1)
    corr = np.zeros((len(dxlist),len(dxlist)))

    image1 = image1.astype(np.float64)
    image2 = image2.astype(np.float64)
    image1 = image1-np.mean(image1)
    image2 = image2-np.mean(image2)

    for i,dx in enumerate(dxlist):
        if dx > 0:
            a1 = image1[dx:,:]
            a2 = image2[:-dx,:]
        elif dx < 0:
            a2 = image2[-dx:,:]
            a1 = image1[:dx,:]
        else:
            a1 = image1.copy()
            a2 = image2.copy()
        for j,dy in enumerate(dxlist):
            if dy > 0:
                a3 = a1[:,dy:]
                a4 = a2[:,:-dy]
            elif dy < 0:
                a4 = a2[:,-dy:]
                a3 = a1[:,:dy]
            else:
                a3 = a1.copy()
                a4 = a2.copy()

            s1 = a3.shape[0]
            s2 = a3.shape[1]
            corr[i,j] = np.sum(a3 * a4) / s1 / s2

    # norm1 = np.mean(image1)
    # norm2 = np.mean(image2)
    # corr = corr # / norm1 / norm2
    return corr

def spatial_corr_manual(image1, image2, mask=0, maxr=100, numr=100):
    if np.all(mask==0):
        a = 1
    else:
        mask = mask/np.sum(mask)
        a = spatial_corr_pre_manual(mask,mask,maxr)

    if numr > maxr:
        numr = maxr

    b = spatial_corr_pre_manual(image1,image2,maxr)
    c = b/a
    r,g = radial_average(c,maxr,numr)

    return r,g

# %% rename files
def rename_file_number(file):
    z = '00000'
    sym = '-_ '
    loc = np.zeros(len(sym))
    for i,s in enumerate(sym):
        loc[i] = file.find(sym[i])
    loc[loc==-1] = np.inf
    minarg = np.argmin(loc)
    loc = loc[minarg]
    num = file[:loc]
    zn = z[:-len(num)]
    fnew = zn + file
    rename(file,fnew)
    return fnew

def print_ttt_str():
    t =  strftime("%Y-%m-%d %H:%M:%S", gmtime())
    return t

def time_string():
    t =  strftime("%Y-%m-%d %H:%M:%S", gmtime())
    return t

# %% for IMS analysis
def load_IMStif_time(file,numtime=1,numcol=1):
    tmeta = tifffile.TiffFile(file)

    s = tmeta.ome_metadata
    d1 = 'DeltaT='
    d2 = 'DeltaTUnit'
    w1 = s.find(d1)
    w2 = s.find(d2)
    time = []
    cnt = 0
    while (w1 > 0 and cnt < 10**6):
        cnt += 1
        time.append(float(s[(w1+8):(w2-2)]))
        s = s[(w2+5):]
        w1 = s.find(d1)
        w2 = s.find(d2)
    if not (numtime == 1 and numcol == 1):
        time = np.array(time).reshape((numtime,numcol))
    time = time/1000/60
    return time

def get_h5_file_info(h5_dataset):
    # Return the resolution levels, time points, channes, z levels, rows, cols, etc from a ims file
    # Pass in an opened f5 file
    # Get a list of all of the resolution options
    resolution_levels = list(h5_dataset)
    resolution_levels.sort(key = lambda x: int(x.split(' ')[-1]))

    # Get a list of the available time points
    time_points = list(h5_dataset[resolution_levels[0]])
    time_points.sort(key = lambda x: int(x.split(' ')[-1]))
    n_time_points = len(time_points)

    # Get a list of the channels
    channels = list(h5_dataset[resolution_levels[0]][time_points[0]])
    channels.sort(key = lambda x: int(x.split(' ')[-1]))
    n_channels = len(channels)

    # Get the number of z levels
    n_z_levels = h5_dataset[resolution_levels[0]][time_points[0]][channels[0]][
                   'Data'].shape[0]
    z_levels = list(range(n_z_levels))

    # Get the plane dimensions
    n_rows, n_cols = h5_dataset[resolution_levels[0]][time_points[0]][channels[0]][
                   'Data'].shape[1:]

    return resolution_levels, time_points, n_time_points, channels, n_channels, n_z_levels, z_levels, n_rows, n_cols

def get_timepoints(h5file):
    # get time each frame started in nanoseconds
    npnts = len(h5file['DataSetTimes']['Time'])
    timestart = np.zeros(npnts)
    for i in range(npnts):
        timestart[i] = h5file['DataSetTimes']['Time'][i][1]
    return timestart

def extract_h5_data(h5file):
    base_data = h5file['DataSet']
    file_info = get_h5_file_info(base_data)

    numtime = len(file_info[1])
    numcol = len(file_info[3])
    numzplanes = len(file_info[-3])
    size1 = file_info[-2]
    size2 = file_info[-1]

    data = np.zeros((numtime,size1,size2,numcol,numzplanes)).astype(np.uint16)

    for i in range(numtime):
        for j in range(numcol):
            for k in range(numzplanes):
                try:
                    d = np.array(base_data[file_info[0][0]][file_info[1][i]][file_info[3][j]]['Data'][k,:,:]).astype(np.uint32)
                    data[i,:,:,j,k] = d
                except:
                    print('  time',i,j,k,' not loading')
                    break
    data = remove_zeros(data)
    return data

def get_ims_data(file):
    if file[-4:] == '.ims':
        h5file = h5py.File(file,'r')
        data = extract_h5_data(h5file)
        times = get_timepoints(h5file) # in nanoseconds
        timestart = times[0]/10**9/60 # in min
        times = (times-times[0])/10**9/60 # in min
        h5file.close()
    else:
        data = 'file not ims'
        times = 'file not ims'
    return data,times

# %%
def bin_matrix(data,bx,by):
    # return a matrix of shape (n,m)
    n = data.shape[0]//bx
    m = data.shape[1]//by
    bs = data.shape[0]//n,data.shape[1]//m  # blocksize averaged over
    return np.reshape(np.array([np.sum(data[k1*bs[0]:(k1+1)*bs[0],k2*bs[1]:(k2+1)*bs[1]]) for k1 in range(n) for k2 in range(m)]),(n,m))

def bin_middle_two(data,bx,by):
    # remove empty rows and columns
    # assumes data is 4D array
    s = data.shape
    for i in range(s[0]):
        for j in range(s[3]):
            dtemp = bin_matrix(data[i,:,:,j],bx,by)
            if i == 0 and j == 0:
                newdata = np.zeros((s[0],dtemp.shape[0],dtemp.shape[1],s[3]))
            newdata[i,:,:,j] = dtemp
    return newdata

def remove_zeros(data):
    # note the if data shape
    # for 4D: only consider middle two dimensions
    # remove empty rows and columns
    s = data.shape
    if len(s) == 4:
        for i in range(s[1]):
            if np.sum(data[:,-1,:,:]) == 0:
                data = np.delete(data,data.shape[1]-1,1)
            else:
                break
        for i in range(s[2]):
            if np.sum(data[:,:,-1,:]) == 0:
                data = np.delete(data,data.shape[2]-1,2)
            else:
                break
    elif len(s) == 3:
        for i in range(s[0]):
            if np.sum(data[-1,:,0]) == 0:
                data = np.delete(data,data.shape[0]-1,0)
            else:
                break
        for i in range(s[1]):
            if np.sum(data[:,-1,0]) == 0:
                data = np.delete(data,data.shape[1]-1,1)
            else:
                break
    elif len(s) == 5: # for 5D IMS files
        for i in range(s[1]): # delete from the end
            if np.sum(data[:,-1,:,:,:]) == 0:
                data = np.delete(data,data.shape[1]-1,axis = 1)
            else:
                break
        for i in range(s[2]):
            if np.sum(data[:,:,-1,:,:]) == 0:
                data = np.delete(data,data.shape[2]-1,axis=2)
            else:
                break
                
    elif len(s) == 2:
        for i in range(s[0]):
            if np.sum(data[-1,:]) == 0:
                data = np.delete(data,data.shape[0]-1,0)
            else:
                break
        for i in range(s[1]):
            if np.sum(data[:,-1]) == 0:
                data = np.delete(data,data.shape[1]-1,1)
            else:
                break
    return data

def remove_maxes(data):
    # assumes data is 4D
    # only consider middle two dimensions
    # remove empty rows and columns
    s = data.shape
    data = np.max(data)-data
    if len(s) == 4:
        for i in range(s[1]):
            if np.sum(data[:,-1,:,:]) == 0:
                data = np.delete(data,data.shape[1]-1,1)
            else:
                break
        for i in range(s[2]):
            if np.sum(data[:,:,-1,:]) == 0:
                data = np.delete(data,data.shape[2]-1,2)
            else:
                break
    elif len(s) == 3:
        # consider the first two dimensions
        for i in range(s[0]):
            if np.sum(data[-1,:,0]) == 0:
                data = np.delete(data,data.shape[0]-1,0)
            else:
                break
        for i in range(s[1]):
            if np.sum(data[:,-1,0]) == 0:
                data = np.delete(data,data.shape[1]-1,1)
            else:
                break
    elif len(s) == 2:
        for i in range(s[0]):
            if np.sum(data[-1,:]) == 0:
                data = np.delete(data,data.shape[0]-1,0)
            else:
                break
        for i in range(s[1]):
            if np.sum(data[:,-1]) == 0:
                data = np.delete(data,data.shape[1]-1,1)
            else:
                break
    data = np.max(data)-data
    return data
    

# %% aligning images and movies 
def match_im_sizes(im1,im2,matchbiggest=False):
    s1 = im1.shape
    s2 = im2.shape
    if len(s1) > 2 or len(s2) > 2:
        print('Images aren''t single frames')
        return 'Images aren''t single frames'
    
    if matchbiggest:
        newwid = np.max([s1[0],s2[0]])
        newhei = np.max([s1[1],s2[1]])
    else:
        newwid = np.min([s1[0],s2[0]])
        newhei = np.min([s1[1],s2[1]])
    
    im1 = cv2.resize(im1, None, fy = newwid/s1[0], fx = newhei/s1[1])
    im2 = cv2.resize(im2, None, fy = newwid/s2[0], fx = newhei/s2[1])
    return im1, im2

def align_movie(data,r=15):
# calculate how much shifting is necessary to align the frames
# keep the full image for each dimension's shift
# data shape = num frames, x-dim, y-dim
# r = range of maximum shift

    shift = np.zeros((data.shape[0],2))
        
    rlist = np.arange(-(r-1),r,1)
    msd = np.zeros((len(rlist),2))
    dgoal = np.mean(data.astype(np.float32),axis=0)
    dgoal = dgoal[r:-r,r:-r]
    dgoal = dgoal - np.mean(dgoal)
    dgoal = dgoal/np.std(dgoal) 

    for frame in range(data.shape[0]):
        for ri,rnow in enumerate(rlist):
            d1 = data[frame,(r+rnow):-(r-rnow),r:-r]
            d2 = data[frame,r:-r,(r+rnow):-(r-rnow)]

            d1 = d1 - np.mean(d1)
            d1 = d1/np.std(d1)
            d2 = d2 - np.mean(d2)
            d2 = d2/np.std(d2)
            
            msd[ri,0] = np.mean((d1-dgoal)**2)
            msd[ri,1] = np.mean((d2-dgoal)**2)
            shift[frame,0] = rlist[np.argmin(msd[:,0])]
            shift[frame,1] = rlist[np.argmin(msd[:,1])]

    shift = shift.astype(np.int16)
    ms = int(np.max(np.abs(shift)))+1 # max shift used
    data2 = np.zeros(data[:,ms:-ms,ms:-ms].shape).astype(np.uint16)
    for i in range(data.shape[0]):
        data2[i,:,:] = data[i,(ms+shift[i,0]):(-ms+shift[i,0]),(ms+shift[i,1]):(-ms+shift[i,1])]

    return data2, shift
    
def high_pass_filter(img,size=50):
    if not size%2:
        size +=1
    kernel = np.ones((size,size),np.float32)/(size*size)
    filtered= cv2.filter2D(img,-1,kernel)
    filtered = img.astype('float32') - filtered.astype('float32')
    filtered = filtered + 127*np.ones(img.shape, np.uint8)
    return filtered

def low_pass_filter(img,size=50):
    if not size%2:
        size +=1
    kernel = np.ones((size,size),np.float32)/(size*size)
    filtered= cv2.filter2D(img,-1,kernel)
    return filtered
    
    
# %% correlation analyises 
def calc_corr_mag(im1,im2):
    im1 = im1-np.mean(im1)
    im2 = im2-np.mean(im2)
    im3 = im1*im2
    c = np.mean(im3)
    return im3,c

def get_G(smooth_range):
    x = np.arange(-smooth_range*3,smooth_range*3,1)
    y = x.copy()
    X,Y = np.meshgrid(x,y)
    G = np.exp(-(X**2+Y**2)/2/smooth_range**2)
    return G

def calc_corr_mag_time(data1,data2,smooth_range):
    # data format = t,x,y
    num_time = data1.shape[0]
    corrmag = np.zeros(num_time)
    # G = get_G(smooth_range)
    for t in range(num_time):
        # im1 = convolve2d(data1[:,:,t],G)
        # im2 = convolve2d(data2[:,:,t],G)
        im1 = data1[t,:,:]
        im2 = data2[t,:,:]

        # im1 = iaf.bin_matrix(im1,16,16)
        # im2 = iaf.bin_matrix(im2,16,16)
        im3,c = calc_corr_mag(im1,im2)
        corrmag[t] = c
    return corrmag

def get_single_cell_masks(imshape,points=[[20,20,10]]):
    points = np.array(points)
    bigarray = np.zeros(imshape)
    x = np.arange(bigarray.shape[0])
    y = np.arange(bigarray.shape[1])
    X,Y = np.meshgrid(x,y)
    Rlist = [bigarray+10**6]
    for pti,pt in enumerate(points):
        X,Y = np.meshgrid(x-pt[0],y-pt[1])
        Rnow = np.sqrt(X**2+Y**2)
        Rnow[Rnow>pt[2]] = 10**6
        Rlist.append(Rnow)
    Rlist = np.array(Rlist)
    for i in range(bigarray.shape[0]):
        for j in range(bigarray.shape[1]):
            bigarray[i,j]=np.argmin(Rlist[:,i,j])
    return bigarray

def pearson_cc(im1,im2,background=0,remove_zeros=True):
    im1 = np.array(im1)
    im2 = np.array(im2)
    im1 = im1.ravel().astype(np.float32)
    im2 = im2.ravel().astype(np.float32)
    
    if remove_zeros:
        keep = np.all([im1>0,im2>0],axis=0)
        im1 = im1[keep]
        im2 = im2[keep]
    
    im1 = im1-background
    im2 = im2-background    

    n = len(im1)
    t1 = n*np.sum(im1*im2)-np.sum(im1)*np.sum(im2)
    t2 = n*np.sum(im1**2)-np.sum(im1)**2
    t3 = n*np.sum(im2**2)-np.sum(im2)**2
    Pcorr = t1/np.sqrt(t2*t3)
    return Pcorr
    
    