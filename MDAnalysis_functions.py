# -*- coding: utf-8 -*-
"""
Standard MDAnalysis functions

Created on Fri Aug 11, 2023
@author: CVKelly



"""
import numpy as np
import MDAnalysis as MDAnalysis
import os, sys, cv2
import nglview as nv

def plot_state(u,
                show_water = True,
                show_tails = True,
                show_tog = True,
                orient = [100,0,0,0,0,0,100,0,0,100,0,0,-70,-70,-42,1]):
                
    # two inputs here
    
    # initialize with ions
    a = u.select_atoms("resname ION")
    if a.n_atoms == 0:
        a = u.select_atoms("all")
    view = nv.show_mdanalysis(a)

    # display proteins
    a = u.select_atoms("not resname SAPI POPC DOPE TOG CL NA W")
    if a.n_atoms > 0:
        t = view.add_component(a)
        t.add_ribbon(color_scheme="atomindex")

    # display water
    if show_water:
        a = u.select_atoms("resname W WP PW")
        if a.n_atoms > 0:
            t = view.add_component(a)
            t.add_point(color = 'cyan') # ,sphereDetail=0.1)
    
    # display lipids
    cols = ["darkcyan","aquamarine","green","gray"]
    selelctions1 = ["POPC","SAPI","DOPE","TOG"] # resnames
    for si,sn in enumerate(selelctions1):
        if show_tails == False and si < 3:
            continue
        elif show_tog == False and si == 3:
            continue
        a = u.select_atoms("resname "+sn+" and not name PO4 NC3") 
        if a.n_atoms > 0:
            t = view.add_component(a)
            t.add_spacefill(color = cols[si])

    # display lipid head groups
    selelctions2 = ["PO4","NC3"] # atom names
    col2 = ["blue","red"]
    for si,sn in enumerate(selelctions2):
        a = u.select_atoms("name "+sn)
        if a.n_atoms > 0:
            t = view.add_component(a)
            t.add_spacefill(color = col2[si])
       
    # final display properties
    view.center()
    view.camera = 'orthographic'
    view._camera_orientation  = orient

    return view

def plot_traj(u,
                show_heads = True,
                show_water = True,
                show_tails = True,
                show_tog = True):

    # initialize with ions
    a = u.select_atoms("resname ION")
    if a.n_atoms == 0:
        a = u.select_atoms("all")
    view = nv.show_mdanalysis(a)

    # display proteins
    a = u.select_atoms("not resname SAPI POPC DOPE TOG CL NA W")
    if a.n_atoms > 0:
        t = view.add_trajectory(a)
        t.add_ribbon(color_scheme="atomindex")

    # display water
    if show_water:
        a = u.select_atoms("resname W WP PW")
        if a.n_atoms > 0:
            t = view.add_trajectory(a)
            t.add_point(color = 'cyan') # ,sphereDetail=0.1)
    
    # display lipids
    if show_tails:
        cols = ["darkcyan","aquamarine","green"]
        selelctions1 = ["POPC","SAPI","DOPE"] # resnames
        for si,sn in enumerate(selelctions1):
            a = u.select_atoms("resname "+sn+" and not name PO4 NC3") 
            if a.n_atoms > 0:
                t = view.add_trajectory(a)
                t.add_spacefill(color = cols[si])

    if show_tog:
        selelctions1 = ["TOG"] # resnames
        cols = ['gray']
        for si,sn in enumerate(selelctions1):
            a = u.select_atoms("resname "+sn+" and not name PO4 NC3") 
            if a.n_atoms > 0:
                t = view.add_trajectory(a)
                t.add_spacefill(color = cols[si])

    # display lipid head groups
    if show_heads:
        selelctions2 = ["PO4","NC3"] # atom names
        col2 = ["blue","red"]
        for si,sn in enumerate(selelctions2):
            a = u.select_atoms("name "+sn)
            if a.n_atoms > 0:
                t = view.add_trajectory(a)
                t.add_spacefill(color = col2[si])

    # final display properties

    view.center()
    view.camera = 'orthographic'
    view._camera_orientation  = orient

    return view


    # available colors: https://www.w3schools.com/colors/color_tryit.asp?hex=7FFFD4
    
    # # consider for nglview:
    ## Options for NGLView

    # # things working and cool
    # view.add_surface()

    # # maybe working? or just not applicable yet?
    # view.add_cartoon()

    # # things not working
    # view.setParameters({cameraType: "orthographic"})
    # view.camera.distance = 3
    # view.orientation = [[1000,0,0], [1000,0,0], [1000,0,0], [1000, 0, 0]]
    # view.rotation = [90,90,90]


def plot_traj_v2(u, 
                show_water = True, 
                show_tails = True,
                show_tag = True,
                minydisplay = 45,
                maxydisplay = 75,
                zoom = 125, # smaller is larger
                orient = -1 ):
    
    if orient == -1:
        # orient = [zoom, 0,0, 0, 0, 2*10**-14, zoom, 0, 0, zoom, 2*10**-14, 0, -70, -70, -42, 1]
        zoom = 125 # smaller is larger
        orient = [zoom, 0,0, 0, 0,0, zoom, 0, 0, zoom, 0, 0, -70, -70, -42, 1] # side view

    ylim =  " and prop y > "+str(minydisplay)+" and prop y < "+str(maxydisplay)
    
    # initialize with ions
    a = u.select_atoms("resname ION"+ylim, updating=True)
    if a.n_atoms == 0:
        a = u.select_atoms("all")
    view = nv.show_mdanalysis(a)

    # display proteins
    # a = u.select_atoms("not resname SAPI POPC DOPE TOG DOPC DOPE DOPS POPI CDL2 CL NA W")
    a = u.select_atoms("protein")
    if a.n_atoms > 0:
        t = view.add_trajectory(a)
        t.add_ribbon(color_scheme="atomindex")
       # print('Protein resnames:',list(np.unique(a.resnames)))

    # display water
    if show_water:
        a = u.select_atoms("resname W WP PW"+ylim, updating=True)
        if a.n_atoms > 0:
            t = view.add_trajectory(a)
            t.add_point(color = 'cyan') # ,sphereDetail=0.1)
          #  print('Water resnames:',list(np.unique(a.resnames)))

        
    # display lipids
    if show_tails:
        cols = ["darkcyan","darkcyan","aquamarine","green","gray","cyan","red","red"]
        selelctions1 = ["POPC","DOPC", "SAPI","DOPE","POPI","DOPS","CLD2","CDL2"] # resnames
        for si,sn in enumerate(selelctions1):
            a = u.select_atoms("resname "+sn+ylim, updating=True)
            if a.n_atoms > 0:
                t = view.add_trajectory(a)
                t.add_spacefill(color = cols[si])
            #    print('Tail resnames:',list(np.unique(a.resnames)))

    if show_tag:
        a = u.select_atoms("resname TOG"+ylim, updating=True)
        if a.n_atoms > 0:
            t = view.add_trajectory(a)
            t.add_spacefill(color = 'gray')       

    # display lipid head groups
    selelctions2 = ["PO4","PO41","PO42","NC3","NH3","C1","C2","C3","CN0"]
    col2 = ["blue","blue","blue","red","red","blue","blue","blue","blue"]
    for si,sn in enumerate(selelctions2):
        a = u.select_atoms("name "+sn+ylim, updating=True)
        if a.n_atoms > 0:
            t = view.add_trajectory(a)
            t.add_spacefill(color = col2[si], radius=0.1)
         #   print('Head resnames:',list(np.unique(a.resnames))

    # final display properties

    view.center()
    view.camera = 'orthographic'

    # orient = [-145,0, -10, 0, 0, 145,0,0,10,0,-145,0,-48,-2,-247] # top view
    view._camera_orientation  = orient

    return view
    
def test_funct():
    return 'test function is working'
    
def loc_to_str(x,y,z,atnum = 0,numspacesforatomnum=5):
    spaces = '       '
    x = str(x)
    xd = x.find('.')
    xs = x[:(xd+4)]

    x = str(y)
    xd = x.find('.')
    ys = x[:(xd+4)]

    x = str(z)
    xd = x.find('.')
    zs = x[:(xd+4)]

    xs = spaces+xs
    ys = spaces+ys
    zs = spaces+zs
    xs = xs[-8:]
    ys = ys[-8:]
    zs = zs[-8:]

    atnum += 1
    atst = str(atnum)
    atst = spaces+atst
    atst = atst[-numspacesforatomnum:]

    return xs,ys,zs,atst,atnum