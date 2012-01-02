# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 00:26:39 2011

@author: Pavlo Shchelokovskyy

comments machine translated from French and not critically analyzed
are marked with FR:
"""

import numpy as np
import scipy as sp

from load import load

def pores(njump=100):
    """Computes pore sizes from binary vesicle images

argin: 
- njump : interval between images that are considered (100 by default, 
          meaning that pore size is computed for slices 1, 101, 201, 301...)

"""
    #load parameters
    names,nslices,Tacq,pixel_size,numero_manip,R,td,tf = load()
    #names - list of names of multipage TIFFs
    #nslices - corresponding number of images in each of multipage TIFFs
    #Tacq - timestep between frames (microsec)
    #pixel_size - physical scale of images (micron/pixel)
    #numero_manip - general number of vesicle studied 
    #               (can have several image sets recorded for the same vesicle)
    #R - radius of the vesicle (micron)
    #td, tf - lower and upper time limits of linear part of R^2*ln(r) (millisec)
    nmanips=len(nslices) #?overall number of image stacks studied
    
    #to start profiling timer
    #TODO: start profiling timer here
    
    #looping through image files (each is a multipage binary TIFF)
    for imanip in range(nmanips):
        #FR:shifting a space to the end of the file name
        pass #some code here
        
        #FR:init of the backup file
        pass #some array init here
        
        #looping through images in the multiimage TIFF
        for islice in range(0, nslices(imanip), njump):
            #FR:we initialize the variables each time through the loop
            pass #many various arrays are inited
            
            #reading the image
            #TODO: read one next image from multipage TIFF
            image = sp.misc.imread(names[imanip]) #temp hack just to define image
            #TODO: here it looks like expanding the image with 2-pixel black borders in all 4 sides
            
            #get image size and find the center
            nl, nc = image.shape()
            xc = floor(nl/2)
            yc = floor(nc/2)
            
            #FR:remove the lower part of the image
            im = image
            #? im[:,:yc-1] = 0
            
            #FR:one landmark or a point there
            #? [x,y]=find(im)
            #TODO: what does Matlab "find" do?
            
            #some debug plotting here
            pass
            
            #--------------------
            # start of analysis
            #--------------------
            
            #NE part
            #set the starting point
            x_NE=trouve_depart_NE(xc,yc,im,nl)
            #FR:init buddylist
            copains_NE=init_liste_copains(x_NE,yc,im)
            #FR: build buddylist
            copains_NE=construit_liste_copains(copains_NE,im)
            #FR:returns the index of friends with the largest theta
            z_NE,theta_NE=theta_max_NE(copains_NE,nl,nc,xc,yc)
            
            
            #NO part
            #set the starting point
            x_NO=trouve_depart_NO(xc,yc,im);
            #FR:init buddylist
            copains_NO=init_liste_copains(x_NO,yc,im);
            #FR: build buddylist
            copains_NO=construit_liste_copains(copains_NO,im);
            #FR:returns the index of friends with the largest theta
            z_NO,theta_NO=theta_min_NO(copains_NO,nl,nc,xc,yc);
            
            #TODO:find coordinates of arc ends
            #[xt_NE,yt_NE]=ind2sub([nl,nc],z_NE)
            #[xt_NO,yt_NO]=ind2sub([nl,nc],z_NO)

            #if there is no pore
            if theta_NE>theta_NO:
                xt_NE=x_NE
                yt_NE=yc
                xt_NO=x_NE
                yt_NO=yc
            
            #TODO:calculate half-distance between arc ends, multiply by pixel size
            #rpore=0.5*norm([(xt_NO-xt_NE) (yt_NO-yt_NE)]')*pixel_size(imanip)
            
            #some debug plotting here
            pass
            
            #--------------------
            # end of analysis
            #--------------------
            
            #we get the results
            #sauv=[sauv;[rpore xt_NO yt_NO xt_NE yt_NE islice]]
        
        #save results of this image and njump to a file
        sp.savetxt()
        
    #TODO:stop and present the timer
    
def trouve_depart_NE(xc,yc,im,nl):
    pass

def trouve_depart_NO(xc,yc,im):
    pass

def init_liste_copains(x_NE,yc,im):
    pass

def ajoute_voisins(x,y,im,cops):
    pass

def construit_liste_copains(copains_NE,im):
    pass

def theta_max_NE(copains_NE,nl,nc,xc,yc):
    pass

def theta_min_NO(copains_NO,nl,nc,xc,yc):
    pass