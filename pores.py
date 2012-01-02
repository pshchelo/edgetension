# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 00:26:39 2011

@author: Pavlo Shchelokovskyy
"""
from __future__ import division

import numpy as np
from scipy.misc import fromimage
 
def pore_edges(imagestack, nofimages, njump=1):
    """Computes pore sizes from binary vesicle images

Images must be rotated with the pore to the right and centered

Input: 
- njump : interval between images that are considered (100 by default, 
          meaning that pore size is computed for slices 1, 101, 201, 301...)

"""
    framenos = range(0, nofimages, njump)    
    edgestop = np.empty((len(framenos),2))
    edgesbottom = np.empty_like(edgestop)
    #looping through images in the multiimage TIFF
    for i, frameno in enumerate(framenos):
        imagestack.seek(frameno)
        image = fromimage(imagestack)
        
        #find the approximate center
        centerx = image.shape[1]//2
        centery = image.shape[0]//2
        
        #extract upper-right quadrant
        top = image[:centery,centerx:]
        ytop, xtop = np.nonzero(top)
        angletop = np.arctan2(centery - ytop, xtop)
        mintop = np.argmin(angletop)
        edgestop[i,:] = ytop[mintop], xtop[mintop]+centerx
        
        #extract lower-right quadrant
        bottom = image[centery:,centerx:]
        ybottom, xbottom = np.nonzero(bottom)
        anglebottom = np.arctan2(ybottom, xbottom)
        minbottom = np.argmin(anglebottom)
        edgesbottom[i,:] = ybottom[minbottom]+centery, xbottom[minbottom]+centerx
    
    return framenos, edgestop, edgesbottom
    
def pore_radii(edgestop, edgesbottom):
    """"""
    #distance between two edges of the pore
    dist = np.sqrt(np.sum((edgesbottom-edgestop)**2, axis=1))
    
    #if the edges are adjacent on y-coordinate - ther is no pore
    nopore = (edgesbottom[:,0]-edgestop[:,0] ==1)
    
    #calculate pore radius and set to 0 when there is no pore
    rpores = np.where(nopore, np.zeros_like(dist), dist) / 2
    return rpores