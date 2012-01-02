# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 00:26:39 2011

@author: Pavlo Shchelokovskyy
"""
from __future__ import division
import os
    
import numpy as np
from scipy.misc import fromimage
from scipy.ndimage import label

from PIL import Image

def count_tiff_frames(tiffimage):
    """Counts frames in a multipage TIFF file"""
    count = 0
    while True:
        try:
            tiffimage.seek(count)
        except EOFError:
            return count
        else:
            count += 1

def load_imagestack(filename):
    """Loads images from multipage TIFF"""
    tif = Image.open(filename)
    size = count_tiff_frames(tif)
    return tif, size
    
def pore_edges(imagestack, nofimages, njump=1):
    """Computes pore sizes from binary vesicle images

Images must be rotated with the pore to the right and centered

Input: 
- njump : interval between images that are considered (100 by default, 
          meaning that pore size is computed for slices 1, 101, 201, 301...)

"""
    framenos = np.arange(0, nofimages, njump)
    edgestop = np.empty((len(framenos),2))
    edgesbottom = np.empty_like(edgestop)

    #looping through images in the multiimage TIFF
    for i, frameno in enumerate(framenos):
        imagestack.seek(frameno)
        image = fromimage(imagestack)
        #find the approximate center
        centerx = image.shape[1]//2
        centery = image.shape[0]//2
        
        s = np.ones((3,3))#structure defining which directions are taken as neighbours

#==============================================================================
# Upper-right quadrant
#==============================================================================
        #extract upper-right quadrant
        top = image[:centery,centerx:]
        #find the innermost point on the left edge of the image
        innertop = np.max(np.nonzero(top[:,0]))
        #label interconnected features/clusters
        toplabeled, noftoplabels = label(top, structure=s)
        #find which cluster corresponds to innermost point on the left edge
        #this is the vesicle
        innertoplabel = toplabeled[innertop, 0]
        #take coordinates comprising this cluster
        ytop, xtop = np.nonzero(toplabeled == innertoplabel)
        #find point with minimal angle
        angletop = np.arctan2(centery - ytop, xtop)
        mintop = np.argmin(angletop)
        #shift coordinates back to frame of the whole image
        edgestop[i,:] = ytop[mintop], xtop[mintop]+centerx
#==============================================================================
# Lower-right quadrant - procedure similar to the upper-right one
#==============================================================================
        bottom = image[centery:,centerx:]
        innerbottom = np.min(np.nonzero(bottom[:,0]))
        bottomlabeled, nofbottomlabels = label(bottom, structure=s)
        innerbottomlabel = bottomlabeled[innerbottom, 0]
        ybottom, xbottom = np.nonzero(bottomlabeled == innerbottomlabel)
        anglebottom = np.arctan2(ybottom, xbottom)
        minbottom = np.argmin(anglebottom)
        edgesbottom[i,:] = ybottom[minbottom]+centery, xbottom[minbottom]+centerx 
    return framenos+1, edgestop, edgesbottom
    
def pore_radii(edgestop, edgesbottom):
    """Computes radii of pore given its edges"""
    
    #distance between two edges of the pore
    dist = np.sqrt(np.sum((edgesbottom-edgestop)**2, axis=1))
    
    #if the edges are adjacent on y-coordinate - there is no pore
    nopore = (edgesbottom[:,0]-edgestop[:,0] ==1)

    #calculate pore radius and set to 0 when there is no pore
    rpores = np.where(nopore, np.zeros_like(dist), dist) / 2

    return rpores
    
def process_image(filename, nskip):
    """Load image, process it and save results"""
    tif, nofimg = load_imagestack(filename)
    framenos, edgestop, edgesbottom = pore_edges(tif, nofimg, njump=nskip)
    rpores = pore_radii(edgestop, edgesbottom)
    
    name, ext = os.path.splitext(filename)
    nameout = name+'_skip%i.txt'%nskip
    out = np.column_stack((rpores, edgestop, edgesbottom, framenos))
    np.savetxt(nameout, out, fmt='%.7e')
    
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Extract pore positions and radii.')
    parser.add_argument('filename', help="Name of the multipage b/w tiff file")
    parser.add_argument('--skip', default=1, type=int, help='Analyse only every SKIPth frame')
    args = parser.parse_args()
    
    process_image(args.filename, args.skip)