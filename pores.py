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
    
def pores(imagestack, nofimages, njump=1):
    """Computes pore positions and sizes from binary vesicle images

Images must be adjusted so that the vesicle approximately centered 
and the pore is in the right half of the image

Input: 
- njump : interval between images that are analysed (i.e. 100 means 
  pore size is computed for slices 1, 101, 201, 301...)

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

        #right half of the image
        rimg = image[:,centerx:]

        #label clusters on the image
        rimglabeled, noflabels = label(rimg, structure=s)
        
        #inner intersection of vertical midsection and top part of the vesicle
        innertop = np.max(np.nonzero(rimg[:centery,0]))
        #define which cluster it belongs to
        innertoplabel = rimglabeled[innertop, 0]
        
        #inner intersection of vertical midsection and bottom part of the vesicle
        innerbottom = np.min(np.nonzero(rimg[centery:,0])) + centery
        #define which cluster it belongs to
        innerbottomlabel = rimglabeled[innerbottom, 0]
        
        #if top and bottom is the same cluster - no pore, just putting edges in the center
        if innertoplabel == innerbottomlabel: 
            edgestop[i,:] = centery, centerx
            edgesbottom[i,:] = centery, centerx
        else:
            #take only one cluster representing top/bottom part of the vesicle
            ytop, xtop = np.nonzero(rimglabeled == innertoplabel)
            ybottom, xbottom = np.nonzero(rimglabeled == innerbottomlabel)
            #find angles between positive x-dir of horizontal midsection and all points
            angletop = np.arctan2(centery - ytop, xtop)
            anglebottom = np.arctan2(centery - ybottom, xbottom)
            
            #if bottom and top clusters are angularly overlapping - no pore
            if max(anglebottom) >= min(angletop):
                edgestop[i,:] = centery, centerx
                edgesbottom[i,:] = centery, centerx
            else:
                #take the point with minimal angle from the top cluster
                mintop = np.argmin(angletop)
                edgestop[i,:] = ytop[mintop], xtop[mintop]+centerx
                #take the point with maximal angle from the bototm cluster
                maxbottom = np.argmax(anglebottom)
                edgesbottom[i,:] = ybottom[maxbottom], xbottom[maxbottom]+centerx 
                
    poreradii = np.sqrt(np.sum((edgesbottom-edgestop)**2, axis=1)) / 2
    return framenos+1, edgestop, edgesbottom, poreradii   
    
def process_image(filename, nskip):
    """Load image, process it and save results"""
    tif, nofimg = load_imagestack(filename)
    framenos, edgestop, edgesbottom, rpores = pores(tif, nofimg, njump=nskip)
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