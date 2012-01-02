# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 00:31:09 2011

@author: Pavlo Shchelokovskyy
"""

from PIL import Image

def count_tiff_frames(tiffimage):
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
    