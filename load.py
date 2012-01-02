# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 00:31:09 2011

@author: Pavlo Shchelokovskyy
"""

from tifffile import TIFFfile

def load_images(filename):
    """Loads images from multipage TIFF"""
    tif = TIFFfile(filename)
    