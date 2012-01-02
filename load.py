# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 00:31:09 2011

@author: Pavlo Shchelokovskyy
"""

from PIL import Image

#class TIFFImage(Image.Image):
#    """"""
#    def __init__(self, *args, **kwargs):
#        Image.Image.__init__(self, *args, **kwargs)
#        self._nofframes = None
#        
#    def _countframes(self):
#        count = 0
#        while True:
#            try:
#                self.seek(count)
#            except EOFError:
#                return count
#            else:
#                count += 1
#    
#    def __len__(self):
#        if self._nofframes:
#            return self._nofframes
#        else:
#            self._nofframes = self._countframes


def count_tiff_frames(tiffimage):
    count = 0
    while True:
        try:
            tiffimage.seek(count)
        except EOFError:
            return count
        else:
            count += 1



def load_images(filename):
    """Loads images from multipage TIFF"""
    tif = Image.open(filename)
    size = count_tiff_frames(tif)
    return tif, size
    