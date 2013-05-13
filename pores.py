#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 15:15:39 2011

@author: Pavlo Shchelokovskyy

Papers referenced in comments:
[Ryham et al 2011] Ryham et al, Biophysical Journal 2011(v101)p2929
[Portet and Dimova 2010] Portet and Dimova, Biophysical Journal 2010(v99)p3264
"""
#TODO: progress indicator for image analysis duration
#FIXME: navbar-zoom is not initialized when opening new from image
#FIXME: AutoZoom is not working properly after having used navbar's pan and zoom
#FIXME: doubleslider indicators are not fitted when too big
from __future__ import division
import os

import numpy as np
from scipy.misc import fromimage
from scipy.ndimage import label
from scipy.stats import linregress
from scipy.optimize import curve_fit

from PIL import Image

import wx
from wx.lib.agw import floatspin as FS

import matplotlib as mplt
mplt.use('WXAgg', warn=False)

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar2
from matplotlib.figure import Figure

DATWILDCARD = "Data files (TXT, CSV, DAT)|*.txt;*.TXT;*.csv;*.CSV;*.dat;*.DAT | All files (*.*)|*.*"
IMGWILDCARD = "TIFF files (TIF, TIFF)|*.tif;*.TIF;*.tiff;*.TIFF | All files (*.*)|*.*"

TXT = 0
IMG = 1
MODES = {}
MODES[TXT] = 'TXT'
MODES[IMG] = 'IMG'

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
    """Loads images from multipage TIFF and counts them"""
    tif = Image.open(filename)
    nframes = count_tiff_frames(tif)
    return tif, nframes
    
def pores(imagestack, nofimages, njump=1):
    """Computes pore positions and sizes from binary vesicle images.

    Images must be adjusted so that the vesicle approximately centered 
    and the pore is in the right half of the image

    Input:
    imagestack -- PIL's multipage TIFF image
    nofimages -- number of frames in the imagestack
    njump -- analyse every njump'th frame (e.g. 100 means 
                  pore size is computed for frames 1, 101, 201, 301...),
                  defaults to 1 (i.e analyse all frames)

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
    return np.column_stack((poreradii, edgestop, edgesbottom, framenos+1)).T
    
def process_image(filename, nskip):
    """Load image, process it and save results
    
    Used in unattended batch processing mode, UI does not uses this
    
    """
    tif, nofimg = load_imagestack(filename)
    data = pores(tif, nofimg, njump=nskip)
    name, ext = os.path.splitext(filename)
    nameout = name+'_skip%i.txt'%nskip
    np.savetxt(nameout, data.T, fmt='%.7e')
    
def rgba_wx2mplt(wxcolour):
    """Convert wx.Colour instance of int tuple of rgba values 
    to range in 0-1 used by matplotlib.
    
    wxcolour -- wx.Colour instance
    
    """
    mpltrgba = []
    wxrgba = wxcolour.Get(includeAlpha=True)
    for item in wxrgba:
        converted = float(item)/255  # float division here!
        mpltrgba.append(converted)
    return tuple(mpltrgba)

class PlotStatusBar(wx.StatusBar):
    
    """Status Bar displaying XY coordinates of Matplotlib's axes."""
    
    def __init__(self, parent):
        wx.StatusBar.__init__(self, parent)
        self.SetFieldsCount(2)
        
    def SetPosition(self, evt):
        """Sets status bar text to current coordinates on the matplotlib plot/subplot
        
        evt -- matplotlib's motion_notify_event
        
        """
        if evt.inaxes:
            x = evt.xdata
            y = evt.ydata
            self.SetStatusText('x = %f, y = %f'%(x, y), 1)


class SimpleToolbar(wx.ToolBar):
    def __init__(self, parent, *buttons):
        """Construct and populate a simple wx.ToolBar 
        
        buttons -- tuple or list of ((Bitmap, shortName, longName, isToggle), Handler)
        
        """
        wx.ToolBar.__init__(self, parent)
        for button in buttons:
            buttonargs, handler = button
            tool = self.AddSimpleTool(-1, *buttonargs)
            self.Bind(wx.EVT_MENU, handler, tool)


class DoubleSlider(wx.Panel):
    
    '''Provides a panel with two sliders
    
    Useful for setting minimum and maximum, limits etc.
    Tryes to replicate the interface of wx.Slider as much as possible.
    Supports minimal gap between slider values.
    
    '''
    
    def __init__(self, parent, id, 
                 value = (1,100), min=1, max=100, gap = None, 
                 pos=wx.DefaultPosition, size=wx.DefaultSize,
                 panelstyle=wx.TAB_TRAVERSAL|wx.NO_BORDER, 
                 style=wx.SL_HORIZONTAL, name=wx.PanelNameStr):
        """Constructor of the DoubleSlider instance.
        
        Arguments:
        parent -- parent wxPython widget;
        id -- wx widget ID to assign to this instance;
        value -- initial values for two sliders, as tuple (dafaults to (1,100));
        min -- the minimal value for both sliders (defaults to 1);
        max -- maximal value for both sliders (defaults to 100);
        gap -- minimal difference between two sliders (defaults to None),
                  if set to some value sliders will move together to prevent
                  the diference between them becoming smaller that gap,
                  if still inside the sliders range, otherwise the sliders will freeze;
        pos -- wxPython position of the instance (defaults to wx.DefaultPosition);
        size -- wxPython size of the instance (defaults to wx.DefaultSize);
        panelstyle -- style to pass to the underlying wx.Panel widget
                            (defaults to wx.TAB_TRAVERSAL|wx.NO_BORDER);
        style -- style to pass to underlying wx.Slider widgets
                    (defaults to wx.SL_HORIZONTAL)
        name -- name 
        
        """
        wx.Panel.__init__(self, parent, id, pos=pos, size=size, style=panelstyle, name=name)
        low, high = value
        self.coupling = False
        if gap != None:
            self.gap = gap
            self.coupling = True
        self.lowslider = wx.Slider(self, -1, low, min, max, style=style)
        self.highslider = wx.Slider(self, -1, high, min, max, style=style)
        self.Bind(wx.EVT_SLIDER, self.OnSlide, self.lowslider)
        self.Bind(wx.EVT_SLIDER, self.OnSlide, self.highslider)
        if style & wx.SL_VERTICAL:
            sizer = wx.BoxSizer(wx.HORIZONTAL)
        else:  # horizontal sliders are default
            sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.lowslider, 0, wx.GROW)
        sizer.Add(self.highslider, 0, wx.GROW)
        self.SetSizer(sizer)
        self.Fit()
        
    def GetLow(self):
        """Value of the first(upper or left) slider"""
        return self.lowslider.GetValue()
    
    def SetLow(self, int):
        return self.lowslider.SetValue(int)
    
    def GetHigh(self):
        """Value of the second(lower or right) slider"""
        return self.highslider.GetValue()
        
    def SetHigh(self, int):
        return self.highslider.SetValue(int)
    
    def GetValue(self):
        """Value of both sliders as a tuple"""
        return self.GetLow(), self.GetHigh()
    
    def SetValue(self, value):
        low, high = value
        self.SetLow(low)
        self.SetHigh(high)
    
    Low = property(GetLow, SetLow)
    High = property(GetHigh, SetHigh)
    Value = property(GetValue, SetValue)
    
    def GetMin(self):
        """Minimal value for both sliders"""
        return self.lowslider.GetMin()
    
    def SetMin(self, int):
        self.lowslider.SetMin(int)
        self.highslider.SetMin(int)
    
    def GetMax(self):
        """Maximum value for both sliders"""
        return self.lowslider.GetMax()
        
    def SetMax(self, int):
        self.lowslider.SetMax(int)
        self.highslider.SetMax(int)
    
    def GetRange(self):
        """Range of both sliders as a tuple"""
        return self.lowslider.GetRange()
        
    def SetRange(self, min, max):
        self.lowslider.SetRange(min, max)
        self.highslider.SetRange(min, max)
    
    Min = property(GetMin, SetMin)
    Max = property(GetMax, SetMax)
    Range = property(GetRange, SetRange)
        
    def GetLineSize(self):
        """The amount of thumb movement when pressing arrow buttons"""
        return self.lowslider.GetLineSize()
        
    def SetLineSize(self, int):
        self.lowslider.SetLineSize(int)
        self.highslider.SetLineSize(int)
    
    def GetPageSize(self):
        """The amount of thumb movement when pressing Page Up/Down buttons"""
        return self.lowslider.GetPageSize()
        
    def SetPageSize(self, int):
        self.lowslider.SetPageSize(int)
        self.highslider.SetPageSize(int)
    
    LineSize = property(GetLineSize, SetLineSize)
    PageSize = property(GetPageSize, SetPageSize)
    
    def SlideLow(self):
        low, high = self.GetValue()
        min, max = self.GetRange()
        if low > max-self.gap:
            self.SetLow(max-self.gap)
            return
        if low > high-self.gap:
            self.SetHigh(low+self.gap)
    
    def SlideHigh(self):
        low, high = self.GetValue()
        min, max = self.GetRange()
        if high < min+self.gap:
            self.SetHigh(min+self.gap)
            return
        if high < low+self.gap:
            self.SetLow(high-self.gap)
    
    def GetCoupling(self):
        return self.coupling
    
    def SetCoupling(self, state):
        self.coupling = state
    
    def GetGap(self):
        if self.coupling:
            return self.gap
        else:
            return None
    
    def SetGap(self, gap):
        self.gap = gap

    def OnSlide(self, inevt):
        if self.coupling and inevt.GetEventObject() == self.lowslider:
            self.SlideLow()
        elif self.coupling and inevt.GetEventObject() == self.highslider:
            self.SlideHigh()
        event = wx.CommandEvent(inevt.GetEventType(), self.GetId()) 
        event.SetEventObject(self) 
        self.GetEventHandler().ProcessEvent(event)
        inevt.Skip()


class TensionsFrame(wx.Frame):
    """Main window.
    
    Displays plot, fitted data and results of the fit,
    plus all the relevant analysis parameters and fitting settings.
    """
    def __init__(self, parent, id, title, filename=None, skip=1):
        wx.Frame.__init__(self, parent, id, title=title)
        
        self.mode = None
        self.basetitle = title
        self.data = np.zeros((6,1))
        self.image = None
        self.nofimg = 0
        self.imagepath = filename
        self.skip = skip
        
        self.models = {'Slow closure (log-linear)':self.fitlinear, 
                       'Fast closure (quadratic)':self.fitquadratic}
        
        self.panel = wx.Panel(self, -1)

        self.toolbar = SimpleToolbar(self, *self.ToolbarData())
        self.SetToolBar(self.toolbar)
        self.toolbar.Realize()

        self.statusbar = PlotStatusBar(self)
        self.SetStatusBar(self.statusbar)

        self.MakeParamsPanel()

        self.MakeImagePanel()

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.imgbox, 1, wx.GROW)

        hbox.Add(self.paramspanel, 0, wx.GROW)

        self.panel.SetSizer(hbox)

        self.SetFrameIcons(wx.ART_TIP, (16,24,32))
        hbox.Fit(self)
        if self.imagepath:
            self.open_images(wx.PyCommandEvent(wx.EVT_BUTTON.typeId, self.GetId()))

    def MakeImagePanel(self):
        """Creates plot with navbar and doubleslider."""
        self.figure = Figure(facecolor = rgba_wx2mplt(self.panel.GetBackgroundColour()))
        self.canvas = FigureCanvas(self.panel, -1, self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.statusbar.SetPosition)
        self.axes = self.figure.add_subplot(111)
        self.axes.set_aspect('auto')

        self.dataplot, = self.axes.plot([], [], 'ro', label = 'Measured')
        self.fitplot, = self.axes.plot([],[], 'b-', label='Fit', lw=2)
        self.lowvline = self.axes.axvline(0, ls='--', c='blue')
        self.upvline = self.axes.axvline(1, ls='--', c='blue')
        self.fittext = self.axes.text(0.95, 0.95, '', 
                                      horizontalalignment='right', 
                                      verticalalignment='top', 
                                      transform=self.axes.transAxes, 
                                      bbox=dict(facecolor='grey', alpha=0.3)
                                      )
        self.axes.set_xlabel('time, s', size='large')

        self.navtoolbar = NavigationToolbar2(self.canvas)
        self.navtoolbar.Realize()
        
        dim = self.data.shape[1]
        if dim == 1:
            dim = 2
        self.slider = DoubleSlider(self.panel, -1, (1, dim), 1, dim, gap=1)
        self.Bind(wx.EVT_SLIDER, self.Draw, self.slider)
        self.lowlabel = wx.StaticText(self.panel, -1, '  %i'%self.slider.GetLow(), 
                                      style=wx.ALIGN_CENTER|wx.ST_NO_AUTORESIZE)
        self.highlabel = wx.StaticText(self.panel, -1, '  %i'%self.slider.GetHigh(), 
                                       style=wx.ALIGN_CENTER|wx.ST_NO_AUTORESIZE)

        self.imgbox = wx.BoxSizer(wx.VERTICAL)
        self.imgbox.Add(self.navtoolbar, 0, wx.GROW)
        self.imgbox.Add(self.canvas, 1, wx.GROW)
        
        labelbox = wx.BoxSizer(wx.VERTICAL)
        labelbox.Add(self.lowlabel, 1, wx.GROW)
        labelbox.Add(self.highlabel, 1, wx.GROW)

        sliderbox = wx.BoxSizer(wx.HORIZONTAL)
        sliderbox.Add(labelbox, 0, wx.GROW)
        sliderbox.Add(self.slider, 1, wx.GROW)

        self.imgbox.Add(sliderbox, 0, wx.GROW)

    def SetFrameIcons(self, artid, sizes):
        """Sets icon bundle for frame icons"""
        ib = wx.IconBundle()
        for size in sizes:
            ib.AddIcon(wx.ArtProvider.GetIcon(artid, size = (size,size)))
        self.SetIcons(ib)

    def ToolbarData(self):
        """Creates data to produce a toolbar"""
        bmpsavetxt = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE, wx.ART_TOOLBAR, (16,16))
        bmpopentxt = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_TOOLBAR, (16,16))
        bmpopenimg = wx.ArtProvider.GetBitmap(wx.ART_NEW, wx.ART_TOOLBAR, (16,16))
        bmpdebug = wx.ArtProvider.GetBitmap(wx.ART_MISSING_IMAGE, wx.ART_TOOLBAR, (16,16))
        return (
                ((bmpopenimg, 'Open Image File', 'Open image file', False),
                 self.OnOpenImg),
                ((bmpopentxt, 'Open Data file', 'Open pores data file', False),
                 self.OnOpenTxt),
                ((bmpsavetxt, 'Save Data File', 'Save pores data file', False),
                 self.OnSaveTxt),               
                ((bmpdebug, 'Debug', 'Show found pores', False),
                 self.OnDebug),               
                )

    def MakeParamsPanel(self):
        """Creates panel with analysis parameters and sets their default values"""
        dim = self.data.shape[1]
        if dim == 1:
            dim = 2
        
        self.paramspanel = wx.Panel(self.panel, -1)
        paramsstbox = wx.StaticBox(self.paramspanel, -1, 'Parameters')
        paramsbox = wx.StaticBoxSizer(paramsstbox, wx.VERTICAL)

        flexsz = wx.FlexGridSizer(cols=2, vgap=5, hgap=5)
        
        labelskip = wx.StaticText(self.paramspanel, -1, 'Take every (frame)')
        self.skipspin = wx.SpinCtrl(self.paramspanel, -1, '1', min=1, max=dim)
        self.Bind(wx.EVT_SPINCTRL, self.OnNewData, self.skipspin)
        flexsz.Add(labelskip, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.skipspin, 1, wx.GROW)

        labelzoom = wx.StaticText(self.paramspanel, -1, 'AutoZoom')
        self.zoomcb = wx.CheckBox(self.paramspanel, -1, )
        self.Bind(wx.EVT_CHECKBOX, self.Draw, self.zoomcb)
        flexsz.Add(labelzoom, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.zoomcb, 1, wx.GROW)

        labelscale = wx.StaticText(self.paramspanel, -1, 'Scale (um/pixel)')
        self.scalespin = FS.FloatSpin(self.paramspanel, -1, 
                                       value = 0.5, min_val=0, max_val=200, 
                                       increment=0.1, digits=4)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.scalespin)
        flexsz.Add(labelscale, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.scalespin, 1, wx.GROW)
        
        labeldiameter = wx.StaticText(self.paramspanel, -1, 'Diameter (pixel)')
        self.diameterspin = FS.FloatSpin(self.paramspanel, -1, 
                                       value = 100., min_val=0, max_val=200, 
                                       increment=1, digits=3)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.diameterspin)
        flexsz.Add(labeldiameter, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.diameterspin, 1, wx.GROW)

        labelfps = wx.StaticText(self.paramspanel, -1, 'Speed (fps)')
        self.fpsspin = FS.FloatSpin(self.paramspanel, -1, 
                                    value=5000., min_val=0, max_val=50000, 
                                    increment=1000, digits=3)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.fpsspin)
        flexsz.Add(labelfps, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.fpsspin, 1, wx.GROW)
        
        # default value is for 112 mM glucose in water
        labelvisc = wx.StaticText(self.paramspanel, -1, 'Viscosity (mPa*s)')
        self.viscspin = FS.FloatSpin(self.paramspanel, -1, 
                                     value=1.052, min_val=0, max_val=2000,
                                     increment=0.01, digits=3)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.viscspin)
        flexsz.Add(labelvisc, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.viscspin, 1, wx.GROW)
        
        # invaiant constant C appearing in [Ryham et al 2011]
        # default value is the one estimated in [Ryham et al 2011] for data of [Portet and Dimova 2010]
        labelfastconst = wx.StaticText(self.paramspanel, -1, 'Fast model constant')
        self.fastconstspin = FS.FloatSpin(self.paramspanel, -1, 
                                     value=8.44, min_val=0, max_val=10, 
                                     increment=0.01, digits=2)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.fastconstspin)
        flexsz.Add(labelfastconst, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.fastconstspin, 1, wx.GROW)
        
        paramsbox.Add(flexsz, 0)
        
        labelmodelchoice = wx.StaticText(self.paramspanel, -1, 'Fitting Model')
        self.modelchoice = wx.Choice(self.paramspanel, -1, 
                                     choices = self.models.keys())
        self.Bind(wx.EVT_CHOICE, self.Draw, self.modelchoice)
        self.modelchoice.SetSelection(1)
        
        paramsbox.Add(labelmodelchoice, 0)
        paramsbox.Add(self.modelchoice, 0)

        self.paramspanel.SetSizer(paramsbox)

    def OnSaveTxt(self, evt):
        """Handler for opening save file dialog"""
        savedlg = wx.FileDialog(self, 'Save pore data', '',
                            'pores.txt', wildcard = DATWILDCARD,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if savedlg.ShowModal() == wx.ID_CANCEL:
            evt.Skip()
            return
        datname = savedlg.GetPath()
        savedlg.Destroy()
        self.save_txt(datname)
        self.datapath = datname
        evt.Skip()
        
    def save_txt(self, path):
        """Saves data file to given path"""
        header = '#poreRadius_px\ty_top\tx_top\ty_bottom\tx_bottom\tframeNo\n'
        with open(path, 'w') as fout:
            fout.write(header)
            np.savetxt(fout, self.data.T, fmt='%g')

    def OnOpenTxt(self, evt):
        """Handler for opening open text file dialog"""
        fileDlg = wx.FileDialog(self, message='Choose Pore Radius file...',
                                 wildcard=DATWILDCARD, style=wx.FD_OPEN)
        if fileDlg.ShowModal() != wx.ID_OK:
            fileDlg.Destroy()
            return
        self.datapath = fileDlg.GetPath()
        fileDlg.Destroy()
        self.open_txt(evt)
        
    def OnOpenImg(self, evt):
        """Handler for opening open image dialog"""
        fileDlg = wx.FileDialog(self, message='Choose Pore Image file...',
                                 wildcard=IMGWILDCARD, style=wx.FD_OPEN)
        if fileDlg.ShowModal() != wx.ID_OK:
            fileDlg.Destroy()
            return
        self.imagepath = fileDlg.GetPath()
        fileDlg.Destroy()
        self.open_images(evt)
     
    def open_txt(self, evt):
        """Opens text file with data"""
        try:
            self.data = np.loadtxt(self.datapath, unpack=1)
        except Exception, e:
            self.OnError(str(e))
            evt.Skip()
            return
        self.mode = TXT
        self.imagepath = ''
        self.SetTitle('%s - %s'%(self.datapath, self.basetitle))
        self.OnNewData(evt)
    
    def open_images(self, evt):
        """Opens image file foe analysis"""
        try:
            self.image, self.nofimg = load_imagestack(self.imagepath)
        except Exception, e:
            self.OnError('Not an appropriate image: %s'%e)
            evt.Skip()
            return
        self.mode = IMG
        self.datapath = ''
        self.SetTitle('%s - %s'%(self.imagepath, self.basetitle))
        self.OnNewData(evt)
    
    def OnNewData(self, evt):
        """Initializes interface after opening image of text data file"""
        # Everything here is quite undefined if played on the maximum limit of skip's
        skip = self.skipspin.GetValue()
        if self.mode == IMG:
            self.data = pores(self.image, self.nofimg, skip)
            self.dataview = self.data
            l = self.nofimg
        else:
            cond = np.arange(self.data.shape[1])%skip == 0
            self.dataview = self.data[:,cond]
            l = self.data.shape[1]
            
        self.skipspin.SetRange(1, l)
        if skip > l:
            self.skipspin.SetValue(l)
        
        maxframe = int(max(self.data[5]))
        minframe = int(min(self.data[5]))
        self.slider.SetMax(maxframe)
        self.slider.SetHigh(maxframe)
        self.slider.SetMin(minframe)
        self.slider.SetLow(minframe)
        self.Draw(evt)
        
    def OnDebug(self, evt):
        """Handler of opening debug window event"""
        debugwindow = PoreDebugFrame(self, -1, self.data, images=self.image, 
                                     datpath=self.datapath, imgpath=self.imagepath)
        debugwindow.Show()
        
    def OnError(self, msg):
        """Handler for error messages.
        
        Displays an error dialog
        @param msg: error message to display (type = string)
        """
        errDlg = wx.MessageDialog(self, msg, "Error!", wx.OK|wx.ICON_ERROR)
        errDlg.ShowModal()
        errDlg.Destroy()

    def Draw(self, evt):
        """Make a fit and draw the plot"""
        
        params = self.getparams()
        
        FPS = params['fps']
        low, high = params['low'], params['high']
        
        self.lowlabel.SetLabel('%i'%low)
        self.highlabel.SetLabel('%i'%high)
#==============================================================================
#         Where FPS is present, it only affects plot display,
#         so that it is in seconds rather than in frame numbers;
#         all the calculations are carried out on frame numbers,
#         as it (most probably) gives less rounding etc artifacts
#==============================================================================
        
        model = self.models[params['model']]
        
        self.lowvline.set_xdata(low/FPS)
        self.upvline.set_xdata(high/FPS)
        
        r = self.dataview[0]
        f = self.dataview[5]
        
        ind = np.nonzero(np.logical_and(f >= low, f <= high))
        x = f[ind]
        y = r[ind]

        model(f, r, x, y, params)

        self.axes.relim()
        self.axes.autoscale_view(tight=False)

        self.canvas.draw()
        evt.Skip()

    def fitlinear(self, f, r, x, y, params):
        """Fit slow pore closure after [Portet and Dimova 2010]"""
        lnr = np.log(r)
        lny = np.log(y)
        self.axes.set_ylabel('$\ln (r_{pore})$', size='x-large')

        scale = params['scale']        
        Rv = params['Dv'] * scale / 2 # now Rv in micrometers
        visc = params['visc']
        FPS = params['fps']
        
        self.toggle_zoom(f/FPS, lnr, x/FPS, lny)
        
        a, b, corrr, p, sterr = linregress(x, lny)
        self.fitplot.set_data(x/FPS, a*x+b)
        
        modelgamma = lambda x: -1.5 * np.pi * visc * Rv*Rv * FPS * x
        gamma =  modelgamma(a)
        gammasterr = np.fabs(modelgamma(sterr))
        

        paramstxt = []
        paramstxt.append('%g FPS'%FPS)
        paramstxt.append('%g $\\mu$m/px'%scale)
        paramstxt.append('$R_v$=%g $\\mu$m'%Rv)
        paramstxt.append('$\\eta_s$=%g mPa*s'%(visc*1000))
        paramstxt.append('frames %i to %i'%(np.min(x), np.max(x)))
        
        self.axes.set_title(', '.join(paramstxt), size='medium')
        
        fitresults = []
        #since Rv is in micrometers and there is Rv**2, gamma is in picoNewtons
        fitresults.append('$\\gamma$=%g$\\pm$%g pN'%(gamma, gammasterr))
        fitresults.append('max $r_{pore}$ = %g $\\mu$m'%(max(r)*scale))
        self.fittext.set_text('\n'.join(fitresults))
    
    def fitquadratic(self, f, r, x, y, params):
        """Fit fast pore closure after [Ryham et al 2011]"""
        fastconst = params['fastconst']
        
        FPS = params['fps']
        visc = params['visc']
        scale = params['scale']
        Rv = params['Dv'] * scale / 2 # now Rv in micrometers
        
        self.axes.set_ylabel('$r_{pore}$', size='x-large')
        self.toggle_zoom(f/FPS, r, x/FPS, y)
        
        quad = lambda y, c, b, a:  c + b*y + a*y*y
        
        popt, pcov = curve_fit(quad, y, x)
        pstd = np.sqrt(np.diag(pcov))

        yfit = np.linspace(min(y), max(y), 100)
        self.fitplot.set_data(quad(yfit, *popt)/FPS, yfit)
        
        gamma = - 0.5 * fastconst * visc *scale*scale * FPS / popt[2] # in pN
        gammasterr = np.fabs(pstd[2]*gamma/popt[2])

        eta_m = 1000*0.25 * popt[1] *fastconst * visc * scale / popt[2] # in nPa*m*s
        eta_msterr = 1000*0.25 *fastconst * visc * scale * np.sqrt(pstd[1]**2 + 
                                (popt[1]*pstd[2]/popt[2])**2) / np.fabs(popt[2])
        
        t_c = popt[0] / FPS
        t_csterr = pstd[0] / FPS
        
        paramstxt = []
        paramstxt.append('%g FPS'%FPS)
        paramstxt.append('%g $\\mu$m/px'%scale)
        paramstxt.append('$C$=%g'%fastconst)
        paramstxt.append('$R_v$=%g $\\mu$m'%Rv)
        paramstxt.append('$\\eta_s$=%g mPa*s'%(visc*1000))
        paramstxt.append('frames %i to %i'%(np.min(x), np.max(x)))
        
        self.axes.set_title(', '.join(paramstxt), size='medium')
        
        fitresults = []
        fitresults.append('$\\gamma$=%f$\\pm$%f pN'%(gamma, gammasterr))
        fitresults.append('$\\eta_m$=%g$\\pm$%g nPa*s*m'%(eta_m, eta_msterr))
        fitresults.append('$t_c$=%g$\\pm$%g s'%(t_c, t_csterr))
        fitresults.append('max $r_{pore}$ = %g $\\mu$m'%(max(r)*scale))

        self.fittext.set_text('\n'.join(fitresults))
    
    def toggle_zoom(self, x, y, subx, suby):
        """Toggles data plot rescale parameters to imitate dynamic zoom"""
        if self.zoomcb.GetValue():
            self.dataplot.set_data(subx, suby)
        else:
            self.dataplot.set_data(x, y)
        self.navtoolbar.home()
        ## self.navtoolbar.update()
        self.navtoolbar.draw()

    def getparams(self):
        """Collects relevant for analisys parameters from UI"""
        params = {}
        params['visc'] = self.viscspin.GetValue() / 1000 #since input value is in mPa*s
        params['Dv'] = self.diameterspin.GetValue() # in pixels
        params['scale'] = self.scalespin.GetValue() # in micrometers per pixel
        params['fps'] = self.fpsspin.GetValue() # in 1/s
        params['fastconst'] = self.fastconstspin.GetValue()
        params['model'] = self.modelchoice.GetStringSelection()
        params['low'], params['high'] = self.slider.GetValue()
        
        return params

class PoreDebugFrame(wx.Frame):
    """Frame displaying image overlayed with found pore edges"""
    def __init__(self, parent, id, data, images=None, datpath='', imgpath=''):
        """"""
        wx.Frame.__init__(self, parent, id)
        
        self.data = data
        self.datpath = datpath
        if images and imgpath:
            self.images = images
            self.imgpath = imgpath
        else:
            fileDlg = wx.FileDialog(self, message='Choose Pore Image file...',
                                 wildcard=IMGWILDCARD, style=wx.FD_OPEN)
            if fileDlg.ShowModal() != wx.ID_OK:
                fileDlg.Destroy()
                self.Close()
            self.imgpath = fileDlg.GetPath()
            fileDlg.Destroy()
            self.images, imgno = load_imagestack(self.imgpath)
        
        title = 'Debug - %s // %s'%(self.imgpath, self.datpath) 
        self.SetTitle(title)
                
        self.panel = wx.Panel(self, -1)

        self.statusbar = PlotStatusBar(self)
        self.SetStatusBar(self.statusbar)
        
        self.frameslider = wx.Slider(self.panel, -1, value=0, 
                                minValue=0, maxValue=data.shape[1]-1)
        self.Bind(wx.EVT_SLIDER, self.OnSlide, self.frameslider)
        
        self.makePlot()
        navtoolbar = NavigationToolbar2(self.canvas)
        navtoolbar.Realize()
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.canvas, 1, wx.GROW)
        vbox.Add(navtoolbar, 0, wx.ALIGN_LEFT|wx.GROW)
        vbox.Add(self.frameslider, 0, wx.GROW)
        
        self.panel.SetSizer(vbox)        
        vbox.Fit(self)
        
    def makePlot(self):
        """creates plot with navbar etc"""
        self.figure = Figure(facecolor = rgba_wx2mplt(self.panel.GetBackgroundColour()))
        self.canvas = FigureCanvas(self.panel, -1, self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.statusbar.SetPosition)
        self.axes = self.figure.add_subplot(111)
        self.axes.set_aspect('equal')
        self.OnSlide(wx.EVT_SLIDER)
        
    def OnSlide(self, evt):
        """Handler of slider events.
        
        Redraws the plot depending on slider position.
        """
        self.axes.clear()
        index = self.frameslider.GetValue()
        r, y1, x1, y2, x2, frame = self.data[:,index]
        self.images.seek(frame-1)
        self.axes.plot((x1, x2),(y1,y2),'yo-', lw=3, ms=5, alpha=0.75)
        #TODO: make clear on what happens with flipping of the image
        self.axes.imshow(self.images.transpose(Image.FLIP_TOP_BOTTOM), 
                         aspect='equal', cmap='gray', interpolation="none")
        title = 'frame %i, Rpore = %g px'%(frame, r)
        self.axes.set_title(title)
        self.canvas.draw()
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Extract pore positions and radii.',
                    epilog="""with --nogui outputs TSV text file named imagename_skipXXX.txt, with 6 columns:
                                    pore radius in pixels,
                                    4 columns for x and y coordinates of pore edges,
                                    corresponding frame number.""")
    parser.add_argument('--file', '-f', default='', 
                        help="Name of the multipage b/w tiff file to process")
    parser.add_argument('--skip', '-s', default=1, type=int, 
                        help='Analyse only every SKIPth frame (default is every frame)')
    parser.add_argument('--test', '-t', type=int, default=0, 
                        help='Run the procedure TEST times and report the minimal of them')
    parser.add_argument('--nogui', '-n', default=False, action='store_true', 
                        help='Run in no GUI mode suitable for batch processing.')
    
    args = parser.parse_args()
    
    if args.test > 0:
        import timeit
        print min(timeit.repeat("process_image('%s', %i)"%(args.file, args.skip), 
                                "from __main__ import process_image", 
                                repeat=args.test, number=1))
    elif args.nogui:
        process_image(args.file, args.skip)
    else:    
        app = wx.App(False)
        frame = TensionsFrame(None, -1, 'Pore Edge Tension', 
                              filename=args.file, skip=args.skip)
        frame.Show()
        app.MainLoop()
