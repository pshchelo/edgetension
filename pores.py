# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 15:15:39 2011

@author: family
"""
from __future__ import division
import os

import numpy as np
from scipy.stats import linregress
from scipy.misc import fromimage
from scipy.ndimage import label

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
    return np.column_stack((poreradii, edgestop, edgesbottom, framenos+1)).T
    
def process_image(filename, nskip):
    """Load image, process it and save results"""
    tif, nofimg = load_imagestack(filename)
    data = pores(tif, nofimg, njump=nskip)
    name, ext = os.path.splitext(filename)
    nameout = name+'_skip%i.txt'%nskip
    np.savetxt(nameout, data.T, fmt='%.7e')
    
def rgba_wx2mplt(wxcolour):
    """
    Convert wx.Colour instance to float tuple of rgba values to range in 0-1 used by matplotlib.
    @param wxcolour: wx.Colour instance
    """
    mpltrgba = []
    wxrgba = wxcolour.Get(includeAlpha=True)
    for item in wxrgba:
        converted = float(item)/255
        mpltrgba.append(converted)
    return tuple(mpltrgba)

class PlotStatusBar(wx.StatusBar):
    '''Status Bar for wxPython VAMP frontend'''
    def __init__(self, parent):
        wx.StatusBar.__init__(self, parent)
        self.SetFieldsCount(2)
        
    def SetPosition(self, evt):
        """
        Set status bar text to current coordinates on the matplotlib plot/subplot
        @param evt: must be a matplotlib's motion_notify_event
        """
        if evt.inaxes:
            x = evt.xdata
            y = evt.ydata
            self.SetStatusText('x = %f, y = %f'%(x, y), 1)


class SimpleToolbar(wx.ToolBar):
    def __init__(self, parent, *buttons):
        """
        Construct and populate a simple wx.ToolBar 
        
        @param buttons: tuple or list of ((Bitmap, shortName, longName, isToggle), Handler)
        
        """
        wx.ToolBar.__init__(self, parent)
        for button in buttons:
            buttonargs, handler = button
            tool = self.AddSimpleTool(-1, *buttonargs)
            self.Bind(wx.EVT_MENU, handler, tool)


class DoubleSlider(wx.Panel):
    '''
    Provides a panel with two sliders to visually set 2 values (i.e. minimum and maximum, limits etc)
    '''
    def __init__(self, parent, id, 
                 value = (1,100), min=1, max=100, gap = None, 
                 pos=wx.DefaultPosition, size=wx.DefaultSize,
                 panelstyle=wx.TAB_TRAVERSAL|wx.NO_BORDER, 
                 style=wx.SL_HORIZONTAL, name=wx.PanelNameStr):
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
        """
        Value of the first(upper or left) slider
        """
        return self.lowslider.GetValue()
    def SetLow(self, int):
        return self.lowslider.SetValue(int)
    def GetHigh(self):
        """
        Value of the second(lower or right) slider
        """
        return self.highslider.GetValue()
    def SetHigh(self, int):
        return self.highslider.SetValue(int)
    def GetValue(self):
        """
        Value of both sliders as a tuple
        """
        return self.GetLow(), self.GetHigh()
    def SetValue(self, value):
        low, high = value
        self.SetLow(low)
        self.SetHigh(high)
    
    Low = property(GetLow, SetLow)
    High = property(GetHigh, SetHigh)
    Value = property(GetValue, SetValue)
    
    def GetMin(self):
        """
        Minimal value for both sliders
        """
        return self.lowslider.GetMin()
    def SetMin(self, int):
        self.lowslider.SetMin(int)
        self.highslider.SetMin(int)
    def GetMax(self):
        """
        Maximum value for both sliders
        """
        return self.lowslider.GetMax()
    def SetMax(self, int):
        self.lowslider.SetMax(int)
        self.highslider.SetMax(int)
    def GetRange(self):
        """
        Range of both sliders as a tuple
        """
        return self.lowslider.GetRange()
    def SetRange(self, min, max):
        self.lowslider.SetRange(min, max)
        self.highslider.SetRange(min, max)
    
    Min = property(GetMin, SetMin)
    Max = property(GetMax, SetMax)
    Range = property(GetRange, SetRange)
        
    def GetLineSize(self):
        """
        The amount of thumb movement when pressing arrow buttons
        """
        return self.lowslider.GetLineSize()
    def SetLineSize(self, int):
        self.lowslider.SetLineSize(int)
        self.highslider.SetLineSize(int)
    def GetPageSize(self):
        """
        The amount of thumb movement when pressing Page Up/Down buttons
        """
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
    def __init__(self, parent, id, title, filename=None, skip=1):
        wx.Frame.__init__(self, parent, id, title=title)

        self.basetitle = title
        self.data = np.zeros((6,1))
        self.image = None
        self.nofimg = 0
        self.imagepath = filename
        self.skip = skip
        
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
        self.figure = Figure(facecolor = rgba_wx2mplt(self.panel.GetBackgroundColour()))
        self.canvas = FigureCanvas(self.panel, -1, self.figure)
        self.canvas.mpl_connect('motion_notify_event', self.statusbar.SetPosition)
        self.axes = self.figure.add_subplot(111)
        self.axes.set_aspect('auto')

        self.dataplot, = self.axes.plot([], [], 'ro', label = 'Measured')
        self.fitplot, = self.axes.plot([],[], 'b-', label='Linear Fit', lw=2)
        self.lowvline = self.axes.axvline(0, ls='--', c='blue')
        self.upvline = self.axes.axvline(1, ls='--', c='blue')

        labelfont = {'fontsize':'large'}
        self.axes.set_xlabel('time, s', fontdict = labelfont)
        self.axes.set_ylabel('$\ln (r_{pore})$')

        navtoolbar = NavigationToolbar2(self.canvas)
        navtoolbar.Realize()
        
        dim = self.data.shape[1]
        if dim == 1:
            dim = 2
        self.slider = DoubleSlider(self.panel, -1, (1, dim), 1, dim, gap=1)
        self.Bind(wx.EVT_SLIDER, self.Draw, self.slider)
        self.lowlabel = wx.StaticText(self.panel, -1, '  %i'%self.slider.GetLow(), style=wx.ALIGN_CENTER|wx.ST_NO_AUTORESIZE)
        self.highlabel = wx.StaticText(self.panel, -1, '  %i'%self.slider.GetHigh(), style=wx.ALIGN_CENTER|wx.ST_NO_AUTORESIZE)

        self.imgbox = wx.BoxSizer(wx.VERTICAL)
        self.imgbox.Add(self.canvas, 1, wx.GROW)
        self.imgbox.Add(navtoolbar, 0, wx.GROW)

        labelbox = wx.BoxSizer(wx.VERTICAL)
        labelbox.Add(self.lowlabel, 1, wx.GROW)
        labelbox.Add(self.highlabel, 1, wx.GROW)

        sliderbox = wx.BoxSizer(wx.HORIZONTAL)
        sliderbox.Add(labelbox, 0, wx.GROW)
        sliderbox.Add(self.slider, 1, wx.GROW)

        self.imgbox.Add(sliderbox, 0, wx.GROW)

    def SetFrameIcons(self, artid, sizes):
        ib = wx.IconBundle()
        for size in sizes:
            ib.AddIcon(wx.ArtProvider.GetIcon(artid, size = (size,size)))
        self.SetIcons(ib)

    def ToolbarData(self):
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

        labelradius = wx.StaticText(self.paramspanel, -1, 'Radius (um)')
        self.radiusspin = FS.FloatSpin(self.paramspanel, -1, 
                                       value = 10., min_val=0, max_val=200, 
                                       increment=1, digits=3)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.radiusspin)
        flexsz.Add(labelradius, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.radiusspin, 1, wx.GROW)

        labelfps = wx.StaticText(self.paramspanel, -1, 'Speed (fps)')
        self.fpsspin = FS.FloatSpin(self.paramspanel, -1, 
                                    value=1000., min_val=0, max_val=50000, 
                                    increment=1000, digits=3)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.fpsspin)
        flexsz.Add(labelfps, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.fpsspin, 1, wx.GROW)

        labelvisc = wx.StaticText(self.paramspanel, -1, 'Viscosity (mPa*s)')
        self.viscspin = FS.FloatSpin(self.paramspanel, -1, 
                                     value=1.0, min_val=0, max_val=2000, 
                                     increment=0.01, digits=3)
        self.Bind(FS.EVT_FLOATSPIN, self.Draw, self.viscspin)
        flexsz.Add(labelvisc, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.viscspin, 1, wx.GROW)

        paramsbox.Add(flexsz, 0)

        self.paramspanel.SetSizer(paramsbox)

    def OnSaveTxt(self, evt):
        
        savedlg = wx.FileDialog(self, 'Save pore data', '',
                            'pores.txt', wildcard = DATWILDCARD,
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if savedlg.ShowModal() == wx.ID_CANCEL:
            evt.Skip()
            return
        datname = savedlg.GetPath()
        savedlg.Destroy()
        header = '#pore radius (pixel)\ty_top\tx_top\ty_bottom\tx_bottom\tframe No\n'
        with open(datname, 'w') as fout:
            fout.write(header)
            np.savetxt(fout, self.data.T, fmt='%g')
        evt.Skip()

    def OnOpenTxt(self, evt):
        fileDlg = wx.FileDialog(self, message='Choose Pore Radius file...',
                                 wildcard=DATWILDCARD, style=wx.FD_OPEN)
        if fileDlg.ShowModal() != wx.ID_OK:
            fileDlg.Destroy()
            return
        self.datapath = fileDlg.GetPath()
        fileDlg.Destroy()
        try:
            self.data = np.loadtxt(self.datapath, unpack=1)
        except Exception, e:
            self.OnError(str(e))
            evt.Skip()
            return
        self.image = None
        self.imagepath = ''
        self.nofimg = 0
        self.init_new_data()
        self.SetTitle('%s - %s'%(self.datapath, self.basetitle))
        self.Draw(evt)
        evt.Skip()
    
    def init_new_data(self):
        maxframe = int(max(self.data[5]))
        minframe = int(min(self.data[5]))
        self.slider.SetMax(maxframe)
        self.slider.SetHigh(maxframe)
        self.slider.SetMin(minframe)
        self.slider.SetLow(minframe)
        
    def OnOpenImg(self, evt):
        fileDlg = wx.FileDialog(self, message='Choose Pore Image file...',
                                 wildcard=IMGWILDCARD, style=wx.FD_OPEN)
        if fileDlg.ShowModal() != wx.ID_OK:
            fileDlg.Destroy()
            return
        self.imagepath = fileDlg.GetPath()
        fileDlg.Destroy()
        self.open_images(evt)
        
    def open_images(self, evt):        
        try:
            self.image, self.nofimg = load_imagestack(self.imagepath)
        except Exception, e:
            self.OnError('Not an appropriate image: %s'%e)
            return
        self.data = np.zeros((6,1))
        self.datapath = ''
        self.SetTitle('%s - %s'%(self.imagepath, self.basetitle))
        self.skipspin.SetRange(1, self.nofimg)
        self.skipspin.SetValue(self.skip)
        self.OnNewData(evt)
        
    def OnDebug(self, evt):
        debugwindow = PoreDebugFrame(self, -1, self.data, images=self.image, 
                                     datpath=self.datapath, imgpath=self.imagepath)
        debugwindow.Show()
    
    def OnNewData(self, evt):
        if self.image:
            self.data = pores(self.image, self.nofimg, self.skipspin.GetValue())
            self.init_new_data()
            self.Draw(evt)
        
    def OnError(self, msg):
        """
        Display an error dialog
        @param msg: error message to display (type = string)
        """
        errDlg = wx.MessageDialog(self, msg, "Error!", wx.ICON_ERROR)
        errDlg.ShowModal()
        errDlg.Destroy()

    def Draw(self, evt):
        visc = self.viscspin.GetValue() / 1000 #since input value is in mPa*s
        Rv = self.radiusspin.GetValue() #in 1/s
        FPS = self.fpsspin.GetValue() # in micrometers
        
        low, high = self.slider.GetValue()
        self.lowlabel.SetLabel('%i'%low)
        self.highlabel.SetLabel('%i'%high)
        
#==============================================================================
#         HINT: where /FPS is present, it only affects plot display,
#         so that it is in seconds instead of frame numbers;
#         all the real calculations are carried out on frame numbers,
#         as it (most probably) gives less rounding etc artifacts
#==============================================================================
        
        self.lowvline.set_xdata(low/FPS)
        self.upvline.set_xdata(high/FPS)
        
        lnr = np.log(self.data[0])
        f = self.data[5]
        
        ind = np.nonzero(np.logical_and(f >= low, f <= high))
        x = f[ind]
        y = lnr[ind]

        if self.zoomcb.GetValue():
            self.dataplot.set_data(x/FPS, y)
        else:
            self.dataplot.set_data(f/FPS, lnr)
        
        a, b, corrr, p, stderr = linregress(x, y)
        self.fitplot.set_data(x/FPS, a*x+b)
        
        modelgamma = lambda x: -1.5 * np.pi * visc * Rv*Rv * FPS * x
        gamma =  modelgamma(a)
        gammastderr = np.fabs(modelgamma(stderr))
        
        #since r is in micrometers and there is r**2, gamma is in picoNewtons
        title1 = '$\\gamma$ = %f $\\pm$ %f pN, frames %i to %i'%(gamma, 
                                                        gammastderr,
                                                        np.min(x), np.max(x))
        title2 = ' %g FPS, $R_v$ = %g$\\mu$m, $\\nu$ = %g Pa*s'%(FPS, Rv, visc)

        self.axes.set_title('\n'.join((title1, title2)))
        self.axes.legend()

        self.axes.relim()
        self.axes.autoscale_view(tight=False)

        self.canvas.draw()
        evt.Skip()

class PoreDebugFrame(wx.Frame):
    """"""
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
        self.axes.clear()
        index = self.frameslider.GetValue()
        r, y1, x1, y2, x2, frame = self.data[:,index]
        self.images.seek(frame-1)
        self.axes.plot((x1, x2),(y1,y2),'yo-', lw=3, ms=5, alpha=0.75)
        #TODO: make clear on what happens with flipping of the image
        self.axes.imshow(self.images.transpose(Image.FLIP_TOP_BOTTOM), aspect='equal', cmap='gray')
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
        app = wx.PySimpleApp(False)
        frame = TensionsFrame(None, -1, 'Pore Edge Tension', 
                              filename=args.file, skip=args.skip)
        frame.Show()
        app.MainLoop()