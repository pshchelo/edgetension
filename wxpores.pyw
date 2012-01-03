# -*- coding: utf-8 -*-
"""
Created on Sun Dec 25 15:15:39 2011

@author: family
"""

import wx

import matplotlib as mplt
mplt.use('WXAgg', warn=False)

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar2
from matplotlib.figure import Figure

import numpy as np
from scipy.stats import linregress

import widgets
import pores

class TensionsFrame(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self, parent, id, title = 'Pore Edge Tension')

        self.data = np.zeros((6,1))
        self.image = None
        
        self.panel = wx.Panel(self, -1)

        self.toolbar = widgets.SimpleToolbar(self, *self.ToolbarData())
        self.SetToolBar(self.toolbar)
        self.toolbar.Realize()

        self.statusbar = widgets.PlotStatusBar(self)
        self.SetStatusBar(self.statusbar)

        self.MakeParamsPanel()

        self.MakeImagePanel()

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.imgbox, 1, wx.GROW)

        hbox.Add(self.paramspanel, 0, wx.GROW)

        self.panel.SetSizer(hbox)

        self.SetFrameIcons(wx.ART_TIP, (16,24,32))
        hbox.Fit(self)

    def MakeImagePanel(self):
        self.figure = Figure(facecolor = widgets.rgba_wx2mplt(self.panel.GetBackgroundColour()))
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
        self.slider = widgets.DoubleSlider(self.panel, -1, (1, dim), 1, dim, gap=1)
        self.Bind(wx.EVT_SLIDER, self.OnSlide, self.slider)
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
        
        labelskip = wx.StaticText(self.paramspanel, -1, 'Skip')
        self.skipspin = wx.SpinCtrl(self.paramspanel, -1, '1', min=1, max=dim)
        self.Bind(wx.EVT_SPINCTRL, self.OnSkip, self.skipspin)
        flexsz.Add(labelskip, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.skipspin, 1, wx.GROW)

        labelzoom = wx.StaticText(self.paramspanel, -1, 'AutoZoom')
        self.zoomcb = wx.CheckBox(self.paramspanel, -1, )
        self.Bind(wx.EVT_CHECKBOX, self.Draw, self.zoomcb)
        flexsz.Add(labelzoom, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.zoomcb, 1, wx.GROW)

#TODO: change all below to agw.floatspin.FloatSpin
        labelradius = wx.StaticText(self.paramspanel, -1, 'Radius')
        self.radiusspin = wx.SpinCtrl(self.paramspanel, -1, '1', min=0, max=1000)
        self.Bind(wx.EVT_SPINCTRL, self.Draw, self.radiusspin)
        flexsz.Add(labelradius, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.radiusspin, 1, wx.GROW)

        labelfps = wx.StaticText(self.paramspanel, -1, 'FPS')
        self.fpsspin = wx.SpinCtrl(self.paramspanel, -1, '1', min=0, max=50000 )
        self.Bind(wx.EVT_SPINCTRL, self.Draw, self.fpsspin)
        flexsz.Add(labelfps, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.fpsspin, 1, wx.GROW)

        labelvisc = wx.StaticText(self.paramspanel, -1, 'Viscosity')
        self.viscspin = wx.SpinCtrl(self.paramspanel, -1, '1', min=0, max=1000)
        self.Bind(wx.EVT_SPINCTRL, self.Draw, self.viscspin)
        flexsz.Add(labelvisc, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALIGN_CENTER_VERTICAL)
        flexsz.Add(self.viscspin, 1, wx.GROW)

        paramsbox.Add(flexsz, 0)

        self.paramspanel.SetSizer(paramsbox)
    
    def OnSlide(self, evt):
        self.lowlabel.SetLabel('%i'%self.slider.GetLow())
        self.highlabel.SetLabel('%i'%self.slider.GetHigh())
        self.Draw(evt)
        evt.Skip()

    def OnSaveTxt(self, evt):
        savedlg = wx.FileDialog(self, 'Save data', self.GetParent().folder,
                            'tensions.dat', wildcard = '*.*',
                            style=wx.FD_SAVE|wx.FD_OVERWRITE_PROMPT)
        if savedlg.ShowModal() == wx.ID_CANCEL:
            evt.Skip()
            return
        datname = savedlg.GetPath()
        savedlg.Destroy()
        np.savetxt(datname, self.data)
        evt.Skip()

    def OnOpenTxt(self, evt):
        fileDlg = wx.FileDialog(self, message='Choose Pore Radius file...',
                                 wildcard='*.*', style=wx.FD_OPEN)
        if fileDlg.ShowModal() != wx.ID_OK:
            fileDlg.Destroy()
            return
        filename = fileDlg.GetPath()
        fileDlg.Destroy()
        
        try:
            self.data = np.loadtxt(filename, unpack=1)
            self.data[0] = np.log(self.data[0])
            maxframe = int(max(self.data[5]))
            minframe = int(min(self.data[5]))
            self.slider.SetMax(maxframe)
            self.slider.SetHigh(maxframe)
            self.slider.SetMin(minframe)
            self.slider.SetLow(minframe)
            self.OnSlide(evt)
        except Exception, e:
            self.OnError(str(e))
        evt.Skip()
        
    def OnOpenImg(self, evt):
        self.OnError('OpenImg not yet implemented')

    def OnDebug(self, evt):
        self.OnError('OnDebug not implemented')
    
    def OnSkip(self, evt):
        self.OnError('OnSkip not yet implemented')
        
    def OnError(self, msg):
        """
        Display an error dialog
        @param msg: error message to display (type = string)
        """
        errDlg = wx.MessageDialog(self, msg, "Error!", wx.ICON_ERROR)
        errDlg.ShowModal()
        errDlg.Destroy()

    def Draw(self, evt):
        low, high = self.slider.GetValue()
        self.lowvline.set_xdata(low)
        self.upvline.set_xdata(high)
        
        tofit = self.data[:,low:high]
        x = tofit[5]
        y = tofit[0]

        if self.zoomcb.GetValue():
            self.dataplot.set_data(x, y)
        else:
            self.dataplot.set_data(self.data[5], self.data[0])
        
        a, b, corrr, p, stderr = linregress(x, y)
        self.fitplot.set_data(x, a*x+b)
        
        v = self.viscspin.GetValue()
        r = self.radiusspin.GetValue()
        f = self.fpsspin.GetValue()
        
        modelgamma = lambda x: -1.5 * np.pi * v * r*r * f * x
        gamma =  modelgamma(a)
        gammastderr = np.abs(modelgamma(stderr))
        
        #since R is in micrometers and there is R**2, gamma is in picoNewtons
        title = '$\\gamma$ = %f $\\pm$ %f pN, from %i to %i'%(gamma, 
                                                        gammastderr, low, high)

        self.axes.set_title(title)
        self.axes.legend()

        self.axes.relim()
        self.axes.autoscale_view(tight=False)

        self.canvas.draw()
        evt.Skip()
        
if __name__ == '__main__':
    app = wx.PySimpleApp()
    frame = TensionsFrame(None, -1)
    frame.Show()
    app.MainLoop()