#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

import wx

class CanvasPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

    def draw(self,x,y):
        self.axes.plot(x,y)


if __name__ == '__main__':
    app = wx.PySimpleApp()
    frame = wx.Frame(None, title='test',size=(800,450))
    panel = CanvasPanel(frame)
    t = np.arange(0,3,0.01)
    s = np.sin(2*np.pi*t)
    panel.draw(t,s)
    
    frame.Show()
    app.MainLoop()
