#!/usr/bin/env python


"""
Given output from PRESTO's single_pulse_search.py routine, plot those events in an interactive session.
"""

__version__ = '0.1'
__revision__ = '1'
__author__ = 'Bradley Meyers\nSteven Tremblay'

import os
import wx
import pandas
import matplotlib
matplotlib.use('WXAgg')
matplotlib.interactive(True)

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg, NavigationToolbar2Wx


## ID numbers for wxPython referencing to objects

ID_OPEN = 10
ID_SAVE = 11
ID_SAVE_AS = 12
ID_EXIT = 13

ID_COLOR_AUTO = 20
ID_COLOR_ADJUST = 21
ID_COLOR_MAP_PAIRED = 2200
ID_COLOR_MAP_SPECTRAL = 2201
ID_COLOR_MAP_BONE = 2202
ID_COLOR_MAP_JET = 2203
ID_COLOR_MAP_EARTH = 2204
ID_COLOR_MAP_HEAT = 2205
ID_COLOR_MAP_NCAR = 2206
ID_COLOR_MAP_RAINBOW = 2207
ID_COLOR_MAP_STERN = 2208
ID_COLOR_MAP_GRAY = 2209
ID_COLOR_INVERT = 23
ID_COLOR_STRETCH_LINEAR = 2400
ID_COLOR_STRETCH_LOG = 2401
ID_COLOR_STRETCH_SQRT = 2402
ID_COLOR_STRETCH_SQRD = 2403
ID_COLOR_STRETCH_ASINH = 2404

ID_HELP = 70
ID_ABOUT = 71





class MainWindow(wx.Frame):
    def __init__(self, parent, id):
        wx.Frame.__init__(self,parent,title="Interactive Single Pulse Plotter",size=(1000,600)) 
        self.dirname = ''
        self.filename = ''
        self.data = None
        self.edited = False


    def render(self):
        self.init_UI()
        self.init_Events()
        self.Show()
        self.SetClientSize((1000,600))


    def init_UI(self):
        self.statusbar = self.CreateStatusBar()  ## Create status bar at bottom of Frame

        font = wx.SystemSettings_GetFont(wx.SYS_SYSTEM_FONT)  ## Set font to system default
        font.SetPointSize(10)                                 ## Set font size to 10

        ## Create MenuBar objects to populate
        menuBar = wx.MenuBar()
        fileMenu = wx.Menu()
        colorMenu = wx.Menu()
        helpMenu = wx.Menu()
        
        ## File menu
        open = wx.MenuItem(fileMenu, ID_OPEN, '&Open')
        save = wx.MenuItem(fileMenu, ID_SAVE, '&Save')
        saveAs = wx.MenuItem(fileMenu, ID_SAVE_AS, 'Save &As')
        exit = wx.MenuItem(fileMenu, ID_EXIT, 'E&xit')
        
        for item in [open, save, saveAs, exit]:
            fileMenu.AppendItem(item)


        ## Help menu
        help = wx.MenuItem(helpMenu, ID_HELP, "plot_singlepulses Handbook")
        helpMenu.AppendItem(help)
        helpMenu.AppendSeparator()
        about = wx.MenuItem(helpMenu, ID_ABOUT, "&About")
        helpMenu.AppendItem(about)


        ## Append sub-menus to the primary menu bar
        menuBar.Append(fileMenu, "&File")
        menuBar.Append(colorMenu, "&Color")
        menuBar.Append(helpMenu, "&Help")

        self.SetMenuBar(menuBar)
        

    def init_Events(self):
        self.Bind(wx.EVT_MENU, self.onOpen,   id=ID_OPEN)
        self.Bind(wx.EVT_MENU, self.onSave,   id=ID_SAVE)
        self.Bind(wx.EVT_MENU, self.onSaveAs, id=ID_SAVE_AS)
        self.Bind(wx.EVT_MENU, self.onExit,   id=ID_EXIT)

        self.Bind(wx.EVT_MENU, self.onAbout,  id=ID_ABOUT)
        



    #####################
    ## EVENT MANAGMENT ##
    #####################
    def onOpen(self, event):
        """
        Open file to read
        """
        if self.edited:
            dlg_warn = wx.MessageDialogue("Current flagged regions have not been saved.\n\n Open new file anyway?", "Confirm Open", style=wx.YES_NO|wx.NO_DEFAULT|wx.ICON_QUESTION)
        
            if dlg_warn.ShowModal() == wx.ID_YES:
                pass
            else:
                return False

        dlg = wx.FileDialog(self, "Choose a file", self.dirname, "", "Singlepulse Files|*.singlepulse|All Files|*.*", wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.filename = dlg.GetFilename()
            self.dirname = dlg.GetDirectory()
            self.data = SinglePulse_GUI(self)

        dlg.Destroy()
        

    def onSave(self, event):
        """
        Save changes to opened file
        """


    def onSaveAs(self, event):
        """
        Save changes to new file
        """
    

    def onExit(self, event):
        """
        Quit plot_singlepulses
        """
        if self.edited:
            dlg_warn = wx.MessageDialogue("Current flagged regions have not been saved.\n\n Exit anyway?", "Confirm Exit", style=wx.YES_NO|wx.NO_DEFAULT|wx.ICON_QUESTION)
                
            if dlg_warn.ShowModal() == wx.ID_YES:
                pass
            else:
                return False

        self.Destroy()  


    def onAbout(self, event):
        """
        Display a very very very brief 'about' window.
        """

        dlg = wx.AboutDialogInfo()

        dlg.SetName('plot_singlepulses')
        dlg.SetVersion(__version__)
        dlg.SetDescription("GUI for interactively viewing PRESTO single_pulse_search.py output")
        #dialog.SetWebSite('')
        dlg.AddDeveloper(__author__)

        wx.AboutBox(dlg)





class SinglePulse_GUI(object):
    def __init__(self, frame):
        self.frame     = frame
        self.filename  = ''
        self.filenames = []

        self.ax1a = None
        self.ax1b = None
        self.ax1c = None
        self.ax2  = None
        self.cmap = matplotlib.cm.get_cmap('cubehelix')
        self.norm = matplotlib.colors.Normalize


    def load_data(self, filename):
        """
        Load in single pulse events detected from PRESTO .singlepulse files"
        """
        self.filename = filename
        print "Extracting events from {0}".format(os.path.basename(self.filename))
        ## put events into a Pandas Data Frame (database-like)
        events = pandas.read_csv(self.filename)
        
        



if __name__ == '__main__':
    print "hello..."
    app = wx.App(0)
    frame = MainWindow(None,-1)
    frame.render()
    app.MainLoop()
