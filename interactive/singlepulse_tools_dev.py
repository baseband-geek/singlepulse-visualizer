#!/usr/bin/python
# DM      Sigma      Time (s)     Sample    Downfact


import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pulsar_tools import disp_delay
import math

from matplotlib.widgets import AxesWidget

class VertSlider(AxesWidget):
    """
    A slider representing a floating point range

    The following attributes are defined
      *ax*        : the slider :class:`matplotlib.axes.Axes` instance

      *val*       : the current slider value

      *vline*     : a :class:`matplotlib.lines.Line2D` instance
                     representing the initial value of the slider

      *poly*      : A :class:`matplotlib.patches.Polygon` instance
                     which is the slider knob

      *valfmt*    : the format string for formatting the slider text

      *label*     : a :class:`matplotlib.text.Text` instance
                     for the slider label

      *closedmin* : whether the slider is closed on the minimum

      *closedmax* : whether the slider is closed on the maximum

      *slidermin* : another slider - if not *None*, this slider must be
                     greater than *slidermin*

      *slidermax* : another slider - if not *None*, this slider must be
                     less than *slidermax*

      *dragging*  : allow for mouse dragging on slider

    Call :meth:`on_changed` to connect to the slider event
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, **kwargs):
        """
        Create a slider from *valmin* to *valmax* in axes *ax*

        *valinit*
            The slider initial position

        *label*
            The slider label

        *valfmt*
            Used to format the slider value

        *closedmin* and *closedmax*
            Indicate whether the slider interval is closed

        *slidermin* and *slidermax*
            Used to constrain the value of this slider to the values
            of other sliders.

        additional kwargs are passed on to ``self.poly`` which is the
        :class:`matplotlib.patches.Rectangle` which draws the slider
        knob.  See the :class:`matplotlib.patches.Rectangle` documentation
        valid property names (e.g., *facecolor*, *edgecolor*, *alpha*, ...)
        """
        AxesWidget.__init__(self, ax)

        self.valmin = valmin
        self.valmax = valmax
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axhspan(valmin, valinit, 0, 1, **kwargs)

        self.vline = ax.axhline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_xticks([])
        ax.set_ylim((valmin, valmax))
        ax.set_yticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(0.5, 1.03, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='center')

        self.valtext = ax.text(0.5, -0.03, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='center')

        self.cnt = 0
        self.observers = {}

        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return

        val = event.ydata
        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val

        self.set_val(val)

    def set_val(self, val):
        xy = self.poly.xy
        xy[1] = 0, val
        xy[2] = 1, val
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw()
        self.val = val
        if not self.eventson:
            return
        for cid, func in self.observers.iteritems():
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed, call *func* with the new
        slider position

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """reset the slider to the initial value if needed"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)

def loadfile(filename):
    if filename==None:
        print "No filename supplied to read..."

    with open(filename) as f:
        data = f.read()

    data = data.split('\n')

    while not data[-1]:
        data = data[:-1]
    return data

def obs_stats(time, flags):
    # Not doing total time correctly, depends on last single pulse detection instead of observation time
    flag_time = 0

    for flag in flags:
        flag_time+=(float(flag.split()[1])-float(flag.split()[0]))
    print "%.2f seconds flagged from %.2f seconds of data (%.2f percent)" % ( flag_time, time[-1], flag_time/time[-1]*100)


def flagfile(basename, max_DM=2097.2, freq_l=0.169615, freq_h=0.200335, padding=3):
    """This function takes in a text file of bad 0 DM times and
    writes out one flagged over the correct de-dispersive smearing
    times, looking for overlaps along the way. There must be a text file named
    basename.bad with rows indicating bad times for this to work. """
    from subprocess import Popen


    bads = np.genfromtxt(basename+'.bad')
    i = 0 # initialize counter for new list
    flags = []
    for bad in bads:
        start = bad[0] - (padding + disp_delay(freq1=freq_l, freq2=freq_h, DM=max_DM)/1000)
        if start < 0:
            start = 0
        stop = bad[1] + padding
        if len(flags) > 0:
            if start <= flags[-1][1]:
                flags[-1][1] = stop
            else:
                flags.append([start, stop])
        else:
            flags.append([start, stop])
    np.savetxt(basename+'.flag', flags, fmt='%d')
    Popen(['flag.sh', basename]).communicate()[0]

def max_nth_percent(n, data):
    """ A function that returns the nth percent top value, planned use is for plotting
    :param n: the percentile value desired for return
    :param data: the iterable object searched through
    :return: nth percent largest value
    """
    import heapq

    data=list(data)
    n=float(n)

    return heapq.nlargest(int(len(data)*(n/100.0)), data)[-1]


def singlepulse_plot(basename=None, DMvTime=1, StatPlots=False, raw = False, threshold=5.0, Downsamps=[20, 20, 20], colormap='cool'	):
    """Plots up the flagged data, should switch to using genfromtxt when I have the time"""

    from matplotlib.widgets import RadioButtons, Slider

#    if raw:
    data_raw=loadfile(basename+'.singlepulse')[1:]
#        flags = False
#    else:
    data_flag = loadfile(basename+'_flagged.singlepulse')[1:]
    flags = loadfile(basename+'.flag')


    DM_raw = [float(row.split()[0]) for row in data_raw if float(row.split()[1]) >= threshold]
    Sigma_raw = [float(row.split()[1]) for row in data_raw if float(row.split()[1]) >= threshold]
    Time_raw = [float(row.split()[2]) for row in data_raw if float(row.split()[1]) >= threshold]
    Sample_raw = [int(row.split()[3]) for row in data_raw if float(row.split()[1]) >= threshold]
    Downfact_raw = [int(row.split()[4]) for row in data_raw if float(row.split()[1]) >= threshold]

    DM_flag = [float(row.split()[0]) for row in data_flag if float(row.split()[1]) >= threshold]
    Sigma_flag = [float(row.split()[1]) for row in data_flag if float(row.split()[1]) >= threshold]
    Time_flag = [float(row.split()[2]) for row in data_flag if float(row.split()[1]) >= threshold]
    Sample_flag = [int(row.split()[3]) for row in data_flag if float(row.split()[1]) >= threshold]
    Downfact_flag = [int(row.split()[4]) for row in data_flag if float(row.split()[1]) >= threshold]



    Sigma_float_raw = [float(value) for value in Sigma_raw]
    Downfact_float_raw = [float(value) for value in Downfact_raw]

    Sigma_float_flag = [float(value) for value in Sigma_flag]
    Downfact_float_flag = [float(value) for value in Downfact_flag]

    DM = DM_flag
    Sigma = Sigma_flag
    Time = Time_flag
    Sample = Sample_flag
    Downfact = Downfact_flag

    Sigma_float = Sigma_float_flag
    Downfact_float = Downfact_float_flag

    fig = plt.figure()

    cm = plt.cm.get_cmap(colormap)

    def onpick(event):
        points = event.artist
        ind = event.ind
        mouseevent = event.mouseevent
        print '\n'
        print "Information for data points around click event %.4f, %.4f:" % (mouseevent.xdata,  mouseevent.ydata)
        for i in ind:
            if ( DM[i] < 150):
                boxcar = Downfact[i]
            elif ( 150<= DM[i] < 823.2 ):
                boxcar = Downfact[i] * Downsamps[0]/10
            elif ( 823.2 <= DM[i] < 1486.2):
                boxcar = Downfact[i] * Downsamps[1]/10
            elif ( 1426.2 <= DM[i] < 2100):
                boxcar = Downfact[i] * Downsamps[2]/10
            print "%.2f seconds, %.2f Sigma event detected at a DM of %.2f with a boxcar of: %d ms" % (Time[i], Sigma[i], DM[i], boxcar)


    if StatPlots:
        ax0 = fig.add_subplot(231)
        plt.hist(Sigma, histtype='step', bins=60)
        ax0.set_xlabel('Signal-to-Noise')
        ax0.set_ylabel('Number of Pulses')

        ax1 = fig.add_subplot(232)
        plt.hist(DM, histtype='step', bins=int(0.5*len(set(DM))))
        ax1.set_xlabel('DM (pc cm^-3)')
        ax1.set_ylabel('Number of Pulses')

        ax2 = fig.add_subplot(233)
        sc0=ax2.scatter(DM, Sigma, c=Downfact_float, vmin=min(Downfact_float), vmax=max_nth_percent(10,Downfact_float), cmap='spectral', picker=1)
        ax2.set_ylabel('Signal-to-Noise')
        ax2.set_xlabel('DM (pc cm^-3)')
        ax2.set_xlim([min(DM), max(DM)])
        ax2.set_ylim([min(Sigma), max(Sigma)])
 #       plt.colorbar(sc0, label="Boxcar")
#        fig.canvas.mpl_connect('pick_event', onpick)

        ax3 = fig.add_subplot(212)
    else:
        ax3 = fig.add_subplot(111)

#	ax3.set_title("Single Pulse Sigma")    
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('DM (pc cm^-3)')
    ax3.set_ylim([min(DM), max(DM)])
    ax3.set_xlim([min(Time), max(Time)])
 #   cm = plt.cm.get_cmap(colormap)
    sc=ax3.scatter(Time,DM, c=Sigma_float, vmin=min(Sigma_float), vmax=max_nth_percent(10,Sigma_float), cmap=cm, picker=1)
#	leg = ax1.legend()
    plt.colorbar(sc, label="Signal-to-Noise")


    rax = plt.axes([0.85, 0.02, 0.13, 0.13])
    radio = RadioButtons(rax, ('Flagged', 'Raw'))
    def datafunc(label):
        plot_argument_dict= {'Flagged':"Time_flag, DM_flag,c=Sigma_float_flag, vmin=min(Sigma_float_flag), vmax=max(Sigma_float_flag)"  , 'Raw':"Time_raw,DM_raw,c=Sigma_float_raw, vmin=min(Sigma_float_raw), vmax=max(Sigma_float_raw)"}
        sc=ax3.scatter(plot_argument_dict[label], cmap=cm, picker=1)
#        ydata = hzdict[label]
#        ax3.set_ydata(ydata)
        plt.draw()
    radio.on_clicked(datafunc)
    axmin = plt.axes([0.85, 0.2, 0.04, 0.3])
    axmax = plt.axes([0.93, 0.2, 0.04, 0.3])
    smin = VertSlider(axmin, 'Min', min(Sigma_float), max(Sigma_float), valinit=min(Sigma_float))
    smax = VertSlider(axmax, 'Max', min(Sigma_float), max(Sigma_float), valinit=max(Sigma_float))
    def update(val):
        sc.vlim([smin.val,smax.val])
        fig.canvas.draw()
#    smin.on_changed(update)
    smax.on_changed(update)
    if not raw:
        for flag in flags:
            flag_area = matplotlib.patches.Rectangle((float(flag.split()[0]), min(DM)), float(flag.split()[1])-float(flag.split()[0]), max(DM)-min(DM),  edgecolor='0', facecolor='0.66')
            ax3.add_patch(flag_area)


    fig.canvas.mpl_connect('pick_event', onpick)

    '''
    ax2 = fig.add_subplot(122)
    ax2.set_title("Single Pulse Boxcar")
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('DM (pc cm^-3)')
    cm = plt.cm.get_cmap('RdYlBu')
    sc2=ax2.scatter(Time,DM, c=Downfact_float, vmin=min(Downfact_float), vmax=max(Downfact_float), cmap=cm)
#	leg = ax1.legend()
    plt.colorbar(sc2)
    if not raw:
        for flag in flags:
            flag_area = matplotlib.patches.Rectangle((float(flag.split()[0]), min(DM)), float(flag.split()[1])-float(flag.split()[0]), max(DM)-min(DM),  edgecolor='0', facecolor='0.66')
            ax2.add_patch(flag_area)
    '''
    fig.suptitle('Single Pulse Search results for '+basename)
    plt.show()

    obs_stats(Time, flags)


def slice(infile, dm=None, timerange=None, sigma=None, downfact=None):
    # Not properly implemented yet

    data = read_singlepulse(infile)

    slices = [None]*5

    slice_map = {'dm':0, 'sigma':1, 'timerange':2, 'sample':3, 'downfact':4}

    DM = [row.split()[0] for row in data]
    Sigma = [row.split()[1] for row in data]
    Time = [row.split()[2] for row in data]
    Sample = [row.split()[3] for row in data]
    Downfact = [row.split()[4] for row in data]

    if dm:
        if type(dm) == type(0) or type(0.0):
            data = [row   for row in data if dm <= row.split()[0]]
        elif type(dm) == type([]):
            data = [row  for row in data if dm[0] <= row.split()[0] <= dm[1]]
    if sigma:
        if type(sigma) == type(0) or type(0.0):
            data = [row for row in data if sigma <= row.split()[1]  ]
        elif type(sigma) == type([]):
            data = [row for row in data if sigma[0] <= row.split()[1] <= sigma[1]]




