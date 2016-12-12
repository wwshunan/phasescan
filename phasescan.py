#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import wx
import numpy as np
import time
import os
import threading
from epics import PV
from epics.ca import CAThread, create_context, destroy_context
import matplotlib
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas  
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg as NavigationToolbar  
from matplotlib.widgets import RectangleSelector
from matplotlib.ticker import MultipleLocator, FuncFormatter
from matplotlib.figure import Figure
from matplotlib.path import Path
import matplotlib.patches as patches
import scipy.constants as C

import pylab  
from matplotlib import pyplot 
from leastsq import getTWPhase

basedir = os.path.abspath(os.path.dirname(__file__))

class MyCustomToolbar(NavigationToolbar):
    ON_CUSTOM_UP = wx.NewId()
    ON_CUSTOM_DOWN = wx.NewId()
    ON_CUSTOM_REMOVE = wx.NewId()

    def __init__(self, plotCanvas, callback):
        NavigationToolbar.__init__(self, plotCanvas)
        self.xs = []
        self.ys = []
        self.callback = callback
        self.AddSeparator()
        up_ico = wx.Image('resources/up.png', wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        self.AddSimpleTool(self.ON_CUSTOM_UP, up_ico,
                           'shift selected points up', 'shift selected points up')
        down_ico = wx.Image('resources/down.png', wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        self.AddSimpleTool(self.ON_CUSTOM_DOWN, down_ico,
                           'shift selected points down', 'shift selected points down')
        delete_ico = wx.Image('resources/delete.png', wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        self.AddSimpleTool(self.ON_CUSTOM_REMOVE, delete_ico,
                           'remove selected points', 'remove selected points')

        wx.EVT_TOOL(self, self.ON_CUSTOM_UP, self.on_custom_up)
        wx.EVT_TOOL(self, self.ON_CUSTOM_DOWN, self.on_custom_down)
        wx.EVT_TOOL(self, self.ON_CUSTOM_REMOVE, self.on_custom_remove)
    
    def set_selected(self, xs, ys):
        self.xs = xs
        self.ys = ys

    def get_selected(self):
        return self.xs, self.ys

    def on_custom_up(self, event):
        self.ys = [y + 360 for y in self.ys]
        self.callback()

    def on_custom_down(self, event):
        self.ys = [y - 360 for y in self.ys]
        self.callback()

    def on_custom_remove(self, event):
        self.xs = []
        self.ys = []
        self.callback()

class FitThread(threading.Thread):
    def __init__(self, master, Win, distance, twPhase, fieldName, slope, x, y, EpeakFactor, start_phase, focus_mode):
        threading.Thread.__init__(self)
        self.master = master
        self.Win = Win
        self.distance = distance
        self.twPhase = twPhase
        self.fieldName = fieldName
        self.slope = slope
        self.x = x
        self.y = y
        self.start_phase = start_phase
        self.EpeakFactor = EpeakFactor
        self.focus_mode = focus_mode

    def run(self):
        rfPhase, energy_gain, amp, e, x_plot, y_plot = getTWPhase(self.x, self.y, self.Win, self.distance, self.twPhase, self.fieldName, self.start_phase, self.slope, self.EpeakFactor, self.focus_mode)
        self.clear_up(rfPhase, energy_gain, amp, x_plot, y_plot)

    def clear_up(self, rfPhase, energy_gain, amp, x_plot, y_plot):
        self.Win += energy_gain
        wx.CallAfter(self.master.display_frame.write_line, '%s\t%s\t%s' % (rfPhase, self.Win, amp))
        wx.CallAfter(self.master.updateGraph, self.master.scan_line, self.x, self.y)
        wx.CallAfter(self.master.updateGraph, self.master.fit_line, x_plot, y_plot)

class WorkThread(CAThread):
    def __init__(self, window, cavity_id, start_phase, stop_phase, phase_step, delay_before_scan, delay_read, num_read ): 
        CAThread.__init__(self)
        self.window = window
        self.cavity_id = cavity_id
        self.start_phase = start_phase
        self.stop_phase = stop_phase
        self.phase_step = phase_step
        self.delay_before_scan = delay_before_scan
        self.delay_read = delay_read
        self.num_read = num_read

        self.timeToQuit = threading.Event()
        self.timeToPause = threading.Event()
        self.timeToQuit.clear()
        self.timeToPause.clear()
        self.pause = False

    def scan(self, cavity_name, cavity_pv_set, cavity_pv_readback, bpm_pv_name):
        x = []
        y = []
        std_errors = []
       
        f = open('%s.%s' % (cavity_name, 'txt'), 'w')
        wx.CallAfter(self.window.slider.SetRange, 0, self.stop_phase - self.start_phase)
        self.cavity_pv = PV(cavity_pv_set)
        self.bpm_pv = PV(bpm_pv_name)
        self.cav_readback = PV(cavity_pv_readback)

        start_phase = self.start_phase

        while ((self.stop_phase - start_phase) * self.phase_step > 0):
            if self.timeToQuit.isSet(): 
                break
            if self.pause:
                self.timeToPause.wait()

            wx.CallAfter(self.window.slider.SetValue, start_phase - self.start_phase)
            wx.CallAfter(self.window.current.SetValue, str(start_phase))
            #if self.window.cavityList[index].startswith("buncher") or self.window.cavityList[index].startswith("cm2"):
            if cavity_name.startswith("buncher"):
                while True:
                    self.cavity_pv.put(start_phase)
                    time.sleep(1)
                    if self.cav_readback.get() and abs(int(self.cav_readback.get()) - int(start_phase)) < 5:
                        break
            #elif cavity_name.startswith("cm2"):
            #    for i in range(4):
            #        self.cavity_pv.put(start_phase)
            #        time.sleep(1)
            else:
                for i in range(2):
                    self.cavity_pv.put(start_phase)
                    time.sleep(1)

            bpm_phases = []
            for i in range(self.num_read):
                bpm_phase = self.bpm_pv.get()
                bpm_phases.append(bpm_phase)
                self.timeToQuit.wait(self.delay_read)
            
            while True:
                rms = np.std(bpm_phases)
                average = np.mean(bpm_phases)
                if rms < 5:
                    break

                abs_errs = [abs(e - average) for e in bpm_phases]
                max_abs_err = max(abs_errs)
                max_err_index = abs_errs.index(max_abs_err)
                bpm_phases.pop(max_err_index)
            
            f.write('%s\t' % start_phase)
            f.write('%s\t' % average)
            f.write('%s\n' % rms)

            x.append(start_phase)
            y.append(average)
            std_errors.append(rms)
            wx.CallAfter(self.window.updateGraph, self.window.scan_line, x, y)
            start_phase += self.phase_step

        f.close()
        return x, y, std_errors

    def data_save(self, x, y, errors):
        self.window.data_changed = True
        self.window.x = x
        self.window.y = y
        self.window.errors = errors

    #def handle_error(self, index):
    #    wx.CallAfter(self.window.handle_error, index)

    def run(self):
        create_context()
        
        cavity_pv_set = self.window.cavity_set_phase[self.cavity_id]
        cavity_pv_readback = self.window.cavity_get_phase[self.cavity_id]
        cavity_name = self.window.cavityList[self.cavity_id]
        bpm_pv_name = self.window.bpm_pv[self.cavity_id]
        x, y, std_errors = self.scan(cavity_name, cavity_pv_set, cavity_pv_readback, bpm_pv_name)
        self.data_save(x, y, std_errors)

        wx.CallAfter(self.window.reset_buttons)
        destroy_context()

class CanvasPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1)
        self.parent = parent
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.axes.margins(0.05)
        self.selected_indexes = []
        self.head = 0
        self.tail = 0
        #self.selected, = self.axes.plot([], [], 'o', ms=12, alpha=0.4,
        #                                color='yellow', visible=False)
        self.scan_line, = self.axes.plot([], [], 'o-')
        self.fit_line, = self.axes.plot([], [], marker='o')
        self.selected, = self.axes.plot([], [], 'rs', visible=False)
        #self.multi_selected, = self.axes.plot([], [], 'ro', ms=12, color='red', visible=False)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.NavigationToolbar = MyCustomToolbar(self.canvas, self.update)

        self.axes.set_autoscale_on(True)
        
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.NavigationToolbar, 0, wx.ALL, 5)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW, 5)
        self.SetSizer(self.sizer)

        #self.canvas.Bind(wx.EVT_ENTER_WINDOW, self.changeCursor)
        #self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('motion_notify_event', parent.updateStatusBar)
        #self.canvas.mpl_connect('motion_notify_event', self.drawRectangle)
        #self.canvas.mpl_connect('button_release_event', self.on_release)
        #self.canvas.mpl_connect('button_release_event', self.on_release)
        #self.canvas.mpl_connect('motion_notify_event', self.on_motion)

        self.clicked = False
        self.patch = None
        self.rs = RectangleSelector(self.axes, self.point_select,
                                    drawtype='box', useblit=True,
                                    button=[1, 3],
                                    minspanx=5, minspany=5,
                                    spancoords='pixels')
        self.canvas.draw()

    #def changeCursor(self, event):
    #    self.canvas.SetCursor(wxc.StockCursor(wx.CURSOR_CROSS))

    def point_select(self, eclick, erelease):
        self.selected_indexes = []
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata
        x_min, x_max = min(x1, x2), max(x1, x2)
        y_min, y_max = min(y1, y2), max(y1, y2)
        #parent.handle_point_select(x1, y1, x2, y2)
        line = self.scan_line
        xdata = line.get_xdata()
        ydata = line.get_ydata()

        for i, (x, y) in enumerate(zip(xdata, ydata)):
            if (x > x_min) and (x < x_max) and (y > y_min) and (y < y_max):
                self.selected_indexes.append(i)

        if self.selected_indexes:
            self.head = self.selected_indexes[0]
            self.tail = self.selected_indexes[-1]
        else:
            self.head = 0
            self.tail = -1
        
        self.selected.set_visible(True)
        self.selected.set_data(xdata[self.head:self.tail+1], ydata[self.head:self.tail+1])
        self.NavigationToolbar.set_selected(xdata[self.head:self.tail+1], ydata[self.head:self.tail+1])
        self.canvas.draw()

    def update(self):
        xs, ys = self.NavigationToolbar.get_selected()
        line = self.scan_line
        xdata = line.get_xdata()
        ydata = line.get_ydata()
        if self.selected_indexes and ys:
            xdata[self.head:self.tail+1] = xs
            ydata[self.head:self.tail+1] = ys
        else:
            xdata = [x for i, x in enumerate(xdata) if i not in self.selected_indexes]
            ydata = [y for i, y in enumerate(ydata) if i not in self.selected_indexes]
        line.set_data(xdata, ydata)
        self.selected.set_data([], [])
        self.selected_indexes = []
        self.autoscale()

    def autoscale(self):
        self.axes.relim()
        self.axes.autoscale_view()
        self.canvas.draw()


        


            
    '''
    def on_pick(self, event, line):

        if event.artist != line:
            return True

        N = len(event.ind)

        if not N:
            return True

        xs = line.get_xdata()
        ys = line.get_ydata()

        # the click locations
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata

        distances = np.hypot(x - xs[event.ind], y - ys[event.ind])

        indmin = distances.argmin()

        dataind = event.ind[indmin]

        self.lastind = dataind
        self.update(xs, ys)

    def update(self, xs, ys):
        if self.lastind is None:
            return

        dataind = self.lastind

        #ax2.cla()
        #ax2.plot(X[dataind])

        #ax2.text(0.05, 0.9, 'mu=%1.3f\nsigma=%1.3f' % (xs[dataind], ys[dataind]),
        #         transform=ax2.transAxes, va='top')
        #ax2.set_ylim(-0.5, 1.5)
        self.selected.set_visible(True)
        self.selected.set_data(xs[dataind], ys[dataind])

        self.canvas.draw()
    '''

'''
    def on_press(self, event):
        self.clicked = True
        self.x_pos = event.xdata
        self.y_pos = event.ydata

    def on_release(self, event):
        self.clicked = False
        if self.patch:
            self.patch.remove()
            self.patch = None
        self.canvas.draw()

    def drawRectangle(self, event):
        if self.clicked:
            x = event.xdata
            y = event.ydata

            verts = [
                (self.x_pos, y), 
                (self.x_pos, self.y_pos), 
                (x, self.y_pos), 
                (x, y),
                (self.x_pos, y),
                ]

            codes = [Path.MOVETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.LINETO,
                     Path.CLOSEPOLY,
                     ]

            path = Path(verts, codes)
            if self.patch:
                self.patch.remove()
            self.patch = patches.PathPatch(path, alpha=0.1)
            self.axes.add_patch(self.patch)
            self.canvas.draw()
    '''

class DisplayFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, -1, 'Display parameters for each cavity')
        panel = wx.Panel(self, -1)
        self.text = wx.TextCtrl(panel, -1, "", style=wx.TE_MULTILINE)

        bsizer = wx.BoxSizer()
        bsizer.Add(self.text, 1, wx.EXPAND)
        panel.SetSizerAndFit(bsizer)

    def write_line(self, text):
        self.text.AppendText(text)
        self.text.AppendText('\n')

class MyFrame(wx.Frame):
    cavity_set_phase = ['LLRF:Buncher1:PHA_SET', 'LLRF:Buncher2:PHA_SET', 'SCRF:CAV7:PHASE:SETPOINT', 'SCRF:CAV8:PHASE:SETPOINT', 'SCRF:CAV9:PHASE:SETPOINT', 'SCRF:CAV10:PHASE:SETPOINT', 'SCRF:CAV11:PHASE:SETPOINT', 'SCRF:CAV12:PHASE:SETPOINT', 'SCRF:CAV13:PHASE:SETPOINT', 'SCRF:CAV14:PHASE:SETPOINT', 'SCRF:CAV15:PHASE:SETPOINT', 'SCRF:CAV16:PHASE:SETPOINT', 'SCRF:CAV17:PHASE:SETPOINT', 'SCRF:CAV18:PHASE:SETPOINT']
    cavity_get_phase = ['LLRF:Buncher1:CAVITY_PHASE', 'LLRF:Buncher2:CAVITY_PHASE', 'SCRF:CAV7:PHASE:SETPOINT', 'SCRF:CAV8:PHASE:SETPOINT', 'SCRF:CAV9:PHASE:SETPOINT', 'SCRF:CAV10:PHASE:SETPOINT', 'SCRF:CAV11:PHASE:SETPOINT', 'SCRF:CAV12:PHASE:SETPOINT', 'SCRF:CAV13:PHASE:SETPOINT', 'SCRF:CAV14:PHASE:SETPOINT', 'SCRF:CAV15:PHASE:SETPOINT', 'SCRF:CAV16:PHASE:SETPOINT', 'SCRF:CAV17:PHASE:SETPOINT', 'SCRF:CAV18:PHASE:SETPOINT']
    bpm_pv = ['Bpm:2-P11', 'Bpm:5-P11', 'Bpm:6-P11', 'Bpm:7-P11', 'Bpm:8-P11', 'Bpm:9-P11', 'Bpm:10-P11', 'Bpm:11-P11', 'Bpm:11-P11', 'Bpm:12-P11', 'Bpm:13-P11', 'Bpm:14-P11', 'Bpm:15-P11', 'Bpm:17-P11']
    distance_cav_bpm = [0.1, 0.15, 0.1026, 0.1026, 0.1026, 0.1026, 0.1026, 1.3026, 0.1026, 0.1026, 0.1026, 0.1026, 0.1026, 0.8045]
    field_names = ['buncher_field.txt', 'buncher_field.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt']
    synch_phases = [-90, -90, -22, -21, -20, -18, -21, -15, -20, -16, -10, -10, 0, 0] 
    slopes = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  
    wildcard = "Phase files (*.txt)|*.txt|All files (*.*)|*.*"
    TOLERANCE = 100

    def __init__(self):
        wx.Frame.__init__(self, None, -1, "PhaseScan1")
        self.panel = wx.Panel(self, -1)
        self.pltPanel = CanvasPanel(self) 

        menuBar = wx.MenuBar()
        file_menu = wx.Menu()
        open_item = file_menu.Append(-1, "Open")
        self.Bind(wx.EVT_MENU, self.open, open_item)
        save_item = file_menu.Append(-1, "Save As")
        self.Bind(wx.EVT_MENU, self.save, save_item)
        menuBar.Append(file_menu, "File")

        graph_menu = wx.Menu()
        graph_fit_item = graph_menu.Append(-1, "Fit")
        graph_clear_item = graph_menu.Append(-1, "Clear")
        self.Bind(wx.EVT_MENU, self.on_fit, graph_fit_item)
        self.Bind(wx.EVT_MENU, self.clear, graph_clear_item)
        menuBar.Append(graph_menu, "Graph")
        self.SetMenuBar(menuBar)

        focus_list = ['Descend', 'Ascend']
        self.focus = wx.RadioBox(self.panel, -1, 'focus', wx.DefaultPosition, wx.DefaultSize, focus_list, 2, wx.RA_SPECIFY_COLS)
        self.cavityList = ['buncher1', 'buncher2', 'cm1-1', 'cm1-2', 'cm1-3', 'cm1-4', 'cm1-5', 'cm1-6', 'cm2-1', 'cm2-2', 'cm2-3', 'cm2-4', 'cm2-5', 'cm2-6']

        self.cavity_name = wx.StaticText(self.panel, -1, 'Cavity')
        #self.end_cavity_name = wx.StaticText(self.panel, -1, 'End Cavity')
        self.cavity = wx.ComboBox(self.panel, -1, 'buncher1', wx.DefaultPosition, wx.DefaultSize, self.cavityList, wx.CB_DROPDOWN)
        #self.end_cavity = wx.ComboBox(self.panel, -1, 'cm2-6', wx.DefaultPosition, wx.DefaultSize, self.cavityList, wx.CB_DROPDOWN)

        self.begin = wx.TextCtrl(self.panel, -1, '-178', size=(50, -1))
        self.current = wx.TextCtrl(self.panel, -1, '0', size=(50, -1))
        self.end = wx.TextCtrl(self.panel, -1, '180', size=(50, -1))
        self.slider = wx.Slider(self.panel, -1, 0, -178, 180, style=wx.SL_HORIZONTAL)

        self.stepLabel = wx.StaticText(self.panel, -1, 'SCAN with step:')
        self.step = wx.TextCtrl(self.panel, -1, '10', size=(50, -1))
        #self.record_label_pre = wx.StaticText(self.panel, -1, 'Fetch data every')
        #self.record_label_suf = wx.StaticText(self.panel, -1, 'steps')
        self.delayLabel = wx.StaticText(self.panel, -1, 'time delay after setting [sec]:')
        self.delay = wx.TextCtrl(self.panel, -1, '0.5', size=(50, -1))

        self.averageRadio = wx.RadioButton(self.panel, -1, 'Average for N read out with T delay')
        self.avg_num_title = wx.StaticText(self.panel, -1, 'N')
        self.avg_num = wx.TextCtrl(self.panel, -1, '4', size=(50, -1))
        self.avg_delay_title = wx.StaticText(self.panel, -1, 'T delay [sec]=')
        self.avg_delay = wx.TextCtrl(self.panel, -1, '1', size=(50, -1))
 
        self.startButton = wx.Button(self.panel, -1, "start")
        self.Bind(wx.EVT_BUTTON, self.OnStart, self.startButton)
        self.pauseButton = wx.Button(self.panel, -1, "pause")
        self.Bind(wx.EVT_BUTTON, self.OnPause, self.pauseButton)
        self.stopButton = wx.Button(self.panel, -1, "stop")
        self.Bind(wx.EVT_BUTTON, self.OnStop, self.stopButton)
        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

        injectEnergyTitle = wx.StaticText(self.panel, -1, 'Win[MeV]')
        self.injectEnergy = wx.TextCtrl(self.panel, -1, '2.1', size=(70, -1))

        self.statusBar = self.CreateStatusBar()

        cavSizer = wx.FlexGridSizer(1, 2, 5, 5)
        cavSizer.AddGrowableCol(1)
        cavSizer.AddMany([(self.cavity_name, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL), (self.cavity, 0, wx.EXPAND)])

        injectSizer = wx.BoxSizer(wx.HORIZONTAL)
        injectSizer.Add(injectEnergyTitle, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        injectSizer.Add(self.injectEnergy, 0, wx.ALL, 2)
        
        rangeSizer = wx.BoxSizer(wx.HORIZONTAL)
        rangeSizer.Add(self.begin, 1, wx.ALL, 2) 
        rangeSizer.Add(self.current, 1, wx.ALL, 2)
        rangeSizer.Add(self.end, 1, wx.ALL, 2)

        stepSizer = wx.BoxSizer(wx.HORIZONTAL)
        stepSizer.Add(self.stepLabel, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        stepSizer.Add(self.step, 0, wx.ALL, 2)

        delaySizer = wx.BoxSizer(wx.HORIZONTAL)
        delaySizer.Add(self.delayLabel, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        delaySizer.Add(self.delay, 0, wx.ALL, 2)

        avgSizer = wx.BoxSizer(wx.HORIZONTAL)
        avgSizer.Add(self.avg_num_title, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        avgSizer.Add(self.avg_num, 0, wx.ALL, 2)
        avgSizer.Add(self.avg_delay_title, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        avgSizer.Add(self.avg_delay, 0, wx.ALL, 2)
        
        btSizer = wx.BoxSizer(wx.HORIZONTAL)
        btSizer.Add(self.startButton, 0, wx.ALL, 2)
        btSizer.Add(self.pauseButton, 0, wx.ALL, 2)
        btSizer.Add(self.stopButton, 0, wx.ALL, 2)

        scanBox = wx.StaticBox(self.panel, -1, 'scan')
        sizer = wx.StaticBoxSizer(scanBox, wx.VERTICAL)
        sizer.Add(self.focus, 0, wx.EXPAND | wx.ALL, 5)
        sizer.Add(injectSizer, 0, wx.EXPAND | wx.ALL, 5)
        sizer.Add(cavSizer, 0, wx.EXPAND | wx.ALL, 5)
        sizer.Add(rangeSizer, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(self.slider, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(stepSizer, 0, wx.ALL, 5)
        sizer.Add(delaySizer, 0, wx.ALL, 5)
        sizer.Add(self.averageRadio, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(avgSizer, 0, wx.ALIGN_RIGHT | wx.ALL, 5)
        sizer.Add(btSizer, 0, wx.ALIGN_RIGHT | wx.ALL, 5)

        self.panel.SetSizer(sizer)
        
        sz = wx.BoxSizer(wx.HORIZONTAL)
        sz.Add(self.panel, 0, wx.ALL, 5)
        sz.Add(self.pltPanel, 1, wx.ALL, 5)
        self.SetSizerAndFit(sz)

        self.pauseButton.Disable()
        self.stopButton.Disable()

        self.display_frame = DisplayFrame()
        self.set_lines()
        #self.pltPanel.canvas.mpl_connect('pick_event', self.on_pick)
        #self.pltPanel.canvas.mpl_connect('pick_event', self.on_pick)
        #self.pltPanel.canvas.mpl_connect('button_press_event', self.on_press)
        #self.pltPanel.canvas.mpl_connect('button_release_event', self.on_release)
        #self.pltPanel.canvas.mpl_connect('motion_notify_event', self.on_motion)


    #def on_pick(self, event):
    #    self.pltPanel.on_pick(event, self.scan_line)

    def save(self, event):
        if self.data_changed:
            dlg = wx.FileDialog(self, "Save data as...", os.getcwd(), style=wx.SAVE | wx.OVERWRITE_PROMPT, wildcard=self.wildcard)
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                if not os.path.splitext(filename)[1]:
                    filename = filename + '.txt'
                self.save_file(filename)
            dlg.Destroy()

    def on_fit(self, event):
        x = self.scan_line.get_xdata()
        y = self.scan_line.get_ydata()
        self.fit(x, y)

    def fit(self, x, y):
        index = self.cavityList.index(self.cavity.GetValue())
        if self.cavityList[index].startswith("buncher"):
            EpeakFactor = 600
        else:
            EpeakFactor = 25
        Win = float(self.injectEnergy.GetValue())
        distance = self.distance_cav_bpm[index]
        twPhase = self.synch_phases[index]
        fieldName = self.field_names[index]
        slope = self.slopes[index]
        focus_mode = self.focus.GetSelection()
        start_phase = x[0] 
        t = FitThread(self, Win, distance, twPhase, fieldName, slope, x, y, EpeakFactor, start_phase, focus_mode) 
        t.start()

    def clear(self, event):
        self.clear_graph()

    def save_file(self, filename):
        if filename:
            f = open(filename, 'w')
            #f.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (self.distance, self.twPhase, self.fieldName, self.phaseStep, self.slope, self.EpeakFactor))
            for e in zip(self.x, self.y, self.errors):
                f.write('%s\t%s\t%s\n' % e)
            f.close()

    def read_fit(self, filename):
        if filename:
            self.clear_graph()
            f = open(filename, 'r')
            data = f.readlines()
            f.close()
            #distance, twPhase, fieldName, step, slope, EpeakFactor = data[0].strip().split()
            x = []
            y = []
            for line in data:
                line_data = line.strip().split()
                x.append(float(line_data[0]))
                y.append(float(line_data[1]))
   
            #Win = float(self.injectEnergy.GetValue())
            #first_phase = float(self.begin.GetValue())

            #rfPhase, energy_gain, amp, e, x_plot, y_plot = getTWPhase(x, y, Win, float(distance), float(twPhase), fieldName, float(step), first_phase, float(slope), float(EpeakFactor))
            self.fit(x, y)
            #self.updateGraph(self.fit_line, x_plot, y_plot)

    def open(self, event):
        dlg = wx.FileDialog(self, "Open phase file...", os.getcwd(), style=wx.OPEN, wildcard=self.wildcard)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            self.read_fit(filename)
        dlg.Destroy()

    def load_config(self):
        config_file = np.loadtxt(os.path.join(basedir, 'config.txt'))
        cavity_set_phase = config_file[:, 0]
        bpm_pv = config_file[:, 1]
        distance_cav_bpm = config_file[:, 2]
        field_names = config_file[:, 3]
        synch_phases = config_file[:, 4]
        slopes = config_file[:, 5]
        return cavity_set_phase, bpm_pv, distance_cav_bpm, field_names, synch_phases, slopes


    def get_ready(self):
        self.data_changed = False
        self.startButton.Disable()
        self.acquire_control_values()
        self.display_frame.Show()
        self.data_list_init()

    def data_list_init(self):
        self.x = []
        self.y = []
        self.errors = []

    def handle_error(self, index):
        wx.MessageBox('Some problem found', 'Error', wx.OK | wx.ICON_ERROR)
        self.injectEnergy.SetValue(self.Win)
        self.cavity.SetValue(self.cavityList[index])
        self.stop_thread()
        self.reset_buttons()

    def stop_thread(self):
        self.thread.timeToQuit.set()

    def reset_buttons(self):
        self.startButton.Enable()
        self.stopButton.Disable()
        self.pauseButton.Disable()
        self.pauseButton.SetLabel('pause')

    def acquire_control_values(self):
        self.cavity_id = self.cavityList.index(self.cavity.GetValue())
        self.start_phase = float(self.begin.GetValue())
        self.stop_phase = float(self.end.GetValue())
        self.phase_step = int(self.step.GetValue())
        self.delay_before_scan = float(self.delay.GetValue())
        self.delay_read = float(self.avg_delay.GetValue())
        self.num_read = int(self.avg_num.GetValue())
        self.focus_value = self.focus.GetSelection()
    
    def set_lines(self):
        #self.scan_line, = self.pltPanel.axes.plot([], [], marker='o', picker=5)
        self.scan_line = self.pltPanel.scan_line
        self.fit_line = self.pltPanel.fit_line
        self.selected = self.pltPanel.selected

    '''
    def resetCanvas(self):
        self.scan_line.set_xdata([])
        self.scan_line.set_ydata([])
        self.fit_line.set_xdata([])
        self.fit_line.set_ydata([])
    '''

    def OnStart(self, event):
        self.get_ready()

        if self.pauseButton.Enabled:
            self.thread.timeToPause.set()
            self.thread.timeToQuit.set()
        else:
            self.pauseButton.Enable()

        self.thread = WorkThread(self, self.cavity_id, self.start_phase, self.stop_phase, self.phase_step, self.delay_before_scan, self.delay_read, self.num_read)
        self.thread.start()

        self.stopButton.Enable()
        self.pauseButton.SetLabel('pause')

    def OnStop(self, event):
        self.reset_buttons()
        self.stop_thread()

    def OnPause(self, event):
        if self.pauseButton.GetLabel() == 'pause':
            self.thread.pause = True
            self.startButton.Enable()
            self.stopButton.Disable()
            self.pauseButton.SetLabel('resume')
        else:
            self.thread.timeToPause.set()
            self.pauseButton.SetLabel('pause')
            self.startButton.Disable()
            self.stopButton.Enable()

    def changeLabel(self, label):
        self.startButton.SetLabel(label)

    def OnCloseWindow(self, event):
        if self.stopButton.Enabled:
            self.OnStop(event)
        self.Destroy()

    def getBpmPhase(self, value):
        self.bpm_phase.SetValue(str(value))

    def updateGraph(self, line, x, y):
        plotPanel = self.pltPanel
        line.set_xdata(x)
        line.set_ydata(y)
        plotPanel.autoscale()

    def updateStatusBar(self, event):
        if event.inaxes:
            x, y = event.xdata, event.ydata
            self.statusBar.SetStatusText(('x= ' + str(x) + ' y=' + str(y)), 0)

    def clear_graph(self):
        self.updateGraph(self.scan_line, [], [])
        self.updateGraph(self.fit_line, [], [])
        self.updateGraph(self.selected, [], [])

if __name__ == '__main__':
    app = wx.App()
    frame = MyFrame()
    frame.Show(True)
    app.MainLoop()


