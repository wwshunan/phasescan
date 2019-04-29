#!/usr/bin/env python
# -*- coding: utf-8 -*- 

import wx
import numpy as np
import time
import os
import re
import threading
from epics import PV, caput
from epics.ca import CAThread, create_context, destroy_context, use_initial_context
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
    ON_CUSTOM_FIT = wx.NewId()

    def __init__(self, plotCanvas, callback, fit):
        NavigationToolbar.__init__(self, plotCanvas)
        self.parent = plotCanvas
        self.xs = []
        self.ys = []
        self.callback = callback
        self.fit = fit
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
        fit_ico = wx.Image('resources/fit.png', wx.BITMAP_TYPE_ANY).ConvertToBitmap()
        self.AddSimpleTool(self.ON_CUSTOM_FIT, fit_ico,
                           'data fit', 'data fit')

        wx.EVT_TOOL(self, self.ON_CUSTOM_UP, self.on_custom_up)
        wx.EVT_TOOL(self, self.ON_CUSTOM_DOWN, self.on_custom_down)
        wx.EVT_TOOL(self, self.ON_CUSTOM_REMOVE, self.on_custom_remove)
        wx.EVT_TOOL(self, self.ON_CUSTOM_FIT, self.on_custom_fit)
    
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

    def on_custom_fit(self, event):
        self.fit()

class FitThread(CAThread):
    def __init__(self, master, Win, distance, twPhase, fieldName, slope, x, y, EpeakFactor, start_phase, focus_mode, f, bpm_harm, bpm_polarity, rf_phase_pv, mass, charge, cavity_name):
        CAThread.__init__(self)
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
        self.f = f
        self.bpm_harm = bpm_harm
        self.bpm_polarity = bpm_polarity
        self.rf_phase_pv = rf_phase_pv
        self.mass = mass
        self.charge = charge
        self.cavity_name = cavity_name

    def run(self):
        rfPhase, energy_gain, amp, e, x_plot, y_plot = getTWPhase(self.x, self.y, self.Win, self.distance, self.twPhase,
                                                                  self.fieldName, self.start_phase, self.slope, self.EpeakFactor, 
                                                                  self.focus_mode, self.f, self.bpm_harm, self.bpm_polarity, 
                                                                  self.mass, self.charge)
        self.clear_up(rfPhase, energy_gain, amp, x_plot, y_plot)

    def clear_up(self, rfPhase, energy_gain, amp, x_plot, y_plot):
        self.Win += energy_gain
        wx.CallAfter(self.master.display.write_line, '{0}\nPhase: {1:.1f}\nEnergy: {2:.3f}\nEpeak: {3:.2f}\n'.format(self.cavity_name, rfPhase, self.Win, amp))
        wx.CallAfter(self.master.updateGraph, self.master.scan_line, self.x, self.y)
        wx.CallAfter(self.master.updateGraph, self.master.fit_line, x_plot, y_plot)
        with open('fit-results.txt', 'a') as f:
            f.write('{0}\t{1}\t{2}\t{3}\n'.format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), self.Win, abs(amp), rfPhase))
        #create_context()
        #if self.rf_phase_pv is not None:
        #    phase_set_pv = PV(self.rf_phase_pv)
        #   phase_set_pv.put(round(rfPhase, 1))
        #destroy_context()

class WorkThread(CAThread):
    def __init__(self, window, cavity_id, start_phase, stop_phase, phase_step, delay_before_scan, delay_read, num_read): 
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

    def scan(self, cavity_name, cavity_pv_set, cavity_pv_readback, bpm_pv_name, cavity_ready):
        x = []
        y = []
        std_errors = []
       
        f = open('%s.%s' % (cavity_name, 'txt'), 'w')
        wx.CallAfter(self.window.slider.SetRange, 0, self.stop_phase - self.start_phase)
        self.cavity_pv = PV(cavity_pv_set)
        if cavity_pv_set.startswith('CM4'):
            cavity_pv_2_str = cavity_pv_set.replace('Setpoint', 'Openset')
            cavity_pv_2 = PV(cavity_pv_2_str)
        self.bpm_pv = PV(bpm_pv_name)
        self.cav_readback = PV(cavity_pv_readback)

        if cavity_name.startswith('cm4'):   #for cm4
            self.start_phase = self.cavity_pv.get() * 360. / 4294967296
            self.stop_phase = self.start_phase + 360
            wx.CallAfter(self.window.slider.SetRange, 0, self.stop_phase - self.start_phase)

        start_phase = self.start_phase    

        while ((self.stop_phase - start_phase) * self.phase_step > 0):
            if self.timeToQuit.isSet(): 
                break
            if self.pause:
                self.timeToPause.wait()

            while True:
                read_values = [PV(x).get() for x in cavity_ready]
                all_ready = all(x > 0.5 for x in read_values)
                if all_ready:
                    break
                time.sleep(1)

            wx.CallAfter(self.window.slider.SetValue, start_phase - self.start_phase)
            wx.CallAfter(self.window.current.SetValue, str(start_phase))
            #if self.window.cavityList[index].startswith("buncher") or self.window.cavityList[index].startswith("cm2"):
            if cavity_name.startswith("buncher"):
                while True:
                    self.cavity_pv.put(start_phase)
                    time.sleep(1)
                    if self.cav_readback.get() and abs(int(self.cav_readback.get()) - int(start_phase)) < 5:
                        break
            elif cavity_name.startswith("cm4"):
                for i in range(1):
                    self.cavity_pv.put(start_phase * 4294967296 / 360.)
                    cavity_pv_2.put(cavity_pv_2.get() - self.phase_step * 4294967296 / 360.)
                    time.sleep(2)
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
        #create_context()
        use_initial_context()
        
        cavity_pv_set = self.window.cavity_set_phase[self.cavity_id]
        cavity_pv_readback = self.window.cavity_get_phase[self.cavity_id]
        cavity_name = self.window.cavityList[self.cavity_id]
        bpm_pv_name = self.window.bpm_pv[self.cavity_id]
        cavity_ready = [self.window.cavity_ready[i] for i in range(self.cavity_id+1)]
        x, y, std_errors = self.scan(cavity_name, cavity_pv_set, cavity_pv_readback, bpm_pv_name, cavity_read)
        self.data_save(x, y, std_errors)

        wx.CallAfter(self.window.reset_buttons)
        #destroy_context()

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
        self.NavigationToolbar = MyCustomToolbar(self.canvas, self.update, self.parent.button_fit)

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
        #self.axes.relim()
        #self.axes.autoscale_view()
        xdata = self.scan_line.get_xdata()
        ydata = self.scan_line.get_ydata()
        if len(xdata) > 1:
            x_min = min(xdata)
            x_max = max(xdata)
            y_min = min(ydata)
            y_max = max(ydata)
            x_width = x_max - x_min
            y_width = y_max - y_min
            x_min = x_min - 0.1 * x_width
            y_min = y_min - 0.1 * y_width
            x_max = x_max + 0.1 * x_width
            y_max = y_max + 0.1 * y_width
            self.axes.set_xlim(x_min, x_max)
            self.axes.set_ylim(y_min, y_max)
        else:
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

class Display(wx.TextCtrl):
    def __init__(self, panel):
        wx.TextCtrl.__init__(self, panel, -1, "", style=wx.TE_MULTILINE, size=(-1, 80))

    def write_line(self, text):
        self.AppendText(text)
        self.AppendText('\n')

class MyFrame(wx.Frame):
    cavity_ready = ['LLRF:Buncher1:Phase_unlock', 'LLRF:Buncher2:Phase_unlock']
    cavity_ready += ['SCRF:CAV1:READY', 'SCRF:CAV2:READY', 'SCRF:CAV3:READY']
    cavity_ready += ['SCRF:CAV4:READY', 'SCRF:CAV5:READY', 'SCRF:CAV6:READY']
    cavity_ready += ['SCRF:CAV7:READY', 'SCRF:CAV8:READY', 'SCRF:CAV9:READY']
    cavity_ready += ['SCRF:CAV10:READY', 'SCRF:CAV11:READY', 'SCRF:CAV12:READY']
    cavity_ready += ['SCRF:CAV13:READY', 'SCRF:CAV14:READY', 'SCRF:CAV15:READY']
    cavity_ready += ['SCRF:CAV16:READY', 'SCRF:CAV17:READY']
    cavity_ready += ['CM4:LLRF:CAV01:LLRF_OK', 'CM4:LLRF:CAV02:LLRF_OK', 'CM4:LLRF:CAV03:LLRF_OK']
    cavity_ready += ['CM4:LLRF:CAV04:LLRF_OK', 'CM4:LLRF:CAV05:LLRF_OK', 'CM4:LLRF:CAV06:LLRF_OK']

    cavity_set_phase = ['LLRF:Buncher1:PHA_SET', 'LLRF:Buncher2:PHA_SET'] 
    cavity_set_phase += ['SCRF:CAV1:PHASE:SETPOINT', 'SCRF:CAV2:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV3:PHASE:SETPOINT', 'SCRF:CAV4:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV5:PHASE:SETPOINT', 'SCRF:CAV6:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV7:PHASE:SETPOINT', 'SCRF:CAV8:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV9:PHASE:SETPOINT', 'SCRF:CAV10:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV11:PHASE:SETPOINT', 'SCRF:CAV12:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV13:PHASE:SETPOINT', 'SCRF:CAV14:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV15:PHASE:SETPOINT', 'SCRF:CAV16:PHASE:SETPOINT'] 
    cavity_set_phase += ['SCRF:CAV17:PHASE:SETPOINT', 'CM4:LLRF:CAV01:Loop_B_Setpoint'] 
    cavity_set_phase += ['CM4:LLRF:CAV02:Loop_B_Setpoint', 'CM4:LLRF:CAV03:Loop_B_Setpoint'] 
    cavity_set_phase += ['CM4:LLRF:CAV04:Loop_B_Setpoint', 'CM4:LLRF:CAV05:Loop_B_Setpoint'] 
    cavity_set_phase += ['CM4:LLRF:CAV06:Loop_B_Setpoint']

    cavity_get_phase = ['LLRF:Buncher1:CAVITY_PHASE', 'LLRF:Buncher2:CAVITY_PHASE'] 
    cavity_get_phase += ['SCRF:CAV1:PHASE:READ', 'SCRF:CAV2:PHASE:READ'] 
    cavity_get_phase += ['SCRF:CAV3:PHASE:READ', 'SCRF:CAV4:PHASE:READ'] 
    cavity_get_phase += ['SCRF:CAV5:PHASE:READ', 'SCRF:CAV6:PHASE:READ'] 
    cavity_get_phase += ['SCRF:CAV7:PHASE:READ', 'SCRF:CAV8:PHASE:READ'] 
    cavity_get_phase += ['SCRF:CAV9:PHASE:READ', 'SCRF:CAV10:PHASE:READ'] 
    cavity_get_phase += ['SCRF:CAV11:PHASE:READ', 'SCRF:CAV12:PHASE:READ']
    cavity_get_phase += ['SCRF:CAV13:PHASE:READ', 'SCRF:CAV14:PHASE:READ'] 
    cavity_get_phase += ['SCRF:CAV15:PHASE:READ', 'SCRF:CAV16:PHASE:READ'] 
    cavity_get_phase += ['SCRF:CAV17:PHASE:READ', 'CM4:LLRF:CAV01:Loop_B_Setpoint'] 
    cavity_get_phase += ['CM4:LLRF:CAV02:Loop_B_Setpoint', 'CM4:LLRF:CAV03:Loop_B_Setpoint']
    cavity_get_phase += ['CM4:LLRF:CAV04:Loop_B_Setpoint', 'CM4:LLRF:CAV05:Loop_B_Setpoint'] 
    cavity_get_phase += ['CM4:LLRF:CAV06:Loop_B_Setpoint']
    #bpm_pv = ['Bpm:2-P11', 'Bpm:5-P11', 'Bpm:6-P11', 'Bpm:7-P11', 'Bpm:8-P11', 'Bpm:9-P11', 'Bpm:10-P11', 'Bpm:11-P11',
    #          'Bpm:11-P11', 'Bpm:12-P11', 'Bpm:13-P11', 'Bpm:14-P11', 'Bpm:15-P11', 'Bpm:16-P11',  'Bpm:16-P11', 
    #          'Bpm:17-P11', 'Bpm:18-P11', 'Bpm:19-P11', 'BI:BPM:COLDBPM01:phase:ai', 'BI:BPM:COLDBPM01:phase:ai', 'BI:BPM:COLDBPM02:phase:ai', 'BI:BPM:COLDBPM03:phase:ai',
    #          'BI:BPM:COLDBPM04:phase:ai', 'BI:BPM:COLDBPM05:phase:ai', 'Bpm:25-P11']
    bpm_pv = ['Bpm:2-P11', 'Bpm:5-P11'] 
    bpm_pv += ['Bpm:6-P11', 'Bpm:7-P11', 'Bpm:8-P11', 'Bpm:9-P11', 'Bpm:10-P11', 'Bpm:11-P11'] 
    bpm_pv += ['Bpm:11-P11', 'Bpm:12-P11', 'Bpm:13-P11', 'Bpm:14-P11', 'Bpm:15-P11', 'Bpm:16-P11']  
    bpm_pv += ['Bpm:16-P11',  'Bpm:17-P11', 'Bpm:18-P11', 'Bpm:19-P11', 'Bpm:20-P11'] 
    bpm_pv += ['Bpm:21-P11', 'Bpm:22-P11', 'Bpm:23-P11', 'Bpm:24-P11', 'Bpm:25-P11', 'Bpm:25-P11']

    #distance_cav_bpm = [0.08134, 0.12984, 0.1026, 0.1026, 0.1026, 0.1026, 0.1026, 1.2426, 0.1026, 0.1026, 0.1026, 0.1026,
    #                    0.1026, 1.3976, 0.1126, 0.1126, 0.1126, 0.1126, 1.3992, 0.795, 0.795, 0.795, 0.795, 1.30655, 0.61855]
    distance_cav_bpm = [0.08134, 0.12984] 
    distance_cav_bpm += [0.1026, 0.1026, 0.1026, 0.1026, 0.1026, 1.2426] 
    distance_cav_bpm += [0.1026, 0.1026, 0.1026, 0.1026, 0.1026, 1.3976] 
    distance_cav_bpm += [0.1017, 0.1027, 0.1126, 0.1048, 1.3992] 
    distance_cav_bpm += [0.795, 0.795, 0.795, 0.795, 1.30655, 0.61855]

    field_names = ['buncher_field.txt', 'buncher_field.txt'] 
    field_names += ['Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt']
    field_names += ['Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt', 'Exyz.txt'] 
    field_names += ['hwr015.txt', 'hwr015.txt', 'hwr015.txt', 'hwr015.txt', 'hwr015.txt'] 
    field_names += ['spoke021.txt', 'spoke021.txt', 'spoke021.txt', 'spoke021.txt', 'spoke021.txt', 'spoke021.txt']
 
    slopes = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]  
    wildcard = "Phase files (*.txt)|*.txt|All files (*.*)|*.*"
    TOLERANCE = 100
    particle_mass_charge = {
        'Proton': {
            'mass': 938.272083,
            'charge': 1
        },
        '3He2+': {
            'mass': 2814.816,
            'charge': 2
        },
        '4He2+': {
            'mass': 3753.07,
            'charge': 2
        },       
    }            

    def __init__(self):
        wx.Frame.__init__(self, None, -1, "PhaseScan")
        self.panel = wx.Panel(self, -1)
        self.pltPanel = CanvasPanel(self) 

        menuBar = wx.MenuBar()
        file_menu = wx.Menu()
        open_item = file_menu.Append(-1, "Open")
        self.Bind(wx.EVT_MENU, self.open, open_item)
        load_item = file_menu.Append(-1, "Load Phase")
        self.Bind(wx.EVT_MENU, self.load, load_item)
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

        focus_list = ['cm1-3', 'cm4']
        self.focus = wx.RadioBox(self.panel, -1, 'focus', wx.DefaultPosition, wx.DefaultSize, focus_list, 2, wx.RA_SPECIFY_COLS)
        self.cavityList = ['buncher1', 'buncher2', 'cm1-1', 'cm1-2', 'cm1-3', 'cm1-4', 'cm1-5', 'cm1-6'] 
        self.cavityList += ['cm2-1', 'cm2-2', 'cm2-3', 'cm2-4', 'cm2-5', 'cm2-6'] 
        self.cavityList += ['cm3-1', 'cm3-2', 'cm3-3', 'cm3-4', 'cm3-5'] 
        self.cavityList += ['cm4-1', 'cm4-2', 'cm4-3', 'cm4-4', 'cm4-5', 'cm4-6']

        self.particleType = ['Proton', '3He2+', '4He2+']
        self.typeName = wx.StaticText(self.panel, -1, 'Type')
        self.type = wx.ComboBox(self.panel, -1, 'Proton', wx.DefaultPosition, wx.DefaultSize, self.particleType, wx.CB_DROPDOWN)

        self.cavity_name = wx.StaticText(self.panel, -1, 'Cavity')
        #self.end_cavity_name = wx.StaticText(self.panel, -1, 'End Cavity')
        self.cavity = wx.ComboBox(self.panel, -1, 'buncher1', wx.DefaultPosition, wx.DefaultSize, self.cavityList, wx.CB_DROPDOWN)
        #self.end_cavity = wx.ComboBox(self.panel, -1, 'cm2-6', wx.DefaultPosition, wx.DefaultSize, self.cavityList, wx.CB_DROPDOWN)

        self.begin = wx.TextCtrl(self.panel, -1, '-178', size=(50, -1))
        self.current = wx.TextCtrl(self.panel, -1, '0', size=(50, -1))
        self.end = wx.TextCtrl(self.panel, -1, '180', size=(50, -1))
        self.slider = wx.Slider(self.panel, -1, 0, -178, 180, style=wx.SL_HORIZONTAL)

        self.stepLabel = wx.StaticText(self.panel, -1, 'SCAN with step:')
        self.step = wx.TextCtrl(self.panel, -1, '20', size=(50, -1))
        #self.record_label_pre = wx.StaticText(self.panel, -1, 'Fetch data every')
        #self.record_label_suf = wx.StaticText(self.panel, -1, 'steps')
        self.delayLabel = wx.StaticText(self.panel, -1, 'time delay after setting [sec]:')
        self.delay = wx.TextCtrl(self.panel, -1, '0.5', size=(50, -1))

        self.averageRadio = wx.RadioButton(self.panel, -1, 'Average for N read out with T delay')
        self.avg_num_title = wx.StaticText(self.panel, -1, 'N')
        self.avg_num = wx.TextCtrl(self.panel, -1, '4', size=(50, -1))
        self.avg_delay_title = wx.StaticText(self.panel, -1, 'T delay [sec]=')
        self.avg_delay = wx.TextCtrl(self.panel, -1, '1', size=(50, -1))

        setPhaseLabel = wx.StaticText(self.panel, -1, 'Set phase')
        self.setPhase = wx.TextCtrl(self.panel, -1, '0', style=wx.TE_PROCESS_ENTER, size=(150, -1))
        rbPhaseLabel = wx.StaticText(self.panel, -1, 'Readback phase')
        self.rbPhase = wx.StaticText(self.panel, -1, '0', style=wx.ALIGN_CENTRE)
        self.rbPhase.SetForegroundColour((255, 0, 0))

        self.bpm_phase_label = wx.StaticText(self.panel, -1, 'BPM phase')
        self.bpm_phase = wx.StaticText(self.panel, -1, '0', style=wx.ALIGN_CENTRE)
        self.bpm_phase.SetForegroundColour((255, 0, 0))
        self.bpm_button = wx.Button(self.panel, -1, 'Get')
        self.Bind(wx.EVT_BUTTON, self.OnGetBpmPhase, self.bpm_button)

        bpmPhaseSizer = wx.BoxSizer(wx.HORIZONTAL)
        bpmPhaseSizer.Add(self.bpm_phase_label, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        bpmPhaseSizer.Add(self.bpm_phase, 1, wx.ALIGN_CENTRE | wx.ALL, 2)
        bpmPhaseSizer.Add(self.bpm_button, 0, wx.ALL, 2)
 
        self.startButton = wx.Button(self.panel, -1, "start")
        self.Bind(wx.EVT_BUTTON, self.OnStart, self.startButton)
        self.pauseButton = wx.Button(self.panel, -1, "pause")
        self.Bind(wx.EVT_BUTTON, self.OnPause, self.pauseButton)
        self.stopButton = wx.Button(self.panel, -1, "stop")
        self.Bind(wx.EVT_BUTTON, self.OnStop, self.stopButton)
        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

        injectEnergyTitle = wx.StaticText(self.panel, -1, 'Win[MeV]')
        self.injectEnergy = wx.TextCtrl(self.panel, -1, '2.12', size=(150, -1))

        #massTitle = wx.StaticText(self.panel, -1, 'Mass[MeV]')
        #self.mass = wx.TextCtrl(self.panel, -1, '938.272083', size=(150, -1))

        #qTitle = wx.StaticText(self.panel, -1, 'Q')
        #self.q = wx.TextCtrl(self.panel, -1, '1', size=(150, -1))

        self.statusBar = self.CreateStatusBar()

        typeSizer = wx.FlexGridSizer(1, 2, 5, 5)
        typeSizer.AddGrowableCol(1)
        typeSizer.AddMany([(self.typeName, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL), (self.type, 0, wx.EXPAND)])

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

        setPhaseSizer = wx.BoxSizer(wx.HORIZONTAL)
        setPhaseSizer.Add(setPhaseLabel, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        setPhaseSizer.Add(self.setPhase, 0, wx.ALL, 2)
     
        rbPhaseSizer = wx.BoxSizer(wx.HORIZONTAL)
        rbPhaseSizer.Add(rbPhaseLabel, 0, wx.ALIGN_CENTRE | wx.ALL, 2)
        rbPhaseSizer.Add(self.rbPhase, 1, wx.ALIGN_CENTRE | wx.ALL, 2)

        btSizer = wx.BoxSizer(wx.HORIZONTAL)
        btSizer.Add(self.startButton, 0, wx.ALL, 2)
        btSizer.Add(self.pauseButton, 0, wx.ALL, 2)
        btSizer.Add(self.stopButton, 0, wx.ALL, 2)

        scanBox = wx.StaticBox(self.panel, -1, 'scan')
        sizer = wx.StaticBoxSizer(scanBox, wx.VERTICAL)
        sizer.Add(self.focus, 0, wx.EXPAND | wx.ALL, 5)
        sizer.Add(typeSizer, 0, wx.EXPAND | wx.ALL, 5)
        sizer.Add(injectSizer, 0, wx.EXPAND | wx.ALL, 5)
        sizer.Add(cavSizer, 0, wx.EXPAND | wx.ALL, 5)
        sizer.Add(rangeSizer, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(self.slider, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(stepSizer, 0, wx.ALL, 5)
        sizer.Add(delaySizer, 0, wx.ALL, 5)
        sizer.Add(self.averageRadio, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(avgSizer, 0, wx.ALIGN_RIGHT | wx.ALL, 5)
        sizer.Add(setPhaseSizer, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(rbPhaseSizer, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(bpmPhaseSizer, 0, wx.EXPAND, wx.ALL, 5)
        sizer.Add(btSizer, 0, wx.ALIGN_RIGHT | wx.ALL, 5)

        self.display = Display(self.panel)
        displayBox = wx.StaticBox(self.panel, -1, 'Result Show')
        displaySizer = wx.StaticBoxSizer(displayBox, wx.VERTICAL)
        displaySizer.Add(self.display, 0, wx.EXPAND | wx.ALL, 5)

        left_sizer = wx.BoxSizer(wx.VERTICAL)
        left_sizer.Add(sizer, 0, wx.EXPAND)
        left_sizer.Add(displaySizer, 0, wx.EXPAND)

        self.panel.SetSizer(left_sizer)
        
        sz = wx.BoxSizer(wx.HORIZONTAL)
        sz.Add(self.panel, 0, wx.ALL, 5)
        sz.Add(self.pltPanel, 1, wx.ALL, 5)
        self.SetSizerAndFit(sz)

        self.pauseButton.Disable()
        self.stopButton.Disable()

        self.set_lines()
        self.setPhase.Bind(wx.EVT_TEXT_ENTER, self.setCavPhase)
        #self.pltPanel.canvas.mpl_connect('pick_event', self.on_pick)
        #self.pltPanel.canvas.mpl_connect('pick_event', self.on_pick)
        #self.pltPanel.canvas.mpl_connect('button_press_event', self.on_press)
        #self.pltPanel.canvas.mpl_connect('button_release_event', self.on_release)
        #self.pltPanel.canvas.mpl_connect('motion_notify_event', self.on_motion)


    #def on_pick(self, event):
    #    self.pltPanel.on_pick(event, self.scan_line)
    
    def setCavPhase(self, event):
        rfPhase = float(self.setPhase.GetValue())
        index = self.cavityList.index(self.cavity.GetValue())
        set_phase_pv = self.cavity_set_phase[index]
        rb_phase_pv = self.cavity_get_phase[index]
        #create_context()
        phase_set_pv = PV(set_phase_pv)
        phase_set_pv.put(round(rfPhase, 1))  
        time.sleep(1)
        phase_get_pv = PV(rb_phase_pv)
        self.rbPhase.SetLabel('{0:.2f}'.format(phase_get_pv.get()))
        #destroy_context()   

    def OnGetBpmPhase(self, event):
        index = self.cavityList.index(self.cavity.GetValue())
        bpm_phases = []
        read_num = 5
        bpm_pv = PV(self.bpm_pv[index])
        for i in range(read_num):
            bpm_phases.append(bpm_pv.get())
            time.sleep(1)
        avg_bpm_phase = np.average(bpm_phases)
        self.bpm_phase_label.SetLabel('BPM {0}'.format(re.findall('\d+', self.bpm_pv[index])[0]))
        self.bpm_phase.SetLabel('{0:.1f}'.format(avg_bpm_phase))

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

    def button_fit(self):
        x = self.scan_line.get_xdata()
        y = self.scan_line.get_ydata()
        self.fit(x, y)

    def fit(self, x, y):
        index = self.cavityList.index(self.cavity.GetValue())
        rf_phase_pv = self.cavity_set_phase[index]
        f = 162.5e6
        #f = 162.4699e6
        bpm_harm = 1
        bpm_polarity = 1
        if self.cavityList[index].startswith("buncher"):
            EpeakFactor = 600
        elif self.cavityList[index].startswith("cm3"):
            if self.cavityList[index].startswith("cm3-5"):
                bpm_harm = 2
                bpm_polarity = -1
            EpeakFactor = 32
        elif self.cavityList[index].startswith("cm4"):
            f = 325e6
            rf_phase_pv = None
            #f = 324.9398e6
            if not (self.cavityList[index].startswith("cm4-6") or self.cavityList[index].startswith("cm4-5")):
                bpm_polarity = -1
            else:
                bpm_harm = 0.5
            EpeakFactor = 25         
        else:
            EpeakFactor = 25
        Win = float(self.injectEnergy.GetValue())
        synch_phases = np.loadtxt('synch-phases/phases.dat')
        distance = self.distance_cav_bpm[index]
        twPhase = synch_phases[index]
        fieldName = self.field_names[index]
        slope = self.slopes[index]
        focus_mode = self.focus.GetSelection()
        start_phase = x[0] 
        type_name = self.type.GetValue()
        mass = self.particle_mass_charge[type_name]['mass']
        charge = self.particle_mass_charge[type_name]['charge']  
        cavity_name = self.cavityList[index]      
        
        t = FitThread(self, Win, distance, twPhase, fieldName, slope, x, y,
                      EpeakFactor, start_phase, focus_mode, f, bpm_harm, bpm_polarity,
                      rf_phase_pv, mass, charge, cavity_name) 
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

    def read_synch_phase(self, filename):
        tracewin_file = open(filename, 'r')
        synch_phase_file = open('synch-phases/phases.dat', 'w')
        tracewin_phases = []
        for line in tracewin_file:
            if not line.strip() or line.startswith(';'):
                continue
            line_split = line.split()
            if line_split[-1].startswith(('hwr', 'cav', 'buncher')):
                tracewin_phases.append(line_split[3])
        tracewin_file.close()
        for phase in tracewin_phases:
            synch_phase_file.write('{0}\n'.format(phase))
        synch_phase_file.close()

        if not tracewin_phases:
            self.statusBar.SetStatusText(('Synchronous phases loaded error'), 0)
        else:
            self.statusBar.SetStatusText(('Synchronous phases loaded successfully'), 0)

    def read_fit(self, filename):
        if filename:
            self.clear_graph()
            f = open(filename, 'r')
            data = f.readlines()
            f.close()
            x = []
            y = []
            for line in data:
                line_data = line.strip().split()
                x.append(float(line_data[0]))
                y.append(float(line_data[1]))
   
            self.fit(x, y)
            #self.updateGraph(self.fit_line, x_plot, y_plot)

    def load(self, event):
        dlg = wx.FileDialog(self, "Open synchronous phase file...", os.getcwd(), style=wx.OPEN, wildcard=self.wildcard)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            self.read_synch_phase(filename)
        dlg.Destroy()

    def open(self, event):
        dlg = wx.FileDialog(self, "Open phase file...", os.getcwd(), style=wx.OPEN, wildcard=self.wildcard)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            self.read_fit(filename)
        dlg.Destroy()

    def get_ready(self):
        self.data_changed = False
        self.startButton.Disable()
        self.acquire_control_values()
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
        self.thread.timeToPause.set()
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

        self.thread = WorkThread(self, self.cavity_id, self.start_phase, self.stop_phase, self.phase_step,
                                 self.delay_before_scan, self.delay_read, self.num_read)
        self.thread.start()

        self.pauseButton.Enable()
        self.stopButton.Enable()
        self.pauseButton.SetLabel('pause')

    def OnStop(self, event):
        self.reset_buttons()
        self.stop_thread()

    def OnPause(self, event):
        if self.pauseButton.GetLabel() == 'pause':
            self.thread.timeToPause.clear()
            self.thread.pause = True
            self.pauseButton.SetLabel('resume')
        else:
            self.thread.timeToPause.set()
            self.pauseButton.SetLabel('pause')

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


