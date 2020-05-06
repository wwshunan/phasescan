#!/usr/bin/env python
import scipy.constants as C
from scipy.constants import c
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import sys

w = 162.5 * 10 ** 6

class ScanDataFit(object):
    def __init__(self, payload):
        self.cav_phases_degree = payload['x']
        self.bpm_phases = payload['y']
        self.Win = payload['Win']
        self.distance = payload['distance']
        self.sync_phase = payload['sync_phase']
        self.field_name = payload['field_name']
        self.start_phase = payload['start_phase']
        self.epk_factor = payload['epk_factor']
        self.focus_mode = payload['focus_mode']
        self.freq = payload['freq']
        self.bpm_harm = payload['bpm_harm']
        self.bpm_polarity = payload['bpm_polarity']
        self.m = payload['m']
        self.q = payload['q']

        data = np.loadtxt(self.field_name)
        z = np.linspace(data[0, 0], data[-1, 0], 1001)
        f = interpolate.interp1d(data[:, 0], data[:, 3], kind='slinear')
        self.Ez = f(z)
        self.dz = (data[-1, 0] - data[0, 0]) / 1000

        if self.focus_mode == 0:
            self.cav_phases_rad = -np.asarray(self.cav_phases_degree) * np.pi / 180
        else:
            self.cav_phases_rad = np.asarray(self.cav_phases_degree) * np.pi / 180

    def single_bpm_fit(self):
        p0 = [1, 0, 0]
        plsq = leastsq(self.residuals, p0)
        error = np.std(self.residuals(plsq[0]))

        sync_phase = self.sync_phase * np.pi / 180
        field_factor = plsq[0][0]
        phase_in = plsq[0][1]

        phase_opt = self.get_entr_phase(sync_phase, field_factor)

        computed_bpm_phases = [self.get_bpm_phases(field_factor, phase_in + delta)
                               for delta in self.cav_phases_rad]

        delta_start = self.cav_phases_rad[0]
        if self.focus_mode == 0:
            rf_phase = (phase_in + delta_start - phase_opt) * 180 / np.pi + self.start_phase
        else:
            rf_phase = -(phase_in + delta_start - phase_opt) * 180 / np.pi  + self.start_phase

        print(phase_in, delta_start, phase_opt)
        exit_energy = self.get_process_params(field_factor, phase_opt)[2]
        rf_phase = phase_wrapping(rf_phase)
        return (rf_phase, exit_energy, field_factor * self.epk_factor, error,
                self.cav_phases_degree, computed_bpm_phases + plsq[0][2])

    def get_bpm_phases(self, field_factor, phase_in):
        W = self.Win
        t = 0
        a = 0
        b = 0
        for i in range(len(self.Ez) - 1):
            phi = phase_in + 2 * np.pi * self.freq * t
            gamma = W / self.m + 1
            beta = (1 - gamma**-2)**0.5
            W = W +  self.q * field_factor * 0.5 * (self.Ez[i] + self.Ez[i+1]) * np.cos(phi) * self.dz
            gammaExit = W / self.m + 1
            betaExit = (1 - gammaExit ** -2) ** 0.5
            #a += c * Ez[i] * sin(phi) * dz
            #b += c * Ez[i] * cos(phi) * dz
            t += self.dz / (0.5 * (beta + betaExit) * c)
        t += self.distance / (betaExit * c)
        #traceWin_phi = actan(a / b)
        return -(self.freq * t) * 180 * 2 * self.bpm_harm * self.bpm_polarity

    def residuals(self, p):
        field_factor, phase_in, offset = p
        err =  np.array(self.bpm_phases) - [
            self.get_bpm_phases(field_factor, phase_in + delta) + offset
            for delta in self.cav_phases_rad]
        return err

    def get_entr_phase(self, sync_phase, field_factor):
        minimum = 0xffffffff
        best_fit_phase = -np.pi

        for phi in np.arange(-np.pi, np.pi, np.pi / 180):
            err = abs(self.get_process_params(field_factor, phi)[0] - sync_phase)
            if (err < minimum):
                minimum = err
                best_fit_phase = phi

        return best_fit_phase

    def get_process_params(self, field_factor, phase_in):
        W = self.Win
        t = 0
        a = 0
        b = 0
        for i in range(len(self.Ez) - 1):
            gamma = W / self.m + 1
            beta = (1 - gamma**-2)**0.5
            phi = phase_in + 2 * np.pi * self.freq * (t + 0.5 * self.dz / ( beta * c))
            W = W +  self.q * field_factor * 0.5 * (self.Ez[i] + self.Ez[i+1]) * np.cos(phi) * self.dz
            gammaExit = W / self.m + 1
            betaExit = (1 - gammaExit ** -2) ** 0.5
            a += self.q * field_factor * 0.5 * (self.Ez[i] + self.Ez[i+1]) * np.sin(phi) * self.dz
            b += self.q * field_factor * 0.5 * (self.Ez[i] + self.Ez[i+1]) * np.cos(phi) * self.dz
            t += self.dz / (0.5 * (beta + betaExit) * c)
        t += self.distance / (betaExit * c)
        entr_phase = np.arctan(a / b)
        if (b < 0 and a < 0):
            entr_phase -= np.pi
        elif (b < 0 and a > 0):
            entr_phase += np.pi

        return entr_phase, a, b

def phase_wrapping(inValue):
    outValue = inValue
    if (abs(inValue) > 180):
        outValue += 180
        while (outValue < 0):
            outValue += 360
        outValue = outValue % 360
        outValue -= 180
    return outValue
