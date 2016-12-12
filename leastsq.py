#!/usr/bin/env python
import scipy.constants as C
from scipy import interpolate
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from pyswarm import pso
import sys

#fieldName = 'Exyz.txt'
#scanPhaseFile = 'HWR.txt'
#injectEnergy = 2.15
#distance = 0.1

w = 162.5 * 10 ** 6


#distance = 0.446
#distance = 0.1026
mass = 938.272

def calTraceWinPhase(Win, c, phase_in, distance, l, dz, Ez):
    W = Win
    t = 0
    a = 0
    b = 0
    for i in range(l - 1):
        gamma = W / mass + 1
        beta = (1 - gamma**-2)**0.5
        phi = phase_in + 2 * C.pi * w * (t + 0.5 * dz / ( beta * C.c))
        W = W +  c * 0.5 * (Ez[i] + Ez[i+1]) * np.cos(phi) * dz
        gammaExit = W / mass + 1
        betaExit = (1 - gammaExit ** -2) ** 0.5
        a += c * 0.5 * (Ez[i] + Ez[i+1]) * np.sin(phi) * dz
        b += c * 0.5 * (Ez[i] + Ez[i+1]) * np.cos(phi) * dz
        t += dz / (0.5 * (beta + betaExit) * C.c)
    t += distance / (betaExit * C.c)
    traceWin_phi = np.arctan(a / b)
    if (b < 0 and a < 0):
        traceWin_phi -= C.pi
    elif (b < 0 and a > 0):
        traceWin_phi += C.pi

    return traceWin_phi, a, b 

def getEntrPhase(twPhase, Win, f, distance, l, dz, Ez):
    minimum = 0xffffffff
    best_fit_phase = -C.pi

    for phi in np.arange(-C.pi, C.pi, C.pi / 180):
        err = abs(calTraceWinPhase(Win, f, phi, distance, l, dz, Ez)[0] - twPhase) 
        if (err < minimum):
            minimum = err
            best_fit_phase = phi

    return best_fit_phase

def energyGain(Win, c, phase_in, distance, l, dz, Ez):
    W = Win
    t = 0
    a = 0
    b = 0
    for i in range(l - 1):
        phi = phase_in + 2 * C.pi * w * t
        gamma = W / mass + 1
        beta = (1 - gamma**-2)**0.5
        W = W +  c * 0.5 * (Ez[i] + Ez[i+1]) * np.cos(phi) * dz
        gammaExit = W / mass + 1
        betaExit = (1 - gammaExit ** -2) ** 0.5
        #a += c * Ez[i] * sin(phi) * dz
        #b += c * Ez[i] * cos(phi) * dz
        t += dz / (0.5 * (beta + betaExit) * C.c)
    t += distance / (betaExit * C.c)
    #traceWin_phi = actan(a / b)
    return -(w * t) * 180 * 2

def residuals(p, y, injectEnergy, distance, l, dz, Ez, x):
    c, phase_in, offset = p
    err =  np.array(y) - [energyGain(injectEnergy, c, phase_in + e, distance, l, dz, Ez) + offset for e in x]
    return err

def getTWPhase(cav_phases, bpm_phases, injectEnergy, distance, twissWinPhase, fieldName, start_phase, slope, EpeakFactor, focus_mode):
    data = np.loadtxt(fieldName)
    z = np.linspace(data[0, 0], data[-1, 0], 3000)
    f = interpolate.interp1d(data[:, 0], data[:, 3], kind='slinear')
    Ez = f(z)
    l = len(z)
    dz = (data[-1, 0] - data[0, 0]) / 3000
    #fitStep = step * slope
    #fitPointNum = len(cav_phases)

    #x = -np.arange(fitPointNum) * fitStep
    if focus_mode == 0:
        x = -np.asarray(cav_phases) * np.pi / 180 * slope
    else:
        x = np.asarray(cav_phases) * np.pi / 180 * slope
    p0 = [1, 0, 0]
    plsq = leastsq(residuals, p0, args=(bpm_phases, injectEnergy, distance, l, dz, Ez, x))
    error = np.std(residuals(plsq[0], bpm_phases, injectEnergy, distance, l, dz, Ez, x))

    twissWinPhase = twissWinPhase * C.pi / 180
    scaleFactor = plsq[0][0]
    xopt = getEntrPhase(twissWinPhase, injectEnergy, scaleFactor, distance, l, dz, Ez)

    y = [energyGain(injectEnergy, plsq[0][0], plsq[0][1] + e, distance, l, dz, Ez) for e in x]

    if focus_mode == 0:
        rfPhase = (plsq[0][1] + x[0] - xopt) * 180 / C.pi / slope + start_phase
    else:
        rfPhase = -(plsq[0][1] + x[0] - xopt) * 180 / C.pi / slope + start_phase

    exit_energy = calTraceWinPhase(injectEnergy, plsq[0][0], xopt, distance, l, dz, Ez)[2]
    rfPhase = phaseWrappingFunction(rfPhase, slope)
    return rfPhase, exit_energy, plsq[0][0] * EpeakFactor, error, cav_phases, y + plsq[0][2]

def phaseWrappingFunction(inValue, slope):
    outValue = inValue
    if (abs(inValue) > 180 / slope):
        outValue += 180 / slope
        while (outValue < 0):
            outValue += 360 / slope
        outValue = outValue % (360 / slope)
        outValue -= 180 / slope
    return outValue
