#!/usr/bin/env python
import time
import sys
import threading
import subprocess
import shlex
import numpy as np

from pcaspy import Driver, SimpleServer

prefix = ''
pvdb = { 
    'LLRF:Buncher1:Loop_State': {},
    'LLRF:Buncher2:Loop_State': {},
    'LLRF:Buncher1:CAVITY_PHASE': {},
    'SCRF:CAV1:READY': {'value': 1},
    'SCRF:CAV2:READY': {'value': 1},
    'SCRF:CAV3:READY': {'value': 1},
    'SCRF:CAV4:READY': {'value': 1},
    'SCRF:CAV5:READY': {'value': 1},
    'SCRF:CAV6:READY': {'value': 1},
    'SCRF:CAV7:READY': {'value': 1},
    'SCRF:CAV8:READY': {'value': 1},
    'SCRF:CAV9:READY': {'value': 1},
    'SCRF:CAV10:READY': {'value': 1},
    'SCRF:CAV11:READY': {'value': 1},
    'SCRF:CAV12:READY':{'value': 1},
    'SCRF:CAV13:READY':{'value': 1},
    'SCRF:CAV14:READY':{'value': 1},
    'SCRF:CAV15:READY':{'value': 1},
    'SCRF:CAV16:READY':{'value': 1},
    'SCRF:CAV17:READY':{'value': 1},
    'SCRF:CAV18:READY':{'value': 1},
    'SCRF:CAV19:READY':{'value': 1},
    'SCRF:CAV20:READY':{'value': 1},
    'SCRF:CAV21:READY':{'value': 1},
    'SCRF:CAV22:READY':{'value': 1},
    'SCRF:CAV23:READY':{'value': 1},
    'LLRF:Buncher1:PHA_SET': {},
    'LLRF:Buncher2:PHA_SET': {},
    'ADS:SCRF:CM01:CAV01:PHASE:WRITE': {},
    'ADS:SCRF:CM01:CAV02:PHASE:WRITE': {},
    'ADS:SCRF:CM01:CAV03:PHASE:WRITE': {},
    'ADS:SCRF:CM01:CAV04:PHASE:WRITE': {},
    'ADS:SCRF:CM01:CAV05:PHASE:WRITE': {},
    'ADS:SCRF:CM01:CAV06:PHASE:WRITE': {},
    'ADS:SCRF:CM02:CAV01:PHASE:WRITE': {},
    'ADS:SCRF:CM02:CAV02:PHASE:WRITE': {},
    'ADS:SCRF:CM02:CAV03:PHASE:WRITE': {},
    'ADS:SCRF:CM02:CAV04:PHASE:WRITE': {},
    'ADS:SCRF:CM02:CAV05:PHASE:WRITE': {},
    'ADS:SCRF:CM02:CAV06:PHASE:WRITE': {},
    'ADS:SCRF:CM03:CAV01:PHASE:WRITE': {},
    'ADS:SCRF:CM03:CAV02:PHASE:WRITE': {},
    'ADS:SCRF:CM03:CAV03:PHASE:WRITE': {},
    'ADS:SCRF:CM03:CAV04:PHASE:WRITE': {},
    'ADS:SCRF:CM03:CAV05:PHASE:WRITE': {},
    'ADS:SCRF:CM03:CAV06:PHASE:WRITE': {},
    'ADS:SCRF:CM04:CAV01:PHASE:WRITE': {},
    'ADS:SCRF:CM04:CAV02:PHASE:WRITE': {},
    'ADS:SCRF:CM04:CAV03:PHASE:WRITE': {},
    'ADS:SCRF:CM04:CAV04:PHASE:WRITE': {},
    'ADS:SCRF:CM04:CAV05:PHASE:WRITE': {},
    'LLRF:Buncher1:VACC': {},
    'LLRF:Buncher2:VACC': {},    
    'SCRF:CAV1:EPIC': {},
    'SCRF:CAV2:EPIC': {},
    'SCRF:CAV3:EPIC': {},
    'SCRF:CAV4:EPIC': {},
    'SCRF:CAV5:EPIC': {},
    'SCRF:CAV6:EPIC': {},
    'SCRF:CAV7:EPIC': {},
    'SCRF:CAV8:EPIC': {},
    'SCRF:CAV9:EPIC': {},
    'SCRF:CAV10:EPIC': {},
    'SCRF:CAV11:EPIC': {},
    'SCRF:CAV12:EPIC': {},
    'SCRF:CAV13:EPIC': {},
    'SCRF:CAV14:EPIC': {},
    'SCRF:CAV15:EPIC': {},
    'SCRF:CAV16:EPIC': {},
    'SCRF:CAV17:EPIC': {},
    'SCRF:CAV18:EPIC': {},
    'SCRF:CAV19:EPIC': {},
    'SCRF:CAV20:EPIC': {},
    'SCRF:CAV21:EPIC': {},
    'SCRF:CAV22:EPIC': {},
    'SCRF:CAV23:EPIC': {},
    'Bpm:2-P11': {},
    'Bpm:3-P11': {},
    'Bpm:5-P11': {},
    'Bpm:6-P11': {},
    'Bpm:7-P11': {},
    'Bpm:8-P11': {},
    'Bpm:9-P11': {},
    'Bpm:10-P11': {},
    'Bpm:11-P11': {},
    'Bpm:12-P11': {},
    'Bpm:13-P11': {},
    'Bpm:14-P11': {},
    'Bpm:15-P11': {},
    'Bpm:16-P11': {},
    'Bpm:17-P11': {},
    'Bpm:18-P11': {},
    'Bpm:19-P11': {},
    'Bpm:20-P11': {},
    'Bpm:21-P11': {},
    'Bpm:22-P11': {},
    'Bpm:23-P11': {},
    'Bpm:24-P11': {},
    'Bpm:25-P11': {},
}

class myDriver(Driver):
    def __init__(self):
        Driver.__init__(self)
        self.data = np.loadtxt('4.28/buncher1.txt', skiprows=1)

    def write(self, reason, value):
        status = True
        # take proper actions
        if reason.startswith('LLRF'):
            idx = np.where(abs(self.data[:, 0]-value)<0.1)
            print(self.data[idx, 1])
            self.setParam('Bpm:2-P11', self.data[idx, 1])
            self.setParam('Bpm:3-P11', self.data[idx, 2])
            self.setParam('LLRF:Buncher1:CAVITY_PHASE', self.data[idx, 0])
        self.setParam(reason, value)
        self.updatePVs()
        return status

if __name__ == '__main__':
    server = SimpleServer()
    server.createPV(prefix, pvdb)

    driver = myDriver()

    while True:
        # process CA transactions
        server.process(0.1)
