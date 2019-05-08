#!/usr/bin/env python
import time
import sys
import numpy as np

from pcaspy import Driver, SimpleServer

#prefix = 'SCRF:'
pvdb = { 
    'SCRF:CAV1:PHASE:SETPOINT': {
    },
    'Bpm:6-P11': {
    },
    'LLRF:Buncher1:Phase_unlock': {},
    'LLRF:Buncher2:Phase_unlock': {},
    'SCRF:CAV1:READY':{},
    'STATUS': {
        'type': 'enum',
        'enums': ['DONE', 'BUSY']
    },
}

class myDriver(Driver):
    def __init__(self, fname):
        Driver.__init__(self)
        self.tid = None 
        self.data = np.loadtxt(fname)

    def write(self, reason, value):
        status = True
        # take proper actions
        if reason == 'SCRF:CAV1:PHASE:SETPOINT':
            idx = abs(self.data[:, 0] - value) < 1
            bpm_phase = self.data[idx, 1]
            self.setParam(reason, value)
            self.setParam('Bpm:6-P11', bpm_phase)
        else:
            self.setParam(reason, value)
        self.updatePVs()

        return status

if __name__ == '__main__':
    server = SimpleServer()
    server.createPV('', pvdb)
    driver = myDriver('cm1-1-bak.txt')

    while True:
        # process CA transactions
        server.process(0.1)
