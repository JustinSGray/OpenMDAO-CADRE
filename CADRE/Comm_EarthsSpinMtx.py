from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_EarthsSpinMtx(Component):

    def __init__(self, n):
        super(Comm_EarthsSpinMtx, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('q_E', Array(iotype='in', shape=(4, self.n)))
        self.add('O_IE', Array(iotype='out', shape=(3, 3, self.n)))

    def linearize(self):
        self.J = self.lib.computerotmtxjacobian(self.n, self.q_E)

    def execute(self, val, sol):
        self.O_IE = self.lib.computerotationmtx(self.n, self.q_E)

    def applyDer(self, arg, result):
        if 'q_E' in arg:
            result['O_IE'] = 0*result['O_IE']
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['O_IE'][u,v,:] += self.J[:,u,v,k] * arg['q_E'][k,:]
        return result

    def applyDerT(self, arg, result):
        if 'O_IE' in arg:
            result['q_E'] = 0*result['q_E']
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['q_E'][k,:] += self.J[:,u,v,k] * arg['O_IE'][u,v,:]
        return result
