from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_EarthsSpin(Component):

    def __init__(self, n):
        super(Comm_EarthsSpin, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        
        self.add('q_E', Array(iotype='out', shape=(4, self.n)))
        self.add('t', Array(iotype='out', shape=(self.n)))

    def linearize(self, val):
        self.dq_dt = self.lib.computejacobianqe(self.n, self.t)

    def execute(self, val, sol):
        self.q_E += self.lib.computeqe(self.n, self.t)

    def applyDer(self, arg, result):
        for k in range(4):
            result['q_E'][k,:] += self.dq_dt[:,k] * arg['t'][:]
        return result

    def applyDerT(self, arg, result):
        for k in range(4):
            result['t'][:] += self.dq_dt[:,k] * arg['q_E'][k,:]
        return result
