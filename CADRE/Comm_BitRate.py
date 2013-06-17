from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_BitRate(Component):

    def __init__(self, n):
        super(Comm_BitRate, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib

        self.add('Dr', Array(iotype='out', shape=(self.n)))

        self.add('P_comm', Array(iotype='in', shape=(self.n)))
        self.add('gain', Array(iotype='in', shape=(self.n)))
        self.add('GSdist', Array(iotype='in', shape=(self.n)))
        self.add('CommLOS', Array(iotype='in', shape=(self.n)))

    def linearize(self):
        response = self.lib.computejacobiandr(self.n, 
                                              self.P_comm, 
                                              self.gain, 
                                              self.GSdist, 
                                              self.CommLOS)
        self.dD_dP, self.dD_dGt, self.dD_dS, self.dD_dLOS = response
        
    def execute(self):
        self.Dr = self.lib.computedr(self.n, 
                                           self.P_comm, 
                                           self.gain, 
                                           self.GSdist, 
                                           self.CommLOS)

    def applyDer(self, arg, result):
        if 'P_comm' in arg:
            result['Dr'][:] = self.dD_dP * arg['P_comm'][:]
        if 'gain' in arg:
            result['Dr'][:] += self.dD_dGt * arg['gain'][:]
        if 'GSdist' in arg:
            result['Dr'][:] += self.dD_dS * arg['GSdist'][:]
        if 'CommLOS' in arg:
            result['Dr'][:] += self.dD_dLOS * arg['CommLOS'][:]
        return result

    def applyDerT(self, arg, result):
        if 'Dr' in arg:
            result['P_comm'][:] = self.dD_dP * arg['Dr'][:]
            result['gain'][:] = self.dD_dGt * arg['Dr'][:]
            result['GSdist'][:] = self.dD_dS * arg['Dr'][:]
            result['CommLOS'][:] = self.dD_dLOS * arg['Dr'][:]
        return result
