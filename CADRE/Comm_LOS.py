from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_LOS(Component):

    def __init__(self, n):
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        
        self.add('CommLOS', Array(iotype='out', shape=(self.n)))
        
        self.add('r_e2g_I', Array(iotype='in', shape=(3, self.n)))
        self.add('r_e2g_E', Array(iotype='in', shape=(3, self.n)))

    def linearize(self, val):
        result = self.lib.computejacobianlos(self.n, self.r_b2g_I[:], 
                                             self.r_e2g_I[:])
        self.dLOS_drb, self.dLOS_dre = result

    def execute(self, val, sol):
        self.CommLOS[:] = self.lib.computelos(self.n, self.r_b2g_I[:], 
                                              self.r_e2g_I[:])

    def applyDer(self, arg, result):
        for k in xrange(3):
            result['CommLOS'][:] = self.dLOS_drb[:,k] * arg['r_b2g_I'][k,:]
            result['CommLOS'][:] += self.dLOS_dre[:,k] * arg['r_e2g_I'][k,:]
        return result
    
    def applyDerT(self, arg, result):
        for k in xrange(3):
            result['r_b2g_I'][k,:] = self.dLOS_drb[:,k] * arg['CommLOS'][:]
            result['r_e2g_I'][k,:] += self.dLOS_dre[:,k] * arg['CommLOS'][:]
        return result
