from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_GSposECI(Component):

    def __init__(self, n):
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib

        self.add('r_e2g_I', Array(iotype='out', shape=(3, self.n)))
        
        self.add('O_IE', Array(iotype='in', shape=(3, 3, self.n)))
        self.add('r_e2g_E', Array(iotype='in', shape=(3, self.n)))

    def linearize(self):
        result = self.lib.computepositionrotdjacobian(self.n, 
                                                      self.r_e2g_E[:], 
                                                      self.O_IE[:])
        self.J1, self.J2 = result

    def execute(self):
        self.r_e2g_I[:] = self.lib.computepositionrotd(self.n, 
                                                      self.r_e2g_E[:], 
                                                       self.O_IE[:])

    def applyDer(self, arg, result):
        if 'O_IE' in arg and 'r_e2g_E' in arg:
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['r_e2g_I'][k,:] += self.J1[:,k,u,v] * arg['O_IE'][u,v,:]
                for j in xrange(3):
                    result['r_e2g_I'][k,:] += self.J2[:,k,j] * arg['r_e2g_E'][j,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_e2g_I' in arg:
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['O_IE'][u,v,:] += self.J1[:,k,u,v] * arg['r_e2g_I'][k,:]
                for j in xrange(3):
                    result['r_e2g_E'][j,:] += self.J2[:,k,j] * arg['r_e2g_I'][k,:]
        return result
