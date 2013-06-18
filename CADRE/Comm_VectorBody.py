from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_VectorBody(Component):

    def __init__(self, n):
        super(Comm_VectorBody, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('r_b2g_B', Array(iotype='out', shape=(3, self.n)))
        
        self.add('r_b2g_I', Array(iotype='in', shape=(3, self.n)))
        self.add('O_BI', Array(iotype='in', shape=(3, 3, self.n)))

    def linearize(self, val):
        result = self.lib.computepositionrotdjacobian(self.n, self.r_b2g_I[:], 
                                                      self.O_BI[:])
        self.J1, self.J2 = result

    def execute(self, val, sol):
        self.r_b2g_B[:] = self.lib.computepositionrotd(self.n, self.r_b2g_I[:], 
                                                       elf.O_BI[:])

    def applyDer(self, arg, result):
        if 'O_BI' in arg and 'r_b2g_I' in arg:
            result['r_b2g_B'] = np.zeros((3, self.n))
            for k in range(3):
                for u in range(3):
                    for v in range(3):
                        result['r_b2g_B'][k,:] += self.J1[:,k,u,v] * arg('O_BI')[u,v,:]
                for j in range(3):
                    result['r_b2g_B'][k,:] += self.J2[:,k,j] * arg('r_b2g_I')[j,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_b2g_B' in arg:
            result['O_BI'] = np.zeros((3, 3, self.n))
            result['r_b2g_I'] = np.zeros((3, self.n))
            for k in range(3):
                for u in range(3):
                    for v in range(3):
                        result['O_BI'][u,v,:] += self.J1[:,k,u,v] * arg['r_b2g_B'][k,:]
                for j in range(3):
                    result['r_b2g_I'][j,:] += self.J2[:,k,j] * arg['r_b2g_B'][k,:]
        return result
