from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_VectorAnt(Component):

    def __init__(self, n):
        super(Comm_VectorAnt, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib

        self.add('r_e2g_A', Array(iotype='out', shape=(3, self.n)))
        
        self.add('r_b2g_B', Array(iotype='in', shape=(3, self.n)))
        self.add('O_AB', Array(iotype='in', shape=(3, 3, self.n)))

    def linearize(self):
        result = self.lib.computepositionrotdjacobian(self.n, self.r_b2g_B[:], 
                                                      self.O_AB[:])
        self.J1, self.J2 = result

    def execute(self):
        self.r_b2g_A[:] = self.lib.computepositionrotd(self.n, self.r_b2g_B[:], 
                                                        self.O_AB[:])

    def applyDer(self, arg, result):
        if 'O_AB' in arg and 'r_b2g_B' in arg:
            result['r_b2g_A'] = np.zeros((3, self.n))
            for k in range(3):
                for u in range(3):
                    for v in range(3):
                        result['r_b2g_A'][k,:] += self.J1[:,k,u,v] * arg['O_AB'][u,v,:]
                for j in range(3):
                    result['r_b2g_A'][k,:] += self.J2[:,k,j] * arg['r_b2g_B'][j,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_e2g_A' in arg:
            result['r_b2g_B'] = np.zeros((3, self.n))
            result['O_AB'] = np.zeros((3, 3, self.n))
            for k in range(3):
                for u in range(3):
                    for v in range(3):
                        result['O_AB'][u,v,:] += self.J1[:,k,u,v] * arg['r_b2g_A'][k,:]
                for j in range(3):
                    result['r_b2g_B'][j,:] += self.J2[:,k,j] * arg['r_b2g_A'][k,:]
        return result
