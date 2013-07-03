from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

from MBI.MBI import MBI

import numpy as np


class BsplineParameters(Component):

    def __init__(self, n):
        super(BsplineParameters, self).__init__()
        self.n = n
        self.add('m', 5)
        m = self.m
        self.add('t1', 0)
        self.add('t2', 48000)        
        self.B = MBI(np.zeros(n), [np.linspace(self.t1,self.t2,n)], [m], [4]).getJacobian(0,0)
        self.Bdot = MBI(np.zeros(n), [np.linspace(self.t1,self.t2,n)], [m], [4]).getJacobian(1,0)
        self.BT = self.B.transpose()
        self.BdotT = self.Bdot.transpose()        
        
        self.add('CP_P_comm', Array(np.zeros((m,)), size=(m,), dtype=float, iotype='in'))
        self.add('CP_gamma', Array(np.zeros((m,)), size=(m,), dtype=float, iotype='in'))
        self.add('CP_Isetpt', Array(np.zeros((12,n)), size=(12,n), dtype=float, iotype='in'))

        self.add('P_comm', Array(np.ones((n,)), size=(n,), dtype=float, iotype='out'))
        self.add('Gamma', Array(0.1*np.ones((n,)), size=(n,), dtype=float, iotype='out'))
        self.add('Isetpt',Array(0.2*np.ones((12,n)), size=(12,n), dtype=float, iotype='out'))

    def execute(self):
        self.P_comm[:] = self.B.dot(self.CP_P_comm[:])
        self.Gamma[:] = self.B.dot(self.CP_gamma[:])
        for k in range(12):
            self.Isetpt[k,:] = self.B.dot(self.CP_Isetpt[k,:])

    def applyDer(self, arg, result):
        result['P_comm'] = np.zeros((self.n,))
        result['Gamma'] = np.zeros((self.n,))
        result['Isetpt'] = np.zeros((12,self.n))
        
        if 'CP_P_comm' in arg:
            result['P_comm'][:] += self.B.dot(arg['CP_P_comm'][:])
        if 'CP_gamma' in arg:
            result['Gamma'][:] += self.B.dot(arg['CP_gamma'][:])
        if 'CP_Isetpt' in arg:
            for k in range(12):
                result['Isetpt'][k,:] += self.B.dot(arg['CP_Isetpt'][k,:])
        return result

    def applyDerT(self, arg, result):
        result['CP_P_comm'] = np.zeros((self.m,))
        result['CP_gamma'] = np.zeros((self.m,))
        result['CP_Isetpt'] = np.zeros((12,self.m))
        
        if 'P_comm' in arg:
            result['CP_P_comm'][:] += self.BT.dot(arg['P_comm'][:])
        if 'Gamma' in arg:
            result['CP_gamma'][:] += self.BT.dot(arg['Gamma'][:])
        if 'Isetpt' in arg:
            for k in range(12):
                result['CP_Isetpt'][k,:] += self.BT.dot(arg['Isetpt'][k,:])
        return result
