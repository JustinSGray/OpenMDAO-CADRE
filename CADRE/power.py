import numpy as np
import scipy.sparse

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import MBI

class Power_CellVoltage( Component ):

    def __init__(self, n=2): 
        super(Power_CellVoltage, self).__init__()

        self.n = n 

        self.add('LOS', Array(np.zeros((n, ), order='F'), size=(n, ), dtype=np.float, iotype="in", 
            desc="Line of Sight over Time"))
        self.add('temperature', Array(np.zeros((5,n, ), order='F'), size=(5,n,), dtype=np.float,
                                      iotype="in"))
        self.add('exposedArea', Array(np.zeros((7,12,n, ), order='F'), size=(7,12,n,), dtype=np.float,
                                      iotype="in"))
        self.add('Isetpt', Array(np.zeros((12,n, ), order='F'), size=(12,n,), dtype=np.float,
                                      iotype="in"))

        self.add('V_sol', Array(np.zeros((12,n, ), order='F'), size=(12,n,), dtype=np.float,
                                      iotype="out"))

        dat = np.genfromtxt('CADRE/data/Power/curve.dat')
        nT, nA, nI = dat[:3]
        T = dat[3:3+nT]
        A = dat[3+nT:3+nT+nA]
        I = dat[3+nT+nA:3+nT+nA+nI]
        V = dat[3+nT+nA+nI:].reshape((nT,nA,nI),order='F')

        self.MBI = MBI.MBI(V,[T,A,I],[6,6,15],[3,3,3])

        self.x = np.zeros((84*self.n,3),order='F')
        self.xV = self.x[:].reshape((self.n,7,12,3),order='F')
        self.dV_dL = np.zeros((self.n,12), order='F')
        self.dV_dT = np.zeros((self.n,12,5), order='F')
        self.dV_dA = np.zeros((self.n,7,12), order='F')
        self.dV_dI = np.zeros((self.n,12), order='F')

    def setx(self):
        for c in range(7):
            for p in range(12):
                i = 4 if p < 4 else (p % 4)
                self.xV[:,c,p,0] = self.temperature[i,:]
                self.xV[:,c,p,1] = self.LOS[:] * self.exposedArea[c,p,:]
                self.xV[:,c,p,2] = self.Isetpt[p,:]
        
    def execute(self): 

        self.setx()
        self.raw = self.MBI.evaluate(self.x)[:,0].reshape((self.n,7,12),order='F')
        for c in range(7):
            for p in range(12):
                self.V_sol[p,:] += self.raw[:,c,p]

    def linearize(self): 

        self.setx()
        self.raw1 = self.MBI.evaluate(self.x,1)[:,0].reshape((self.n,7,12),order='F')
        self.raw2 = self.MBI.evaluate(self.x,2)[:,0].reshape((self.n,7,12),order='F')
        self.raw3 = self.MBI.evaluate(self.x,3)[:,0].reshape((self.n,7,12),order='F')
        self.dV_dL[:] = 0.0
        self.dV_dT[:] = 0.0
        self.dV_dA[:] = 0.0
        self.dV_dI[:] = 0.0
        for c in range(7):
            for p in range(12):
                i = 4 if p < 4 else (p % 4)
                self.dV_dL[:,p] += self.raw2[:,c,p] * self.exposedArea[c,p,:]
                self.dV_dT[:,p,i] += self.raw1[:,c,p]
                self.dV_dA[:,c,p] += self.raw2[:,c,p] * self.LOS[:]
                self.dV_dI[:,p] += self.raw3[:,c,p]

    def applyDer(self, arg, result):

        if not result['V_sol'] :
            result['V_sol'] = np.zeros( (12, self.n,) )

        if 'LOS' in arg: 
            for p in range(12):
                result['V_sol'][p,:] += self.dV_dL[:,p] * arg['LOS'][:]

        if 'temperature' in arg: 
            for p in range(12):
                i = 4 if p < 4 else (p % 4)
                result['V_sol'][p,:] += self.dV_dT[:,p,i] * arg['temperature'][i,:]

        if 'Isetpt' in arg: 
            for p in range(12):
                result['V_sol'][p,:] += self.dV_dI[:,p] * arg['Isetpt'][p,:]

        if 'exposedArea' in arg: 
            for p in range(12):
                for c in range(7):
                    result['V_sol'][p,:] += self.dV_dA[:,c,p] * arg['exposedArea'][c,p,:]

        return result   

    def applyDerT(self, arg, result): 

        if not result['LOS'] :
            result['LOS'] = np.zeros( (self.n,) )

        if not result['temperature'] :
            result['temperature'] = np.zeros( (5,self.n,) )

        if not result['exposedArea'] :
            result['exposedArea'] = np.zeros( (7,12,self.n,) )

        if not result['Isetpt'] :
            result['Isetpt'] = np.zeros( (12, self.n,) )

        if 'V_sol' in arg: 
            for p in range(12):
                i = 4 if p < 4 else (p % 4)
                result['LOS'][:] += self.dV_dL[:,p] * arg['V_sol'][p,:]
                result['temperature'][i,:] += self.dV_dT[:,p,i] * arg['V_sol'][p,:]
                result['Isetpt'][p,:] += self.dV_dI[:,p] * arg['V_sol'][p,:]
                for c in range(7):
                    result['exposedArea'][c,p,:] += self.dV_dA[:,c,p] * arg['V_sol'][p,:]

        return result

class Power_SolarPower( Component ):

    def __init__(self, n=2): 
        super(Power_SolarPower, self).__init__()

        self.n = n 

        self.add('Isetpt', Array(np.zeros((12,n, ), order='F'), size=(12,n,), dtype=np.float,
                                      iotype="in"))
        self.add('V_sol', Array(np.zeros((12,n, ), order='F'), size=(12,n,), dtype=np.float,
                                      iotype="in"))

        self.add('P_sol', Array(np.zeros((n, ), order='F'), size=(n,), dtype=np.float,
                                      iotype="out"))

    def execute(self): 

        for p in range(12):
            self.P_sol[:] += self.V_sol[p,:] * self.Isetpt[p,:]

    def applyDer(self, arg, result):

        if not result['P_sol'] :
            result['P_sol'] = np.zeros( (self.n,) )

        if 'V_sol' in arg: 
            for p in range(12):
                result['P_sol'][:] += arg['V_sol'][p,:] * self.Isetpt[p,:]

        if 'Isetpt' in arg: 
            for p in range(12):
                result['P_sol'][:] += arg['Isetpt'][p,:] * self.V_sol[p,:]

        return result   

    def applyDerT(self, arg, result): 

        if not result['V_sol'] :
            result['V_sol'] = np.zeros( (12, self.n,) )

        if not result['Isetpt'] :
            result['Isetpt'] = np.zeros( (12, self.n,) )

        if 'P_sol' in arg:
            for p in range(12):
                result['V_sol'][p,:] += arg['P_sol'][:] * self.Isetpt[p,:]
                result['Isetpt'][p,:] += self.V_sol[p,:] * arg['P_sol'][:]

        return result

class Power_Total( Component ):

    def __init__(self, n=2): 
        super(Power_Total, self).__init__()

        self.n = n 

        self.add('P_sol', Array(np.zeros((n, ), order='F'), size=(n,), dtype=np.float,
                                      iotype="in"))
        self.add('P_comm', Array(np.zeros((n, ), order='F'), size=(n,), dtype=np.float,
                                      iotype="in"))
        self.add('P_RW', Array(np.zeros((3,n, ), order='F'), size=(3,n,), dtype=np.float,
                                      iotype="in"))

        self.add('P_bat', Array(np.zeros((n, ), order='F'), size=(n,), dtype=np.float,
                                      iotype="out"))

    def execute(self): 
        self.P_bat[:] += self.P_sol[:] - 5*self.P_comm[:] - np.sum(self.P_RW[:], 0) - 2.0

    def applyDer(self, arg, result):

        # for k in range(3):
        #     res('P_bat')[:] -= arg('P_RW')[k,:]


        if not result['P_bat'] :
            result['P_bat'] = np.zeros( (self.n,) )

        if 'P_sol' in arg: 
            result['P_bat'][:] += arg['P_sol'][:]

        if 'P_comm' in arg: 
            result['P_bat'][:] -= 5 * arg['P_comm'][:]

        if 'P_RW' in arg: 
            for k in range(3):
                result['P_bat'][:] -= arg['P_RW'][k,:]

        return result   

    def applyDerT(self, arg, result): 

        if not result['P_sol'] :
            result['P_sol'] = np.zeros( (self.n,) )
        if not result['P_comm'] :
            result['P_comm'] = np.zeros( (self.n,) )
        if not result['P_RW'] :
            result['P_RW'] = np.zeros( (3,self.n,) )

        if 'P_bat' in arg:
            result['P_sol'][:] += arg['P_bat'][:]
            result['P_comm'][:] -= 5*arg['P_bat'][:]
            for k in range(3):
                result['P_RW'][k,:] -= arg['P_bat'][:]

        return result

