from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np
import MBI


class Comm_GainPattern(Component):

    def __init__(self, n, rawG):
        super(Comm_GainPattern, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('gain', Array(iotype='out', shape=(self.n)))
        
        self.add('azimuthGS', Array(iotype='in', shape=(self.n)))
        self.add('elevationGS', Array(iotype='in', shape=(self.n)))

        pi = np.pi
        az = np.linspace(0,2*pi,361)
        el = np.linspace(0,2*pi,361)

        self.MBI = MBI.MBI(rawG,[az,el],[15,15],[4,4])
        self.x = numpy.zeros((self.n,2),order='F')

    def setx(self):
        self.x[:,0] = self.azimuthGS
        self.x[:,1] = self.elevationGS

    def linearize(self):
        result = self.lib.fixangles(self.n, self.azimuthGS[:], 
                                    self.elevationGS[:])
        self.azimuthGS[:], self.elevationGS[:] = result
        self.x[:,0] = self.azimuthGS
        self.x[:,1] = self.elevationGS
        self.dg_daz = self.MBI.evaluate(self.x,1)[:,0]
        self.dg_del = self.MBI.evaluate(self.x,2)[:,0]

    def execute(self):
        result = self.lib.fixangles(self.n, self.azimuthGS[:], 
                                    self.elevationGS[:])
        self.azimuthGS[:], self.elevationGS[:]  = result
        self.x[:,0] = self.azimuthGS
        self.x[:,1] = self.elevationGS
        self.gain = self.MBI.evaluate(self.x)[:,0]

    def applyDer(self, arg, result):
        result['gain'][:] = np.zeros(self.n)
        if 'azimuthGS' in arg:
            result['gain'][:] += self.dg_daz * arg['azimuthGS'][:]
        if 'elevationGS' in arg:
            result['gain'][:] += self.dg_del * arg['elevationGS'][:]
        return result

    def applyDerT(self, arg, result):
        if 'gain' in arg:
            result['azimuthGS'][:] = self.dg_daz * arg['gain'][:]
            result['elevationGS'][:] = self.dg_del * arg['gain'][:]
        return result
