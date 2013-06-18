from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np
import scipy.sparse


class Comm_VectorSpherical(Component):

    def __init__(self, n):
        super(Comm_VectorSpherical, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('azimuthGS', Array(iotype='out', shape=(self.n)))
        self.add('elevationGS', Array(iotype='out', shape=(self.n)))
        
        self.add('r_b2g_A', Array(iotype='in', shape=(3, self.n)))

    def linearize(self):
        result = self.lib.computepositionsphericaljacobian(self.n, 3*self.n, 
                                                           self.r_b2g_A[:])
        self.Ja1, self.Ji1, self.Jj1, self.Ja2, self.Ji2, self.Jj2  = result
        self.J1 = scipy.sparse.csc_matrix((self.Ja1,(self.Ji1,self.Jj1)), 
                                          shape=(self.n,3*self.n))
        self.J2 = scipy.sparse.csc_matrix((self.Ja2,(self.Ji2,self.Jj2)), 
                                          shape=(self.n,3*self.n))
        self.J1T = self.J1.transpose()
        self.J2T = self.J2.transpose()

    def execute(self):
        azimuthGS, elevationGS = self.lib.computepositionspherical(self.n, 
                                                                   self.r_b2g_A[:])
        self.azimuthGS[:] = azimuthGS
        self.elevationGS[:] = elevationGS

    def applyDer(self, arg, result):
        if 'r_b2g_A' in arg:
            r_b2g_A = arg('r_b2g_A')[:].reshape((3*self.n),order='F')
            result['azimuthGS'][:] = self.J1.dot(r_b2g_A)
            result['elevationGS'][:] = self.J2.dot(r_b2g_A)
        return result

    def applyDerT(self, arg, result):
        if 'azimuthGS' in arg and elevationGS in arg:
            azimuthGS = arg['azimuthGS'][:]
            elevationGS = arg['elevationGS'][:]
            result['r_b2g_A'][:] = (self.J1T.dot(azimuthGS) + 
                                    self.J2T.dot(elevationGS)).reshape((3, self.n),order='F')
        return result
