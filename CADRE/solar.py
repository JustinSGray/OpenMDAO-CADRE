from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

from MBI import MBI

import numpy as np


class Solar_ExposedArea(Component):

    '''
    Exposed area calculation for a given solar cell
       p: panel ID [0,11]
       c: cell ID [0,6]
       a: fin angle [0,90]
       z: azimuth [0,360]
       e: elevation [0,180]
       LOS: line of sight with the sun [0,1]
    '''
    finAngle = Float(0., iotype="in", copy=None)
    def __init__(self, n):
        super(Solar_ExposedArea, self).__init__()
        self.n = n
        self.add('raw1', Array(dtype=np.float))
        self.add('raw2', Array(dtype=np.float))
        raw1 = np.genfromtxt('CADRE/data/Solar/Area10.txt')
        raw2 = np.loadtxt("CADRE/data/Solar/Area_all.txt")
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib

        self.add('nc', 7)
        self.add('np', 12)


        self.add('azimuth', Array(np.zeros((5,)), size=(5,), dtype=np.float, iotype='in'))
        self.add('elevation', Array(np.zeros((5,)), size=(5,), dtype=np.float, iotype='in'))
        
        self.add('exposedArea', Array(np.zeros((self.nc,self.np,self.n)), size=(self.nc,self.np,self.n), 
                                      dtype=np.float, iotype='out', low=-5e-3, high=1.834e-1))

        self.add('na', 10)
        self.add('nz', 73)
        self.add('ne', 37)
        
        angle = np.zeros(self.na)
        azimuth = np.zeros(self.nz)
        elevation = np.zeros(self.ne)

        index = 0
        for i in range(self.na):
            angle[i] = raw1[index]
            index += 1
        for i in range(self.nz):
            azimuth[i] = raw1[index]
            index += 1
        index -= 1
        azimuth[self.nz-1] = 2*np.pi
        for i in range(self.ne):
            elevation[i] = raw1[index]
            index += 1

        angle[0] = 0.0
        angle[-1] = np.pi/2.0
        azimuth[0] = 0.0
        azimuth[-1] = 2*np.pi
        elevation[0] = 0.0
        elevation[-1] = np.pi

        counter = 0
        data = np.zeros((self.na, self.nz, self.ne, self.np*self.nc))
        for p in range(self.np):
            for c in range(self.nc):
                index = 119
                for i in range(self.na):
                    for j in range(self.nz):
                        for k in range(self.ne):
                            data[i,j,k,counter] = raw2[7*p+c][index]
                            index += 1
                counter += 1
        self.MBI = MBI(data,[angle,azimuth,elevation],[4,10,8],[4,4,4])
        self.x = np.zeros((self.n,3),order='F')
        self.exposedArea = np.zeros((self.nc,self.np,self.n),order='F')
        self.Js = [None for i in range(3)]
    
    def setx(self):
        result = self.lib.fixangles(self.n, self.azimuth[:], self.elevation[:])
        self.x[:,0] = self.finAngle
        self.x[:,1] = result[0]
        self.x[:,2] = result[1]
    

    def linearize(self):
        self.setx()
        for i in range(3):
            self.Js[i] = self.MBI.evaluate(self.x,1+i)

    def execute(self):
        self.setx()
        self.P = self.MBI.evaluate(self.x)
        for c in range(7):
            for p in range(12):
                self.exposedArea[c,p,:] += self.P[:,7*p+c]

    def applyDer(self, arg, result):
        result['exposedArea'] = np.zeros((self.nc,self.np,self.n))
        
        for c in range(7):
            for p in range(12):
                if 'finAngle' in arg:
                    result['exposedArea'][c,p,:] += self.Js[0][:,7*p+c]*arg['finAngle']
                if 'azimuth' in arg:
                    result['exposedArea'][c,p,:] += self.Js[1][:,7*p+c]*arg['azimuth'][:]
                if 'elevation' in arg:
                    result['exposedArea'][c,p,:] += self.Js[2][:,7*p+c]*arg['elevation'][:]
        return result

    def applyDerT(self, arg, result):
        result['finAngle'] = 0.0
        result['azimuth'] = np.zeros((5,))
        result['elevation'] = np.zeros((5,))
        for c in range(7):
            for p in range(12):
                if 'finAngle' in arg:
                    result['finAngle'] += np.dot(self.Js[0][:,7*p+c], arg['exposedArea'][c,p,:])
                if 'azimuth' in arg:
                    result['azimuth'][:] += self.Js[1][:,7*p+c]*arg['exposedArea'][c,p,:]
                if 'elevation' in arg:
                    result['elevation'][:] += self.Js[2][:,7*p+c]*arg['exposedArea'][c,p,:]
        return result
