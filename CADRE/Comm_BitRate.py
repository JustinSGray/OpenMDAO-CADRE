from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_BitRate(Component):

    def __init__(self, n):
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
        self.Dr += self.lib.computedr(self.n, 
                                           self.P_comm, 
                                           self.gain, 
                                           self.GSdist, 
                                           self.CommLOS)

    def applyDer(self):
        self.Dr += self.dD_dP * self.P_comm
        self.Dr += self.dD_dGt * self.gain
        self.Dr += self.dD_dS * self.GSdist
        self.Dr += self.dD_dLOS * self.CommLOS

    def applyDerT(self):
        self.P_comm += self.dD_dP * self.Dr
        self.gain += self.dD_dGt * self.Dr
        self.GSdist += self.dD_dS * self.Dr
        self.CommLOS += self.dD_dLOS * self.Dr
