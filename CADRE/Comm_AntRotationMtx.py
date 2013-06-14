from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_AntRotationMtx(Component):

    def __init__(self, n):
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        self.n = n
        self.add('O_AB', Array(iotype='out', shape=(3, 3, self.n)))
        self.add('q_A', Array(iotype='in', shape=(4, self.n)))

        self.J = np.empty((self.n, 3, 3, 4))

    def linearize(self):
        self.J = self.lib.computerotmtxjacobian(self.n, self.q_A)

    def execute(self):
        self.O_AB += self.lib.computerotationmtx(self.n, self.q_A)

    def applyDer(self):
        for u in xrange(3):
            for v in xrange(3):
                for k in xrange(4):
                    self.O_AB[u,v,:] += self.J[:,u,v,k] * self.q_A[k,:]

    def applyDerT(self):
        for u in xrange(3):
            for v in xrange(3):
                for k in xrange(4):
                    self.q_A[k,:] += self.J[:,u,v,k] *  self.O_AB[u,v,:]
