from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_AntRotation(Component):
    antAngle = Float(0, iotype="in", copy=None)

    def __init__(self, n):
        self.add('q_A', Array(iotype='out', shape=(4,n)))
        self.add('dq_dt', Array(iotype='out', shape=(4,n)))
        self.n = n

    def linearize(self):
        rt2 = np.sqrt(2)
        self.dq_dt[0,:] = - np.sin(self.antAngle/2.) / 2.
        self.dq_dt[1,:] = np.cos(self.antAngle/2.) / rt2 / 2.
        self.dq_dt[2,:] = - np.cos(self.antAngle/2.) / rt2 / 2.
        self.dq_dt[3,:] = 0.0

    def execute(self):
        rt2 = np.sqrt(2)
        self.q_A[0,:] += np.cos(self.antAngle/2.)
        self.q_A[1,:] += np.sin(self.antAngle/2.) / rt2
        self.q_A[2,:] += - np.sin(self.antAngle/2.) / rt2
        self.q_A[3,:] += 0.0

    def applyDer(self):
        self.q_A += self.antAngle * self.dq_dt

    def applyDerT(self):
        self.antAngle += np.sum([self.dq_dt * np.sum(self.q_A[i,:]) for i in xrange(4)])
