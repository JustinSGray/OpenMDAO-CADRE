from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_Distance(Component):

    def __init__(self, n):
        super(Comm_Distance, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        self.add('GSdist', Array(iotype='out', shape=(self.n)))
        self.add('r_b2g_A', Array(iotype='out', shape=(3, self.n)))

    def linearize(self):
        self.J = self.lib.computejacobiand(self.n, self.r_b2g_A)

    def execute(self):
        self.GSdist[:] = self.lib.computed(self.n, self.r_b2g_A)

    def applyDer(self, arg, result):
        if 'r_b2g_A' in arg:
            result['GSdist'][:]  = np.zeros(self.n)
            for k in xrange(3):
                result['GSdist'][:] += self.J[:,k] * arg['r_b2g_A'][k,:]
        return result

    def applyDerT(self, arg, result):
        if 'GSdist' in arg:
            result['r_b2g_A'] = np.zeros((3, self.n))
            for k in xrange(3):
                result['r_b2g_A'][k,:] += self.J[:,k] * arg['GSdist'][:]
        return result
