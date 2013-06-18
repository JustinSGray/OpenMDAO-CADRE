from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np


class Comm_VectorECI(Component):

    def __init__(self, n):
        super(Comm_VectorECI, self).__init__()
        self.n = n
        
        self.add('r_b2g_I', Array(iotype='out', shape=(3, self.n)))
        
        self.add('r_e2g_I', Array(iotype='in', shape=(3, self.n)))
        self.add('r_e2b_I', Array(iotype='in', shape=(6, self.n)))

    def execute(self):
        self.r_b2g_I = self.r_e2g_I[:] - self.r_e2b_I[:3,:]

    def applyDer(self, arg, result):
        if 'r_e2g_I' in arg and 'r_e2b_I' in arg:
            result['r_b2g_I'][:] = arg['r_e2g_I'][:]
            result['r_b2g_I'][:] -= arg['r_e2b_I'][:3,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_b2g_I' in arg:
            result['r_e2g_I'][:] = arg('r_b2g_I')[:]
            result['r_e2b_I'][:3,:] -= arg('r_b2g_I')[:]
        return result
