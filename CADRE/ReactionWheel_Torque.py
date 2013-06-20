from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import numpy as np


class ReactionWheel_Torque(Component):

    def __init__(self, n):
        super(ReactionWheel_Torque, self).__init__()
        self.n = n
        
        self.add('T_tot', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='out'))
        
    def execute(self):
        self.T_RW[:] = self.T_tot[:]

    def applyDer(self, arg, result):
        if not result['T_RW']:
            result['T_RW'] = np.zeros((3,self.n))
        
        if 'T_tot' in arg:
            result['T_RW'][:] += arg['T_tot'][:]
        return result

    def applyDerT(self, arg, result):
        if not result['T_tot']:
            result['T_tot'] = np.zeros((3,self.n))
        
        if 'T_tot' in arg:
            result['T_tot'][:] += arg['T_RW'][:]
        return result
