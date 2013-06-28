from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import rk4

import numpy as np


class ReactionWheel_Dynamics(rk4.RK4):
    
    def __init__(self, n_times, time_step=.01):
        super(ReactionWheel_Dynamics, self).__init__()
        self.time_step = time_step
        self.lib = __import__('CADRE.lib.RWLib').lib.RWLib

        #self.add('J_RW', Array(2.8e-5*np.ones((1,)), size=(1,), dtype=np.float, iotype='in'))
        self.add('J_RW', Float(2.8e-5, iotype='in'))
        
        self.add('w_B', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='in'))#MAY NEED TO BE SIZE (n_times,)
        
        self.add('w_RW', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='out'))
        self.add('w_RW0', Array(np.zeros((3,)), size=(3,), dtype=np.float, iotype='in'))        
        
        self.state_var = 'w_RW'
        self.init_state_var = 'w_RW0'#THIS MAY NEED TO BE w_B OR T_RW
        self.external_vars = ['w_B', 'T_RW']#CHANGE BASED ON PREVIOUS LINE
       
    def f_dot(self, external, state):
        return external[0]
        
    
    def df_dy(self, external, state):
        return np.array([[0.]])
    
    def df_dx(self, external, state):
        return np.array([[1.]])
    
    def compute(self, val):
        self.dat[:3,:] = val('w_B')[:]
        self.dat[3:,:] = val('T_RW')[:] / self.J_RW

    def applyJext(self, arg, result):
        result['w_RW'] = np.zeros((self.w_RW.shape))
        for k in range(3):
            for j in range(3):
                if 'w_B' in arg:
                    result['w_RW'][k,1:] += self.Jx[1:,j,k] * arg['w_B'][j,:-1]
                if 'T_RW' in arg:
                    result['w_RW'][k,1:] += self.Jx[1:,j+3,k] * arg['T_RW'][j,:-1] / self.J_RW                    
        return result

    def applyJextT(self, arg, result):
        if 'w_RW' in arg:
            for k in range(3):
                for j in range(3):
                    result['w_B'][j,:-1] += self.Jx[1:,j,k] * arg['w_RW'][k,1:]
                    result['T_RW'][j,:-1] += self.Jx[1:,j+3,k] * arg['w_RW'][k,1:] / self.J_RW
        return result