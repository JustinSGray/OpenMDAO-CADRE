from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import rk4

import numpy as np


class ReactionWheel_Dynamics(rk4.RK4):
    
    def __init__(self, n_times, time_step=.01):
        super(ReactionWheel_Dynamics, self).__init__()
        self.time_step = time_step
        
        self.add('w_B', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='in'))#MAY NEED TO BE SIZE (n_times,)
        
        self.add('w_RW', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='out'))
        self.add('w_RW0', Array(np.zeros((3,)), size=(3,), dtype=np.float, iotype='in'))        
        
        self.state_var = 'w_RW'
        self.init_state_var = 'w_RW0'
        self.external_vars = ['w_B', 'T_RW']
       
        self.jy = np.zeros((3, 3))

        self.djy_dx = np.zeros((3, 3, 3))
        self.djy_dx[:,:,0] = [[0,0,0],[0,0,-1],[0,1,0]]
        self.djy_dx[:,:,1] = [[0,0,1],[0,0,0],[-1,0,0]]
        self.djy_dx[:,:,2] = [[0,-1,0],[1,0,0],[0,0,0]]

        

        self.J_RW = 2.8e-5 #unit conversion of some kind

    def f_dot(self, external, state):
        self.jy[0, :] = [0, -external[2], external[1]]
        self.jy[1, :] = [external[2], 0, -external[1]]
        self.jy[2, :] = [-external[1], external[0], 0]

        #TODO sort out unit conversion here with T_RW
        return (-external[3:]/2.8e-5 - self.jy.dot(state))
        
    
    def df_dy(self, external, state):
        return -self.jy
    
    def df_dx(self, external, state):
        self.jx = np.zeros((3,6))

        for i in xrange(3):
            self.jx[:,i] = -self.djy_dx[:,:,i].dot(state)
            self.jx[i,i+3] = -1

        return self.jx

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
            result['w_B'] = np.zeros(self.w_B.shape)
            result['T_RW'] = np.zeros(self.T_RW.shape)
            for k in range(3):
                for j in range(3):
                    result['w_B'][j,:-1] += self.Jx[1:,j,k] * arg['w_RW'][k,1:]
                    result['T_RW'][j,:-1] += self.Jx[1:,j+3,k] * arg['w_RW'][k,1:] / self.J_RW
        return result