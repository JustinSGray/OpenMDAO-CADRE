from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import numpy as np

import rk4


class ReactionWheel_Motor(Component):
         
    def __init__(self, n):
        super(ReactionWheel_Motor, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.RWLib').lib.RWLib
         
        self.add('J_RW', 2.8e-5)            
        
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        self.add('w_B', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        self.add('w_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))         
        
        self.add('T_m', Array(np.ones((3,n)), size=(3,n), dtype=np.float, iotype='out'))

    def linearize(self): 
        self.dT_dTm, self.dT_dwb, self.dT_dh = self.lib.computejacobiant(self.n, 
                                                                         self.T_RW[:], 
                                                                         self.w_B[:], 
                                                                         self.J_RW * self.w_RW[:])
        self.dT_dTm, self.dT_dwb, self.dT_dh = self.lib.computejacobiant(self.n, 
                                                                         self.T_RW[:], 
                                                                         self.w_B[:], 
                                                                         self.J_RW * self.w_RW[:])

    def execute(self):
        self.T_m = self.lib.computet(self.n, self.T_RW[:], self.w_B[:], self.J_RW * self.w_RW[:])
        
    def applyDer(self, arg, result):
        result['T_m'] = np.zeros((3,self.n))
        
        for k in range(3):
            for j in range(3):
                if 'T_RW' in arg:
                    result['T_m'][k,:] += self.dT_dTm[:,k,j] * arg['T_RW'][j,:]
                if 'w_B' in arg:
                    result['T_m'][k,:] += self.dT_dwb[:,k,j] * arg['w_B'][j,:]
                if 'w_RW' in arg:
                    result['T_m'][k,:] += self.dT_dh[:,k,j] * arg['w_RW'][j,:] * self.J_RW
        return result
                    
    def applyDerT(self, arg, result):
        result['T_RW'] = np.zeros((3,self.n))
        result['w_B'] = np.zeros((3,self.n))        
        result['w_RW'] = np.zeros((3,self.n))
        
        for k in range(3):
            for j in range(3):
                if 'T_m' in arg:
                    result['T_RW'][j,:] += self.dT_dTm[:,k,j] * arg['T_m'][k,:]
                    result['w_B'][j,:] += self.dT_dwb[:,k,j] * arg['T_m'][k,:]
                    result['w_RW'][j,:] += self.dT_dh[:,k,j] * arg['T_m'][k,:] * self.J_RW        
        return result
    
    
class ReactionWheel_Power(Component):

    def __init__(self, n):
        super(ReactionWheel_Power, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.RWLib').lib.RWLib
        
        self.add('w_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))

        self.add('P_RW', Array(np.ones((3,n)), size=(3,n), dtype=np.float, iotype='out')) 
             
    def linearize(self):
        self.dP_dw, self.dP_dT = self.lib.computejacobianp(self.n, 
                                                           self.w_RW[:], 
                                                           self.T_RW[:])   

    def execute(self):
        self.P_RW[:] = self.lib.computep(self.n, self.w_RW[:], self.T_RW[:])       

    def applyDer(self, arg, result):
        result['P_RW'] = np.zeros((3,self.n))
        
        for k in range(3):
            if 'w_RW' in arg:
                result['P_RW'][k,:] += self.dP_dw[:,k] * arg['w_RW'][k,:]
            if 'T_RW' in arg:
                result['P_RW'][k,:] += self.dP_dT[:,k] * arg['T_RW'][k,:]
        return result

    def applyDerT(self, arg, result):
        result['w_RW'] = np.zeros((3,self.n))
        result['T_RW'] = np.zeros((3,self.n))
        
        for k in range(3):
            if 'P_RW' in arg:
                result['w_RW'][k,:] += self.dP_dw[:,k] * arg['P_RW'][k,:]
                result['T_RW'][k,:] += self.dP_dT[:,k] * arg['P_RW'][k,:]
        return result


class ReactionWheel_Torque(Component):

    def __init__(self, n):
        super(ReactionWheel_Torque, self).__init__()
        self.n = n
        
        self.add('T_tot', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='out'))
        
    def execute(self):
        self.T_RW[:] = self.T_tot[:]

    def applyDer(self, arg, result):
        result['T_RW'] = np.zeros((3,self.n))
        
        if 'T_tot' in arg:
            result['T_RW'][:] += arg['T_tot'][:]
        return result

    def applyDerT(self, arg, result):
        if not result['T_tot']:
            result['T_tot'] = np.zeros((3,self.n))
        
        if 'T_RW' in arg:
            result['T_tot'][:] += arg['T_RW'][:]
        return result


class ReactionWheel_Dynamics(rk4.RK4):
    
    def __init__(self, n_times):
        super(ReactionWheel_Dynamics, self).__init__()
        #self.time_step = time_step
        self.lib = __import__('CADRE.lib.RWLib').lib.RWLib
        
        self.add('w_B', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3,n_times)), size=(3,n_times), dtype=np.float, iotype='in'))
        
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
        self.jy[0, :] = [0., -external[2], external[1]]
        self.jy[1, :] = [external[2], 0., -external[0]]
        self.jy[2, :] = [-external[1], external[0], 0.]

        #TODO sort out unit conversion here with T_RW
        return (-external[3:]/2.8e-5 - self.jy.dot(state))
        
    
    def df_dy(self, external, state):
        return -self.jy
    
    def df_dx(self, external, state):
        self.jx = np.zeros((3,6))

        for i in xrange(3):
            self.jx[:,i] = -self.djy_dx[:,:,i].dot(state)
            self.jx[i,i+3] = -1 / self.J_RW   

        return self.jx

    