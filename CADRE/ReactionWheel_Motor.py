from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import numpy as np


class ReactionWheel_Motor(Component):
         
    def __init__(self, n):
        super(ReactionWheel_Motor, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.RWLib').lib.RWLib
         
        self.add('J_RW', 2.8e-5)            
        
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        self.add('w_B', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        self.add('w_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))         
        
        self.add('T_m', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='out'))

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
        if not result['T_m']:
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
        if not result['T_RW']:
            result['T_RW'] = np.zeros((3,self.n))
        if not result['w_B']:
            result['w_B'] = np.zeros((3,self.n))        
        if not result['w_RW']:
            result['w_RW'] = np.zeros((3,self.n))
        
        for k in range(3):
            for j in range(3):
                if 'T_RW' in arg:
                    result['T_RW'][j,:] += self.dT_dTm[:,k,j] * arg['T_m'][k,:]
                if 'w_B' in arg:
                    result['w_B'][j,:] += self.dT_dwb[:,k,j] * arg['T_m'][k,:]
                if 'w_RW' in arg:
                    result['w_RW'][j,:] += self.dT_dh[:,k,j] * arg['T_m'][k,:] * self.J_RW
        return result