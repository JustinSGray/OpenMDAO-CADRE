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
                if 'T_RW' in arg:
                    result['T_RW'][j,:] += self.dT_dTm[:,k,j] * arg['T_m'][k,:]
                if 'w_B' in arg:
                    result['w_B'][j,:] += self.dT_dwb[:,k,j] * arg['T_m'][k,:]
                if 'w_RW' in arg:
                    result['w_RW'][j,:] += self.dT_dh[:,k,j] * arg['T_m'][k,:] * self.J_RW        
        return result
    
    
class ReactionWheel_Power(Component):

    def __init__(self, n):
        super(ReactionWheel_Power, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.RWLib').lib.RWLib
        
        self.add('w_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))

        self.add('P_RW', Array(np.ones((3,n)), size=(3,n), dtype=np.float, iotype='out')) #CHANGE TO ZEROS()??
             
    def linearize(self):
        self.dP_dw, self.dP_dT = self.lib.computejacobianp(self.n, 
                                                           self.w_RW[:], 
                                                           self.T_RW[:])
        #print "dP_dw", self.dP_dw, "dP_dT", self.dP_dT        

    def execute(self):
        self.P_RW[:] = self.lib.computep(self.n, self.w_RW[:], self.T_RW[:])
        #print "P_RW", self.P_RW        

    def applyDer(self, arg, result):
        result['P_RW'] = np.zeros((3,self.n))
        
        for k in range(3):
            if 'w_RW' in arg:
                result['P_RW'][k,:] += self.dP_dw[:,k] * arg['w_RW'][k,:]
            if 'T_RW' in arg:
                result['P_RW'][k,:] += self.dP_dT[:,k] * arg['T_RW'][k,:]
        #print "w_RW", arg['w_RW'], "T_RW", arg['T_RW'], "P_RW", result['P_RW']                
        return result

    def applyDerT(self, arg, result):
        result['w_RW'] = np.zeros((3,self.n))
        result['T_RW'] = np.zeros((3,self.n))
        
        for k in range(3):
            if 'w_RW' in arg:
                result['w_RW'][k,:] += self.dP_dw[:,k] * arg['P_RW'][k,:]
            if 'T_RW' in arg:
                result['T_RW'][k,:] += self.dP_dT[:,k] * arg['P_RW'][k,:]
        #print "w_RW", result['w_RW'], "T_RW", result['T_RW']        
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
        
        if 'T_tot' in arg:
            result['T_tot'][:] += arg['T_RW'][:]
        return result
