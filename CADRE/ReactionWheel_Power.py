from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import numpy as np

class ReactionWheel_Power(Component):

    def __init__(self, n):
        super(ReactionWheel_Power, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.RWLib').lib.RWLib

        self.add('w_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        self.add('T_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='in'))
        
        self.add('P_RW', Array(np.zeros((3,n)), size=(3,n), dtype=np.float, iotype='out'))

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
            #print 'T_RW',self.T_RW,'\nP_RW',self.P_RW,'\nw_RW',self.w_RW 
            #print self.dP_dw[:,k] * arg['w_RW'][k,:]
            #print self.dP_dT[:,k] * arg['T_RW'][k,:]            
        #print result['P_RW']
        #exit()
        return result

    def applyDerT(self, arg, result):
        #if not result['w_RW']:
        result['w_RW'] = np.zeros((3,self.n))
        #if not result['T_RW']:
        result['T_RW'] = np.zeros((3,self.n))        
        
        for k in range(3):
            if 'w_RW' in arg:
                result['w_RW'][k,:] += self.dP_dw[:,k] * arg['P_RW'][k,:]
            if 'T_RW' in arg:
                result['T_RW'][k,:] += self.dP_dT[:,k] * arg['P_RW'][k,:]
        return result
