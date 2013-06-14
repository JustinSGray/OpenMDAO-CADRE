import numpy as np

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import KS

class BatteryPower(Component): 

    def __init__(self, n): 
        super(BatteryPower, self).__init__()

        self.n = n 

        self.sigma = 1e-10
        self.eta = 0.99
        self.Cp = 2900*0.001*3600
        self.IR = 0.9
        self.T0 = 293
        self.alpha = log(1/1.1**5)

        self.add('SOC', Array(zeros((1, n)), size=(1, n), dtype=Float, iotype="in"))
        self.add('temperature', Array(zeros((5, n)), size=(n, ), dtype=Float, iotype="in"))
        self.add('P_bat', Array(zeros((n, )), size=(n, ), dtype=Float, iotype="in"))


        self.add('I_bat', Array(zeros((n, )), size=(n, ), dtype=Float, iotype="out"), 
            desc="Batter Current over Time")

    def execute(self): 
        self.exponential = (2 - np.exp(self.alpha*(self.temperaturese-self.T0)/self.T0))
        self.voc = 3 + (np.exp(self.SOC)-1) / (np.e-1)
        self.V = self.IR * voc * self.exponential
        self.I_bat = self.P_bat/V

    def linearize(self): 
        #dI_dP
        dV_dvoc = self.IR * self.exponential
        dV_dT = -self.V*self.alpha/self.T0
        dVoc_dSOC = np.exp(self.SOC) / (np.e-1)

        self.dI_dP = 1.0 / self.V
        tmp = -self.P_bat/(self.V**2)
        self.dI_dT = tmp * dV_dT
        self.dI_dSOC = tmp * dV_dvoc * dVoc_dSOC

    def applyDer(self, arg, result):
        if 'P_bat' in arg: 
            result['I_bat'] = self.dI_dP * arg['P_bat'] 
        if 'temperature' in arg:  
            result['I_bat'] += self.dI_dT * arg['temperature'][4,:]
        if 'SOC' in arg:     
            result['I_bat'] += self.dI_dSOC * arg['SOC'][0,:]

            
        


class BatteryConstraints(Component): 

    ConCh = Float(0.0, iotype="out", desc="Constraint on charging rate")
    ConDs = Float(0.0, iotype="out", desc="Constraint on dischargin rate")
    ConS0 = Float(0.0, iotype="out", desc="Constraint on minimum state of charge")
    ConS1 = Float(0.0, iotype="out", desc="Constraint on maximum state of charge")

    def __init__(self, n=2): 
        super(BatteryConstraints, self).__init__()
        self.n = n

        self.rho = 50
        self.Imin = -10.0
        self.Imax = 5.0
        self.SOC0 = 0.2
        self.SOC1 = 1.0

        self.add('I_bat', Array(np.zeros((n,)), size=(n,), iotype="in", 
            desc="Battery Current over time"))
        self.add('SOC', Array(np.zeros((1, n)), size=(1, n), iotype="in", 
            desc="Battery State of Charge over time"))

        self.KS_ch = KS.KSfunction()
        self.KS_ds = KS.KSfunction()
        self.KS_s0 = KS.KSfunction()
        self.KS_s1 = KS.KSfunction()

    def execute(self):
        self.ConCh = self.KS_ch.compute(self.I_bat - self.Imax, self.rho)
        self.ConDs = self.KS_ds.compute(self.Imin - self.I_bat, self.rho)
        self.ConS0 = self.KS_s0.compute(self.SOC0 - self.SOC, self.rho)
        self.ConS1 = self.KS_s1.compute(self.SOC - self.SOC1, self.rho)

    def linearize(self): 
        
        self.dCh_dg, self.dCh_drho = self.KS_ch.derivatives()
        self.dDs_dg, self.dDs_drho = self.KS_ds.derivatives()
        self.dS0_dg, self.dS0_drho = self.KS_s0.derivatives()
        self.dS1_dg, self.dS1_drho = self.KS_s1.derivatives()
        
    def applyDer(self, arg, result):
        if 'I_bat' in arg: 
            result['ConCh'] = np.dot(self.dCh_dg, arg['I_bat'])
            result['ConDs'] = -1 * np.dot(self.dDs_dg, arg['I_bat'])
        if 'SOC' in arg:
            result['ConS0'] = -1* np.dot(self.dS0_dg, arg['SOC'][0,:])
            result['ConS1'] = np.dot(self.dS1_dg, arg['SOC'][0,:])

        return result

    def applyDerT(self, arg, result): 
        if 'ConCh' in arg: 
            result['I_bat'] = np.dot(self.dCh_dg, arg['ConCh'])
        if 'ConDs' in arg: 
            result['I_bat'] += -1 * np.dot(self.dDs_dg, arg['ConDs'])
        if 'ConS0' in arg: 
            result['SOC'] = -1 * np.dot(self.dS0_dg, arg['ConS0'])
        if 'ConS1' in arg: 
            result['SOC'] += np.dot(self.dS1_dg, arg['ConS1']) 

        return result
            


        
