import numpy as np

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import KS

class BatteryConstraints(Component): 

    ConCh = Float(0.0, iotype="out", desc="Constraint on charging rate")
    ConDS = Float(0.0, iotype="out", desc="Constraint on dischargin rate")
    ConS0 = Float(0.0, iotype="out", desc="Constraint on minimum state of charge")
    ConS1 = Float(0.0, iotype="out", desc="Constraint on maximum state of charge")

    def __init__(self, n=2): 
        self.n = n

        self.rho = 50
        self.Imin = -10.0
        self.Imax = 5.0
        self.SOC0 = 0.2
        self.SOC1 = 1.0

        self.add('I_bat', Array(np.zeros((n,)), size=(n,), iotype="in", 
            desc="Battery Current over time"))
        self.add('SOC', Array(np.zeros((n,)), size=(n,), iotype="in", 
            desc="Battery State of Charge over time"))

        self.KS_ch = KS.KSfunction()
        self.KS_ds = KS.KSfunction()
        self.KS_s0 = KS.KSfunction()
        self.KS_s1 = KS.KSfunction()

    def execute(self):

        self.ConCh = self.KS_ch(self.I_bat - self.Imax, self.rho)
        self.ConDS = self.KS_DS(self.Imin - self.I_bat, self.rho)
        self.ConS0 = self.KS_S0(self.SOC0 - self.SOC, self.rho)
        self.ConS1 = self.KS_S1(self.SOC - self.SOC1, self.rho)

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
            result['ConS0'] = -1 * np.dot(self.dS0_dg, arg['SOC'])
            result['ConS1'] = np.dot(self.dS1_dg, arg['SOC'])

    def applyDerT(self, arg, result): 
        if 'ConCh' in arg: 
            result['I_bat'] += np.dot(self.dCh_dg, arg['ConCh'])
        if 'ConDs' in arg: 
            result['I_bat'] += -1 * np.dot(self.dDs_dg, arg['ConDs'])
        if 'ConS0' in arg: 
            result['SOC'] += -1 * np.dot(self.dS0_dg, arg['ConS0'])
        if 'ConS1' in arg: 
            result['SOC'] += np.dot(self.dS1_dg, arg['ConS1']) 


        
