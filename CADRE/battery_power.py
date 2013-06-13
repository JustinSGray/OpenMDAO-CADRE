from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import numpy as np

def BatteryPower(Component):


    def __init__(self, n): 
        super(BatteryPower, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.BatteryLib').lib.BatteryLib


        self.add(zeros((n,)),'P_bat',Array(dtype=Float, size=(n,), iotype="in"))
        self.add(zeros((n,)),'temperature',Array(dtype=Float, size=(n,),iotype="in"))
        self.add(zeros((n,)),'SOC',Array(dtype=Float, size=(n,),iotype="in"))
        
        self.add(zeros((n,)),'I_bat',Array(dtype=Float, size=(n,), iotype="out"))


    def linearize(self): 
        DI = self.lib.computejacobiani(self.n, self.P_bat, self.temperature, self.SOC)
        self.dI_dP, self.dI_dT, self.dI_dSOC = DI

    def execute(self):
        self.I_bat += self.lib.computei(self.n, self.P_bat,self.temperature,self.SOC)

    def applyDer(self, arg, result):
        result['I_bat'] += self.dI_dP * arg['P_bat']
        result['I_bat'] += self.dI_dT * arg['temperature']
        result['I_bat'] += self.dI_dSOC * arg['SOC']    

        return result

    def applyDerT(self, arg, result):
        result['P_bat'] += self.dI_dP * arg['I_bat']
        result['temperature'] += self.dI_dT * arg['I_bat']
        result['SOC'] += self.dI_dSOC * arg['I_bat'] 

        return result   
