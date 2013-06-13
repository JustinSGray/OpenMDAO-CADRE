import numpy as np

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

class KSfunction(Component): 
    """Aggregates a number of functions to a single value via the 
    Kreisselmeier-Steinhauser Function""" 

    rho = Float(.1, iotype="in", desc="hyperparameter for the KS function")
    KS = Float(0, iotype="out", desc="value of the aggregate KS function")

    def __init__(self, n=2): 
        super(KS, self).__init__()

        self.n = n

        self.add('g',Array(zeros((n,)), size=(n,1), dtype=Float, iotype="in", 
            desc="array of function values to be aggregated"))

    def execute(self): 
        self.g_max = np.max(self.g)
        self.g_diff = self.g-self.g_max
        self.exponents = np.exp(self.rho * self.g_diff)
        self.summation = np.sum(self.exponents)
        self.KS = self.g_max + 1.0/self.rho * np.log(self.summation)

    def linearize(self): 
        """linearize around the last executed point""" 

        #use g_max, exponsnte, summation from last executed point
        dsum_dg = self.rho*self.exponents
        dKS_dsum = 1.0/self.rho/self.summation
        self.dKS_dg = dKS_dsum * dsum_dg

        dsum_drho = np.sum(self.g_diff*self.exponents)
        self.dKS_drho = dKS_dsum * dsum_drho


    def provideDer(self): 
        ins = ('g','rho')
        outs = ('KS', )
        J = np.hstack((self.dks_dg, [self.dKs_drho,]))

