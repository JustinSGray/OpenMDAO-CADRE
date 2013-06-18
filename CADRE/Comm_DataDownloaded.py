import numpy as np

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

import rk4


class Comm_DataDownloaded(rk4.RK4):

    def __init__(self, n_times, time_step=.01):
        super(Comm_DataDownloaded, self).__init__()
        self.time_step = time_step
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        
        self.add('Data', Array(iotype='out', shape=(1, n_times)))
        self.add('Dr', Array(iotype='in', shape=(n_times)))
        
        self.state_var = "Data"
        
        self.external_vars = ["Dr"]

    def f_dot(self, external, state): 
        return external[0]

    def df_dy(self, external, state): 
        return np.array([[0.]])

    def df_dx(self, external, state):   
        return np.array([[1.]])
