import numpy as np

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

from rk4 import RK4


#constants
m_f = 0.4
m_b = 2.0
cp_f = 0.6e3
cp_b = 2.0e3
A_T = 2.66e-3

alpha_c = 0.9
alpha_r = 0.2
eps_c = 0.87
eps_r = 0.88

q_sol = 1360.0
K = 5.6704e-8

class ThermalTemperature(RK4): 

    def __init__(self, n_times, time_step=.01): 
        super(ThermalTemperature, self).__init__()

        self.add("temperature", Array(np.zeros((5,n_times)), shape=(5,n_times), dtype=np.float, 
            iotype="out", desc="temperature for the 4 fins and body over time")
        )

        self.add("T0", Array(273*np.ones((5,)), shape=(5,), dtype=np.float, 
            iotype="in", desc="initial temperatures for the 4 fins and the body")
        )

        self.add("exposedArea", Array(np.zeros((7,12,n_times)), size=(7,12,n_times), dtype=np.float, 
            iotype="in", desc="exposed area for each solar cell")
        )

        self.add("cellInstd", Array(np.ones((7,12)), size=(7,12), dtype=np.float, 
            iotype="in", desc="Cell/Radiator indication", low=0, high=1)
        )

        self.add("LOS", Array(np.zeros((n_times, )), size=(n_times, ), dtype=np.float, 
            iotype="in", desc="Line of sight to the sun", low=0, high=1)
        )

        self.add("P_comm", Array(np.ones((n_times, )), size=(n_times, ), dtype=np.float, 
            iotype="in", desc="Power required by the communication system", low=0, high=1)
            )

        self.state_var = "temperature"
        self.init_state_var = "T0"
        self.external_vars = ["exposedArea","LOS","P_comm"]
        self.fixed_external_vars = ["cellInstd",]


    def f_dot(self, external, state): 

        f = np.zeros((5, ))

        exposedArea = external[:84]
        LOS = external[84]
        P_comm = external[85]
        cellInstd = external[86:]

        for i,(A_exp,w) in enumerate(zip(exposedArea,cellInstd)): 
            p_i = i/7
            c_i = i%7

            if p_i < 4: #body
                f_i = 4
                m = m_b
                cp = cp_b

            else: #fin
                f_i = (p_i+1)%4
                m = m_f
                cp = cp_f

            alpha = alpha_c*w + alpha_r*(1-w)
            eps = eps_c*w + eps_r*(1-w)

            f[f_i] += alpha * q_sol * A_exp * LOS / m / cp
            f[f_i] -= eps * K * A_T * state[f_i]**4 / m / cp

        f[4] += 4.0 * P_comm / m_b / cp_b

        return f








