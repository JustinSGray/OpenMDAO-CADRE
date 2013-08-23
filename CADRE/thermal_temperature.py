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
K = 5.67051e-8 #from Thermal_Temperature0.py -> lower error
#K = 5.6704e-8 #from Thermal_Temperature.py and in original implementation

class ThermalTemperature(RK4): 

    def __init__(self, n_times): 
        super(ThermalTemperature, self).__init__()
        #self.time_step=time_step


        self.add("temperature", Array(np.zeros((5, n_times)), shape=(5, n_times), dtype=np.float,
            iotype="out", desc="temperature for the 4 fins and body over time", low=50, high =400)
        ) 

        self.add("T0", Array(273.*np.ones((5,)), shape=(5,), dtype=np.float,
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
        self.external_vars = ["exposedArea", "LOS", "P_comm"]
        self.fixed_external_vars = ["cellInstd",]
        
        # implementation of fixTemps from Thermal_Temperature.f90
        for i in range (0, n_times):
            for k in range (0, 5):
                self.temperature[k, i] = self.T0[k]
                if self.temperature[k, i] < 0:
                    self.temperature[k, i] = 0.
        
    def f_dot(self, external, state):
        
        # revised implementation from ThermalTemperature.f90
        exposedArea = external[:84]
        LOS = external[84]
        P_comm = external[85]
        cellInstd = external[86:]
        
        f = np.zeros((5, ))
        
        # Panels
        for p in range(0, 12): 
            
            # Body
            if p < 4: 
                f_i = 4
                m = m_b
                cp = cp_b
                
            # Fin
            else: 
                f_i = (p+1)%4
                m = m_f
                cp = cp_f
                
            # Cells    
            for c in range(0, 7):
                
                idat = p + c*12 
                A_exp = exposedArea[idat]
                w = cellInstd[idat]
                
                alpha = alpha_c*w + alpha_r*(1-w)
                eps = eps_c*w + eps_r*(1-w)
                
                f[f_i] += alpha * q_sol * A_exp * LOS / m / cp
                f[f_i] -= eps * K * A_T * state[f_i]**4 / m / cp

        f[4] += 4.0 * P_comm / m_b / cp_b
        
        return f
        
    
    def df_dy(self, external, state):
        
        # revised implementation from ThermalTemperature.f90
        exposedArea = external[:84]
        LOS = external[84]
        cellInstd = external[86:]
        
        dfdy = np.zeros((5,5))
        for p in range(0,12): #panels
            if p < 4: #body #lowest at p<8
                f_i = 4
                m = m_b
                cp = cp_b
            else: #fin
                f_i = (p+1)%4
                m = m_f
                cp = cp_f
            for c in range (0,7): #cells
                idat = p + c*12 #lowest at p-2
                A_exp = exposedArea[idat]
                w = cellInstd[idat]
                alpha = alpha_c*w + alpha_r*(1-w)
                eps = eps_c*w + eps_r*(1-w)

                dfdy[f_i, f_i] -= 4.0 * eps * K * A_T * state[f_i]**3 / m / cp

        return dfdy

    def df_dx(self, external, state):
        
        # revised implementation from ThermalTemperature.f90
        exposedArea = external[:84]
        LOS = external[84]
        P_comm = external[85]
        cellInstd = external[86:]
        
        dfdx = np.zeros((5,170))
        for p in range(0,12): #panels
            if p < 4: #body #lowest at p<8
                f_i = 4
                m = m_b
                cp = cp_b
            else: #fin
                f_i = (p+1)%4
                m = m_f
                cp = cp_f
            for c in range (0,7): #cells
                idat = p + c*12 #lowest at p-2
                A_exp = exposedArea[idat]
                w = cellInstd[idat]
                iA = idat
                iw = 86 + idat

                alpha = alpha_c*w + alpha_r*(1-w)
                eps = eps_c*w + eps_r*(1-w)
                
                dalpha_dw = alpha_c - alpha_r
                deps_dw = eps_c - eps_r
    
                dfdx[f_i, iA] += alpha * q_sol * LOS / m / cp
                dfdx[f_i, iw] += dalpha_dw * q_sol * A_exp * LOS / m / cp
                dfdx[f_i, iw] -= deps_dw * K * A_T * state[f_i]**4 / m / cp
                dfdx[f_i, 84] += alpha * q_sol * A_exp / m / cp
                #print f_i, iw, dfdx[f_i, iw]

        dfdx[4, 85] += 4.0 / m_b / cp_b

        return dfdx







