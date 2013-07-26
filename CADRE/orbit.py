from openmdao.lib.datatypes.api import Float, Array
from openmdao.main.api import Component
import numpy as np

import rk4



mu = 398600.44
Re = 6378.137
J2 = 1.08264e-3
J3 = -2.51e-6
J4 = -1.60e-6

C1 = -mu
C2 = -1.5*mu*J2*Re**2
C3 = -2.5*mu*J3*Re**3
C4 = 1.875*mu*J4*Re**4

class Orbit_Dynamics(rk4.RK4): 

    def __init__(self, n_times, time_step=28.819212808539024): 
        super(Orbit_Dynamics, self).__init__()
        self.time_step = time_step

        self.add('r_e2b_I', Array(1000*np.ones((6,n_times)), size=(6, n_times), 
            dtype=np.float, iotype="out"))

        self.add('r_e2b_I0', Array(np.zeros((6,)), size=(6,), iotype="in", 
            dtype=np.float))

        self.state_var = 'r_e2b_I'
        self.init_state_var = 'r_e2b_I0'

        self.dfdx = np.zeros((6, 1))

    def f_dot(self, external, state):

        x = state[0]
        y = state[1]
        z = state[2] if abs(state[2]) > 1e-15 else 1e-5
        r = (x**2 + y**2 + z**2)**.5

        T2 = 1 - 5*z**2/r**2
        T3 = 3*z - 7*z**3/r**2
        T4 = 1 - 14*z**2/r**2 + 21*z**4/r**4
        T3z = 3*z - 0.6*r**2/z
        T4z = 4 - 28.0/3.0*z**2/r**2

        f_dot = np.zeros((6,))
        f_dot[0:3] = state[3:]
        f_dot[3:] = state[0:3]*(C1/r**3 + C2/r**5*T2 + C3/r**7*T3 + C4/r**7*T4)
        f_dot[5] += z*(C2/r**5*2 + C3/r**7*T3z + C4/r**7*T4z)

        return f_dot

    def df_dy(self, external, state): 
        
        x = state[0]
        y = state[1]
        z = state[2] if abs(state[2]) > 1e-15 else 1e-5
        r = (x**2 + y**2 + z**2)**0.5

        T2 = 1 - 5*z**2/r**2
        T3 = 3*z - 7*z**3/r**2
        T4 = 1 - 14*z**2/r**2 + 21*z**4/r**4
        T3z = 3*z - 0.6*r**2/z
        T4z = 4 - 28.0/3.0*z**2/r**2

        drdx = x/r
        drdy = y/r
        drdz = z/r
        
        T2 = 1 - 5*z**2/r**2
        T3 = 3*z - 7*z**3/r**2
        T4 = 1 - 14*z**2/r**2 + 21*z**4/r**4
        T3z = 3*z - 0.6*r**2/z
        T4z = 4 - 28.0/3.0*z**2/r**2
        
        dT2_dx = (10*z**2)/(r**3)*drdx
        dT2_dy = (10*z**2)/(r**3)*drdy
        dT2_dz = (10*z**2)/(r**3)*drdz - 10.*z/r**2

        dT3_dx = 14*z**3/r**3*drdx
        dT3_dy = 14*z**3/r**3*drdy
        dT3_dz = 14*z**3/r**3*drdz - 21.*z**2/r**2 + 3

        dT4_dx = 28*z**2/r**3*drdx - 84.*z**4/r**5*drdx
        dT4_dy = 28*z**2/r**3*drdy - 84.*z**4/r**5*drdy
        dT4_dz = 28*z**2/r**3*drdz - 84.*z**4/r**5*drdz - 28*z/r**2 + 84*z**3/r**4

        dT3z_dx = -1.2*r/z*drdx
        dT3z_dy = -1.2*r/z*drdy
        dT3z_dz = -1.2*r/z*drdz + 0.6*r**2/z**2 + 3

        dT4z_dx = 56.0/3.0*z**2/r**3*drdx
        dT4z_dy = 56.0/3.0*z**2/r**3*drdy
        dT4z_dz = 56.0/3.0*z**2/r**3*drdz - 56.0/3.0*z/r**2 
        
        eye = np.identity(3)
        
        dfdy = np.zeros((6,6))
        dfdy[0:3,3:] += eye
        
        dfdy[3:,:3] = dfdy[3:,:3] + eye*(C1/r**3 + C2/r**5*T2 + C3/r**7*T3 + C4/r**7*T4)
        dfdy[3:,0] = dfdy[3:,0] + drdx*state[:3]*(-3*C1/r**4 - 5*C2/r**6*T2 - 7*C3/r**8*T3 - 7*C4/r**8*T4)
        dfdy[3:,1] = dfdy[3:,1] + drdy*state[:3]*(-3*C1/r**4 - 5*C2/r**6*T2 - 7*C3/r**8*T3 - 7*C4/r**8*T4)
        dfdy[3:,2] = dfdy[3:,2] + drdz*state[:3]*(-3*C1/r**4 - 5*C2/r**6*T2 - 7*C3/r**8*T3 - 7*C4/r**8*T4)
        dfdy[3:,0] = dfdy[3:,0] + state[:3]*(C2/r**5*dT2_dx + C3/r**7*dT3_dx + C4/r**7*dT4_dx)
        dfdy[3:,1] = dfdy[3:,1] + state[:3]*(C2/r**5*dT2_dy + C3/r**7*dT3_dy + C4/r**7*dT4_dy)
        dfdy[3:,2] = dfdy[3:,2] + state[:3]*(C2/r**5*dT2_dz + C3/r**7*dT3_dz + C4/r**7*dT4_dz)
        dfdy[5,0] = dfdy[5,0] + drdx*z*(-5*C2/r**6*2 - 7*C3/r**8*T3z - 7*C4/r**8*T4z)
        dfdy[5,1] = dfdy[5,1] + drdy*z*(-5*C2/r**6*2 - 7*C3/r**8*T3z - 7*C4/r**8*T4z)
        dfdy[5,2] = dfdy[5,2] + drdz*z*(-5*C2/r**6*2 - 7*C3/r**8*T3z - 7*C4/r**8*T4z)
        dfdy[5,0] = dfdy[5,0] + z*(C3/r**7*dT3z_dx + C4/r**7*dT4z_dx)
        dfdy[5,1] = dfdy[5,1] + z*(C3/r**7*dT3z_dy + C4/r**7*dT4z_dy)
        dfdy[5,2] = dfdy[5,2] + z*(C3/r**7*dT3z_dz + C4/r**7*dT4z_dz)
        dfdy[5,2] = dfdy[5,2] + (C2/r**5*2 + C3/r**7*T3z + C4/r**7*T4z)


        return dfdy

    def df_dx(self, external, state): 
        
        return self.dfdx




class Orbit_Initial(Component):

    altPerigee = Float(500., iotype="in", copy=None)
    altApogee = Float(500., iotype="in", copy=None)
    RAAN = Float(66.279, iotype="in", copy=None)
    Inc = Float(82.072, iotype="in", copy=None)
    argPerigee = Float(0, iotype="in", copy=None)
    trueAnomaly = Float(337.987, iotype="in", copy=None)
    
    def __init__(self):
        super(Orbit_Initial, self).__init__()
        self.add('r_e2b_I0', Array(np.ones((6,)), size=(6,), dtype=np.float, iotype='out'))
        
    def compute(self, altPerigee, altApogee, RAAN, Inc, argPerigee, trueAnomaly):
        Re=6378.137
        mu=398600.44
        
        def S(v):
            S = np.zeros((3,3),complex)
            S[0,:] = [0, -v[2], v[1]]
            S[1,:] = [v[2], 0, -v[0]]
            S[2,:] = [-v[1], v[0], 0]
            return S
            
        def getRotation(axis,angle):
            R = np.eye(3,dtype=complex) + S(axis)*np.sin(angle) + (1 - np.cos(angle)) * (np.outer(axis,axis) - np.eye(3,dtype=complex))
            return R
            
        d2r = np.pi/180.0
        r_perigee = Re + altPerigee
        r_apogee = Re + altApogee
        e = (r_apogee-r_perigee)/(r_apogee+r_perigee)
        a = (r_perigee+r_apogee)/2
        p = a*(1-e**2)
        h = np.sqrt(p*mu)

        rmag0 = p/(1+e*np.cos(d2r*trueAnomaly))
        r0_P = np.array([rmag0*np.cos(d2r*trueAnomaly), rmag0*np.sin(d2r*trueAnomaly), 0], complex)
        v0_P = np.array([-np.sqrt(mu/p)*np.sin(d2r*trueAnomaly), np.sqrt(mu/p)*(e+np.cos(d2r*trueAnomaly)), 0], complex)

        O_IP = np.eye(3, dtype=complex)
        O_IP = np.dot(O_IP, getRotation(np.array([0,0,1]),RAAN*d2r))
        O_IP = np.dot(O_IP, getRotation(np.array([1,0,0]),Inc*d2r))
        O_IP = np.dot(O_IP, getRotation(np.array([0,0,1]),argPerigee*d2r))

        r0_ECI = np.dot(O_IP,r0_P)
        v0_ECI = np.dot(O_IP,v0_P)

        return r0_ECI, v0_ECI
        
    def linearize(self):
        h = 1e-16
        ih = complex(0,h)
        v = np.zeros(6,complex)
        v[:] = [self.altPerigee, self.altApogee, self.RAAN, self.Inc, self.argPerigee, self.trueAnomaly]
        self.J = np.zeros((6,6))
        for i in range(6):
            v[i] += ih
            r0_ECI, v0_ECI = self.compute(v[0], v[1], v[2], v[3], v[4], v[5])
            v[i] -= ih
            self.J[:3,i] = r0_ECI.imag/h
            self.J[3:,i] = v0_ECI.imag/h
    
    def execute(self):
        import pdb; pdb.set_trace()
        r0_ECI, v0_ECI = self.compute(self.altPerigee, self.altApogee, self.RAAN, self.Inc, self.argPerigee, self.trueAnomaly)
        self.r_e2b_I0[:3] = r0_ECI.real
        self.r_e2b_I0[3:] = v0_ECI.real
        
    def applyDer(self, arg, result):
        if not result['r_e2b_I0']:
            result['r_e2b_I0'] = np.zeros(6)

        J = self.J
        result['r_e2b_I0'][:3] = J[:3,0]*arg['altPerigee'] + J[:3,1]*arg['altApogee'] + J[:3,2]*arg['RAAN'] + \
        J[:3,3]*arg['Inc'] + J[:3,4]*arg['argPerigee'] + J[:3,5]*arg['trueAnomaly']
        result['r_e2b_I0'][3:] = J[3:,0]*arg['altPerigee'] + J[3:,1]*arg['altApogee'] + J[3:,2]*arg['RAAN'] + \
        J[3:,3]*arg['Inc'] + J[3:,4]*arg['argPerigee'] + J[3:,5]*arg['trueAnomaly']
    
        return result
        
    def applyDerT(self, arg, result):
        J = self.J
        if 'r_e2b_I0' in arg:
            result['altPerigee'] = sum(J[:3,0]*arg['r_e2b_I0'][:3]) + sum(J[3:,0]*arg['r_e2b_I0'][3:])
            result['altApogee'] = sum(J[:3,1]*arg['r_e2b_I0'][:3]) + sum(J[3:,1]*arg['r_e2b_I0'][3:])
            result['RAAN'] = sum(J[:3,2]*arg['r_e2b_I0'][:3]) + sum(J[3:,2]*arg['r_e2b_I0'][3:])
            result['Inc'] = sum(J[:3,3]*arg['r_e2b_I0'][:3]) + sum(J[3:,3]*arg['r_e2b_I0'][3:])
            result['argPerigee'] = sum(J[:3,4]*arg['r_e2b_I0'][:3]) + sum(J[3:,4]*arg['r_e2b_I0'][3:])
            result['trueAnomaly'] = sum(J[:3,5]*arg['r_e2b_I0'][:3]) + sum(J[3:,5]*arg['r_e2b_I0'][3:])
        return result
