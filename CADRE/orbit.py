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

    def __init__(self, n_times, time_step=.01): 
        super(Orbit_Dynamics, self).__init__()
        self.time_step = time_step

        self.add('r_e2b_I', Array(1000*np.ones((6,n_times)), size=(6, n_times), 
            dtype=np.float, iotype="out"))

        self.add('r_e2b_I0', Array(np.zeros((6,)), size=(6,), iotype="in", 
            dtype=np.float))

        self.state_var = 'r_e2b_I'
        self.init_state_var = 'r_e2b_I0'

        self.dfdx = np.zeros((6, ))

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
        r = (x**2 + y**2 + z**2)**.5

        T2 = 1 - 5*z**2/r**2
        T3 = 3*z - 7*z**3/r**2
        T4 = 1 - 14*z**2/r**2 + 21*z**4/r**4
        T3z = 3*z - 0.6*r**2/z
        T4z = 4 - 28.0/3.0*z**2/r**2

        drdx = x/r
        drdy = y/r
        drdz = z/r

    def df_dx(self, external, state): 
        
        return self.dfdx




class Orbit_Initial(Component):
    altPerigee = Float(0, iotype="in", copy=None)
    altApogee = Float(0, iotype="in", copy=None)
    RAAN = Float(0, iotype="in", copy=None)
    Inc = Float(0, iotype="in", copy=None)
    argPerigee = Float(0, iotype="in", copy=None)
    trueAnomaly = Float(0, iotype="in", copy=None)
    
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
