from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component
import numpy as np
import scipy.sparse
import MBI

import rk4


class Comm_DataDownloaded(rk4.RK4):

    def __init__(self, n_times):
        super(Comm_DataDownloaded, self).__init__()
        #self.time_step = time_step
        
        self.add('Data0', Array([0.0], iotype='in', shape=(1,)))
        self.add('Data', Array(np.zeros((1, n_times)), iotype='out', 
                               shape=(1, n_times)))
        
        self.add('Dr', Array(np.zeros(n_times), iotype='in', shape=(n_times,)))
        
        
        self.state_var = "Data"
        self.init_state_var = "Data0"
        self.external_vars = ["Dr"]

        self.dfdy = np.array([[0.]])
        self.dfdx = np.array([[1.]])

    def f_dot(self, external, state):
        return external[0]

    def df_dy(self, external, state):
        return self.dfdy

    def df_dx(self, external, state):
        return self.dfdx

   

class Comm_AntRotation(Component):
    antAngle = Float(0., iotype="in", copy=None)

    def __init__(self, n):
        super(Comm_AntRotation, self).__init__()
        self.add('q_A', Array(np.zeros((4, n)), iotype='out', shape=(4,n)))
        self.dq_dt = np.zeros(4)
        self.n = n

    def linearize(self):
        rt2 = np.sqrt(2)
        self.dq_dt[0] = - np.sin(self.antAngle/2.) / 2.
        self.dq_dt[1] = np.cos(self.antAngle/2.) / rt2 / 2.
        self.dq_dt[2] = - np.cos(self.antAngle/2.) / rt2 / 2.
        self.dq_dt[3] = 0.0

    def execute(self):
        rt2 = np.sqrt(2)
        self.q_A[0,:] = np.cos(self.antAngle/2.)
        self.q_A[1,:] = np.sin(self.antAngle/2.) / rt2
        self.q_A[2,:] = - np.sin(self.antAngle/2.) / rt2
        self.q_A[3,:] = 0.0

    def applyDer(self, arg, result):
        if 'antAngle' in arg:
            result['q_A'] = np.zeros((4, self.n))
            for k in xrange(4):
                result['q_A'][k,:] += self.dq_dt[k] * arg['antAngle']
        return result

    def applyDerT(self, arg, result):
        if 'q_A' in arg:
            result['antAngle'] = 0.
            for k in xrange(4):
                result['antAngle'] += self.dq_dt[k] * np.sum(arg['q_A'][k,:])
        return result


class Comm_AntRotationMtx(Component):

    def __init__(self, n):
        super(Comm_AntRotationMtx, self).__init__()
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        self.n = n
        self.add('O_AB', Array(np.zeros((3, 3, self.n)), 
                               iotype='out', shape=(3, 3, self.n)))
        self.add('q_A', Array(np.zeros((4, self.n)), iotype='in', 
                              shape=(4, self.n)))

        self.J = np.empty((self.n, 3, 3, 4))

    def linearize(self):
        self.J = self.lib.computerotmtxjacobian(self.n, self.q_A)

    def execute(self):
        self.O_AB = self.lib.computerotationmtx(self.n, self.q_A)

    def applyDer(self, arg, result):
        if 'q_A' in arg:
            result['O_AB'] = np.zeros((3, 3, self.n))
            for u in xrange(3):
                for v in xrange(3):
                    for k in xrange(4):
                        result['O_AB'][u,v,:] += self.J[:,u,v,k] * arg['q_A'][k,:]
        return result  

    def applyDerT(self, arg, result):
        if 'O_AB' in arg:
            result['q_A'] = np.zeros((4, self.n))
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['q_A'][k,:] += self.J[:,u,v,k] * arg['O_AB'][u,v,:]
        return result
    
    
class Comm_BitRate(Component):

    def __init__(self, n):
        super(Comm_BitRate, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib

        self.add('Dr', Array(np.zeros(self.n), iotype='out', shape=(self.n,)))

        self.add('P_comm', Array(np.zeros(self.n), iotype='in', shape=(self.n,)))
        self.add('gain', Array(np.zeros(self.n), iotype='in', shape=(self.n,)))
        self.add('GSdist', Array(np.zeros(self.n), iotype='in', shape=(self.n,)))
        self.add('CommLOS', Array(np.zeros(self.n), iotype='in', shape=(self.n,)))

    def linearize(self):
        response = self.lib.computejacobiandr(self.n, 
                                              self.P_comm, 
                                              self.gain, 
                                              self.GSdist, 
                                              self.CommLOS)
        self.dD_dP, self.dD_dGt, self.dD_dS, self.dD_dLOS = response
        
    def execute(self):
        self.Dr = self.lib.computedr(self.n, 
                                           self.P_comm, 
                                           self.gain, 
                                           self.GSdist, 
                                           self.CommLOS)

    def applyDer(self, arg, result):
        if 'P_comm' in arg:
            result['Dr'] = self.dD_dP * arg['P_comm'][:]
        if 'gain' in arg:
            result['Dr'] += self.dD_dGt * arg['gain'][:]
        if 'GSdist' in arg:
            result['Dr'] += self.dD_dS * arg['GSdist'][:]
        if 'CommLOS' in arg:
            result['Dr'] += self.dD_dLOS * arg['CommLOS'][:]
        return result

    def applyDerT(self, arg, result):
        if 'Dr' in arg:
            result['P_comm'] = self.dD_dP * arg['Dr']
            result['gain'] = self.dD_dGt * arg['Dr']
            result['GSdist'] = self.dD_dS * arg['Dr']
            result['CommLOS'] = self.dD_dLOS * arg['Dr']
        return result
    
    
class Comm_Distance(Component):

    def __init__(self, n):
        super(Comm_Distance, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        self.add('GSdist', Array(np.zeros(self.n), iotype='out', 
                                 shape=(self.n,)))
        self.add('r_b2g_A', Array(np.zeros((3, self.n)), iotype='in', 
                                  shape=(3, self.n)))

    def linearize(self):
        self.J = self.lib.computejacobiand(self.n, self.r_b2g_A)

    def execute(self):
        self.GSdist = self.lib.computed(self.n, self.r_b2g_A)

    def applyDer(self, arg, result):
        if 'r_b2g_A' in arg:
            result['GSdist']  = np.zeros(self.n)
            for k in xrange(3):
                result['GSdist'] += self.J[:,k] * arg['r_b2g_A'][k,:]
        return result

    def applyDerT(self, arg, result):
        if 'GSdist' in arg:
            result['r_b2g_A'] = np.zeros((3, self.n))
            for k in xrange(3):
                result['r_b2g_A'][k,:] += self.J[:,k] * arg['GSdist'][:]
        return result


class Comm_EarthsSpin(Component):

    def __init__(self, n):
        super(Comm_EarthsSpin, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        
        self.add('q_E', Array(np.zeros((4, self.n)), iotype='out', 
                              shape=(4, self.n)))
        self.add('t', Array(np.zeros(self.n), iotype='in', shape=(self.n,)))

    def linearize(self):
        self.dq_dt = self.lib.computejacobianqe(self.n, self.t)

    def execute(self):
        self.q_E = self.lib.computeqe(self.n, self.t)

    def applyDer(self, arg, result):
        if 't' in arg:
            result['q_E'] = np.zeros((4, self.n))
            for k in range(4):
                result['q_E'][k,:] += self.dq_dt[:,k] * arg['t'][:]
        return result

    def applyDerT(self, arg, result):
        if 'q_E' in arg:
            result['t'] = np.zeros(self.n)
            for k in range(4):
                result['t'][:] += self.dq_dt[:,k] * arg['q_E'][k,:]
        return result


class Comm_EarthsSpinMtx(Component):

    def __init__(self, n):
        super(Comm_EarthsSpinMtx, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('q_E', Array(np.zeros((4, self.n)), iotype='in', 
                              shape=(4, self.n)))
        self.add('O_IE', Array(np.zeros((3, 3, self.n)), iotype='out', 
                               shape=(3, 3, self.n)))

    def linearize(self):
        self.J = self.lib.computerotmtxjacobian(self.n, self.q_E)

    def execute(self):
        self.O_IE = self.lib.computerotationmtx(self.n, self.q_E)

    def applyDer(self, arg, result):
        if 'q_E' in arg:
            result['O_IE'] = np.zeros((3, 3, self.n))
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['O_IE'][u,v,:] += self.J[:,u,v,k] * arg['q_E'][k,:]
        return result

    def applyDerT(self, arg, result):
        if 'O_IE' in arg:
            result['q_E'] = np.zeros((4, self.n))
            for u in range(3):
                for v in range(3):
                    for k in range(4):
                        result['q_E'][k,:] += self.J[:,u,v,k] * arg['O_IE'][u,v,:]
        return result


class Comm_GainPattern(Component):

    def __init__(self, n):
        super(Comm_GainPattern, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('gain', Array(np.zeros(n), iotype='out', shape=(n,)))
        
        self.add('azimuthGS', Array(np.zeros(n), iotype='in', shape=(self.n,)))
        self.add('elevationGS', Array(np.zeros(n), iotype='in', 
                                      shape=(self.n,)))
        
        rawGdata = np.genfromtxt('CADRE/data/Comm/Gain.txt')
        rawG = (10**(rawGdata/10.0)).reshape((361,361),order='F')
        
        pi = np.pi
        az = np.linspace(0,2*pi,361)
        el = np.linspace(0,2*pi,361)

        self.MBI = MBI.MBI(rawG,[az,el],[15,15],[4,4])
        self.x = np.zeros((self.n,2),order='F')

    def setx(self):
        self.x[:,0] = self.azimuthGS
        self.x[:,1] = self.elevationGS

    def linearize(self):
        result = self.lib.fixangles(self.n, self.azimuthGS, self.elevationGS)
        self.azimuthGS, self.elevationGS = result
        self.x[:,0] = self.azimuthGS
        self.x[:,1] = self.elevationGS
        self.dg_daz = self.MBI.evaluate(self.x,1)[:,0]
        self.dg_del = self.MBI.evaluate(self.x,2)[:,0]

    def execute(self):
        result = self.lib.fixangles(self.n, self.azimuthGS, self.elevationGS)
        self.x[:,0] = result[0]
        self.x[:,1] = result[1]
        self.gain = self.MBI.evaluate(self.x)[:,0]

    def applyDer(self, arg, result):
        result['gain'] = np.zeros(self.n)
        if 'azimuthGS' in arg:
            result['gain'] += self.dg_daz * arg['azimuthGS'][:]
        if 'elevationGS' in arg:
            result['gain'] += self.dg_del * arg['elevationGS'][:]
        return result

    def applyDerT(self, arg, result):
        if 'gain' in arg:
            result['azimuthGS'] = self.dg_daz * arg['gain'][:]
            result['elevationGS'] = self.dg_del * arg['gain'][:]
        return result


class Comm_GSposEarth(Component):
    
    lon = Float(0, iotype="in")
    lat = Float(0, iotype="in")
    alt = Float(0, iotype="in")
    
    def __init__(self, n):
        super(Comm_GSposEarth, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        self.add('r_e2g_E', Array(np.zeros((3, self.n)), iotype='out', 
                                  shape=(3, self.n)))

    def linearize(self):
        result = self.lib.computejacobiangs(self.lon, self.lat, self.alt)
        self.dr_dlon, self.dr_dlat, self.dr_dalt = result

    def execute(self):
        self.r_e2g_E = self.lib.computegs(self.n, self.lon, self.lat, self.alt)

    def applyDer(self, arg, result):
        result['r_e2g_E'] = np.zeros((3, self.n))
        if 'lon' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dlon[k] * arg['lon']
        if 'lat' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dlat[k] * arg['lat']
        if 'alt' in arg:
            for k in xrange(3):
                result['r_e2g_E'][k,:] += self.dr_dalt[k] * arg['alt']
        return result

    def applyDerT(self, arg, result):
        if 'r_e2g_E' in arg:
            result['lon'] = 0.
            result['lat'] = 0.
            result['alt'] = 0.
            for k in xrange(3):
                result['lon'] += self.dr_dlon[k] * np.sum(arg['r_e2g_E'][k,:])
                result['lat'] += self.dr_dlat[k] * np.sum(arg['r_e2g_E'][k,:])
                result['alt'] += self.dr_dalt[k] * np.sum(arg['r_e2g_E'][k,:])
        return result


class Comm_GSposECI(Component):

    def __init__(self, n):
        super(Comm_GSposECI, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib

        self.add('r_e2g_I', Array(np.zeros((3, self.n)), iotype='out', 
                                  shape=(3, self.n)))
        
        self.add('O_IE', Array(np.zeros((3, 3, self.n)), iotype='in', 
                               shape=(3, 3, self.n)))
        self.add('r_e2g_E', Array(np.zeros((3, self.n)), iotype='in', 
                                  shape=(3, self.n)))

    def linearize(self):
        result = self.lib.computepositionrotdjacobian(self.n, 
                                                      self.r_e2g_E, 
                                                      self.O_IE)
        self.J1, self.J2 = result

    def execute(self):
        self.r_e2g_I[:] = self.lib.computepositionrotd(self.n, 
                                                      self.r_e2g_E, 
                                                       self.O_IE)

    def applyDer(self, arg, result):
        if 'O_IE' in arg and 'r_e2g_E' in arg:
            result['r_e2g_I'] = np.zeros((3, self.n))
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['r_e2g_I'][k,:] += self.J1[:,k,u,v] * arg['O_IE'][u,v,:]
                for j in xrange(3):
                    result['r_e2g_I'][k,:] += self.J2[:,k,j] * arg['r_e2g_E'][j,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_e2g_I' in arg:
            result['O_IE'] = np.zeros((3, 3, self.n))
            result['r_e2g_E'] = np.zeros((3, self.n))
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['O_IE'][u,v,:] += self.J1[:,k,u,v] * arg['r_e2g_I'][k,:]
                for j in xrange(3):
                    result['r_e2g_E'][j,:] += self.J2[:,k,j] * arg['r_e2g_I'][k,:]
        return result


class Comm_LOS(Component):

    def __init__(self, n):
        super(Comm_LOS, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.CommLib').lib.CommLib
        
        self.add('CommLOS', Array(np.zeros(n), iotype='out', shape=(self.n,)))
        
        self.add('r_b2g_I', Array(np.zeros((3, n)), iotype='in', 
                                  shape=(3, self.n)))
        self.add('r_e2g_I', Array(np.zeros((3, n)), iotype='in', 
                                  shape=(3, self.n)))

    def linearize(self):
        result = self.lib.computejacobianlos(self.n, self.r_b2g_I, self.r_e2g_I)
        self.dLOS_drb, self.dLOS_dre = result

    def execute(self):
        self.CommLOS = self.lib.computelos(self.n, self.r_b2g_I, self.r_e2g_I)

    def applyDer(self, arg, result):
        if 'r_b2g_I' in arg:
            result['CommLOS'] = np.zeros(self.n)
            for k in xrange(3):
                result['CommLOS'] += self.dLOS_drb[:,k] * arg['r_b2g_I'][k,:]
                result['CommLOS'] += self.dLOS_dre[:,k] * arg['r_e2g_I'][k,:]
        return result
    
    def applyDerT(self, arg, result):
        if 'CommLOS' in arg:
            result['r_b2g_I'] = np.zeros((3,self.n))
            result['r_e2g_I'] = np.zeros((3,self.n))
            for k in xrange(3):
                result['r_b2g_I'][k,:] += self.dLOS_drb[:,k] * arg['CommLOS']
                result['r_e2g_I'][k,:] += self.dLOS_dre[:,k] * arg['CommLOS']
        return result


class Comm_VectorAnt(Component):

    def __init__(self, n):
        super(Comm_VectorAnt, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib

        self.add('r_b2g_A', Array(np.zeros((3, n)), iotype='out', shape=(3, n)))
        
        self.add('r_b2g_B', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        self.add('O_AB', Array(np.zeros((3, 3, n)), iotype='in', 
                               shape=(3, 3, n)))

    def linearize(self):
        result = self.lib.computepositionrotdjacobian(self.n, self.r_b2g_B, 
                                                      self.O_AB)
        self.J1, self.J2 = result

    def execute(self):
        self.r_b2g_A = self.lib.computepositionrotd(self.n, self.r_b2g_B, 
                                                        self.O_AB)

    def applyDer(self, arg, result):
        if 'O_AB' in arg and 'r_b2g_B' in arg:
            result['r_b2g_A'] = np.zeros((3, self.n))
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['r_b2g_A'][k,:] += self.J1[:,k,u,v] * arg['O_AB'][u,v,:]
                for j in xrange(3):
                    result['r_b2g_A'][k,:] += self.J2[:,k,j] * arg['r_b2g_B'][j,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_b2g_A' in arg:
            result['r_b2g_B'] = np.zeros((3, self.n))
            result['O_AB'] = np.zeros((3, 3, self.n))
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['O_AB'][u,v,:] += self.J1[:,k,u,v] * arg['r_b2g_A'][k,:]
                for j in xrange(3):
                    result['r_b2g_B'][j,:] += self.J2[:,k,j] * arg['r_b2g_A'][k,:]
        return result


class Comm_VectorBody(Component):

    def __init__(self, n):
        super(Comm_VectorBody, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('r_b2g_B', Array(np.zeros((3, n)), iotype='out', shape=(3, n)))
        
        self.add('r_b2g_I', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        self.add('O_BI', Array(np.zeros((3, 3, n)), iotype='in', 
                               shape=(3, 3, n)))

    def linearize(self):
        result = self.lib.computepositionrotdjacobian(self.n, self.r_b2g_I, 
                                                      self.O_BI)
        self.J1, self.J2 = result

    def execute(self):
        self.r_b2g_B = self.lib.computepositionrotd(self.n, self.r_b2g_I, 
                                                       self.O_BI)

    def applyDer(self, arg, result):
        if 'O_BI' in arg and 'r_b2g_I' in arg:
            result['r_b2g_B'] = np.zeros((3, self.n))
            for k in range(3):
                for u in range(3):
                    for v in range(3):
                        result['r_b2g_B'][k,:] += self.J1[:,k,u,v] * arg['O_BI'][u,v,:]
                for j in range(3):
                    result['r_b2g_B'][k,:] += self.J2[:,k,j] * arg['r_b2g_I'][j,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_b2g_B' in arg:
            result['O_BI'] = np.zeros((3, 3, self.n))
            result['r_b2g_I'] = np.zeros((3, self.n))
            for k in range(3):
                for u in range(3):
                    for v in range(3):
                        result['O_BI'][u,v,:] += self.J1[:,k,u,v] * arg['r_b2g_B'][k,:]
                for j in range(3):
                    result['r_b2g_I'][j,:] += self.J2[:,k,j] * arg['r_b2g_B'][k,:]
        return result
    

class Comm_VectorECI(Component):

    def __init__(self, n):
        super(Comm_VectorECI, self).__init__()
        self.n = n
        
        self.add('r_b2g_I', Array(np.zeros((3, n)), iotype='out', shape=(3, n)))
        
        self.add('r_e2g_I', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        self.add('r_e2b_I', Array(np.zeros((6, n)), iotype='in', shape=(6, n)))

    def linearize(self):
        return

    def execute(self):
        self.r_b2g_I = self.r_e2g_I - self.r_e2b_I[:3,:]

    def applyDer(self, arg, result):
        if 'r_e2g_I' in arg and 'r_e2b_I' in arg:
            result['r_b2g_I'] = arg['r_e2g_I']
            result['r_b2g_I'] -= arg['r_e2b_I'][:3,:]
        return result

    def applyDerT(self, arg, result):
        if 'r_b2g_I' in arg:
            result['r_e2b_I'] = np.zeros((6, self.n))
            result['r_e2g_I'] = arg['r_b2g_I']
            result['r_e2b_I'][:3,:] -= arg['r_b2g_I']
        return result


class Comm_VectorSpherical(Component):

    def __init__(self, n):
        super(Comm_VectorSpherical, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib
        
        self.add('azimuthGS', Array(np.zeros(n), iotype='out', shape=(n,)))
        self.add('elevationGS', Array(np.zeros(n), iotype='out', shape=(n,)))
        
        self.add('r_b2g_A', Array(np.zeros((3, n)), iotype='in', 
                                  shape=(3, self.n)))

    def linearize(self):
        result = self.lib.computepositionsphericaljacobian(self.n, 3*self.n, 
                                                           self.r_b2g_A)
        self.Ja1, self.Ji1, self.Jj1, self.Ja2, self.Ji2, self.Jj2  = result
        self.J1 = scipy.sparse.csc_matrix((self.Ja1,(self.Ji1,self.Jj1)), 
                                          shape=(self.n,3*self.n))
        self.J2 = scipy.sparse.csc_matrix((self.Ja2,(self.Ji2,self.Jj2)), 
                                          shape=(self.n,3*self.n))
        self.J1T = self.J1.transpose()
        self.J2T = self.J2.transpose()

    def execute(self):
        azimuthGS, elevationGS = self.lib.computepositionspherical(self.n, 
                                                                   self.r_b2g_A)
        self.azimuthGS = azimuthGS
        self.elevationGS = elevationGS

    def applyDer(self, arg, result):
        if 'r_b2g_A' in arg:
            r_b2g_A = arg['r_b2g_A'].reshape((3*self.n),order='F')
            result['azimuthGS'] = self.J1.dot(r_b2g_A)
            result['elevationGS'] = self.J2.dot(r_b2g_A)
        return result

    def applyDerT(self, arg, result):
        if 'azimuthGS' in arg and 'elevationGS' in arg:
            azimuthGS = arg['azimuthGS']
            elevationGS = arg['elevationGS']
            result['r_b2g_A'] = (self.J1T.dot(azimuthGS) + 
                                    self.J2T.dot(elevationGS)).reshape((3, self.n),order='F')
        return result
