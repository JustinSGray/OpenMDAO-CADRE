from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component, Assembly
import numpy as np


class Attitude_Angular(Component):

    def __init__(self, n=2):
        super(Attitude_Angular, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.AttitudeLib').lib.AttitudeLib

        self.add('w_B', Array(np.zeros((3, n)), iotype='out', shape=(3,n)))
        
        self.add('O_BI', Array(np.zeros((3, 3, n)), iotype='in', 
                               shape=(3, 3, n)))
        self.add('Odot_BI', Array(np.zeros((3, 3, n)), iotype='in', 
                                  shape=(3, 3, n)))

    def linearize(self):
        #self.dw_dO, self.dw_dOdot = self.lib.computejacobianw(self.n, self.O_BI,self.Odot_BI)

        for i in range(0,self.n):
            self.dw_dOdot[i,0,2,:] = self.O_BI[1,:,i]
            self.dw_dO[i,0,1,:] = self.Odot_BI[2,:,i]
                    
            self.dw_dOdot[i,1,0,:] = self.O_BI[2,:,i]
            self.dw_dO[i,1,2,:] = self.Odot_BI[0,:,i]
                    
            self.dw_dOdot[i,2,1,:] = self.O_BI[0,:,i]
            self.dw_dO[i,2,0,:] = self.Odot_BI[1,:,i]


    def execute(self):
        #self.w_B = self.lib.computew(self.n, self.O_BI, self.Odot_BI)
        for i in range(0,self.n):
            self.w_B[0,i] = np.dot(self.Odot_BI[2,:,i] , self.O_BI[1,:,i])
            self.w_B[1,i] = np.dot(self.Odot_BI[0,:,i] , self.O_BI[2,:,i])
            self.w_B[2,i] = np.dot(self.Odot_BI[1,:,i] , self.O_BI[0,:,i])

    def applyDer(self, arg, result):
        if 'O_BI' in arg and 'Odot_BI' in arg:
            result['w_B'] = np.zeros((3, self.n))
            for k in xrange(3):
                for i in xrange(3):
                    for j in xrange(3):
                        result['w_B'][k,:] += self.dw_dO[:,k,i,j] * arg['O_BI'][i,j,:]
                        result['w_B'][k,:] += self.dw_dOdot[:,k,i,j] * arg['Odot_BI'][i,j,:]
        return result

    def applyDerT(self, arg, result):
        if 'w_B' in arg:
            result['O_BI'] = np.zeros((3, 3, self.n))
            result['Odot_BI'] = np.zeros((3, 3, self.n))
            for k in xrange(3):
                for i in xrange(3):
                    for j in xrange(3):
                        result['O_BI'][i,j,:] += self.dw_dO[:,k,i,j] * arg['w_B'][k,:]
                        result['Odot_BI'][i,j,:] += self.dw_dOdot[:,k,i,j] * arg['w_B'][k,:]
        return result


class Attitude_AngularRates(Component):

    def __init__(self, n=2, h = 0.01):
        super(Attitude_AngularRates, self).__init__()
        self.n = n
        self.add('h', Float(28.8, iotype='in'))

        self.add('wdot_B', Array(np.zeros((3, n)), iotype='out', shape=(3, n)))
        
        self.add('w_B', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        
    def linearize(self):
        return
        
    def execute(self):
        for k in xrange(3):
            self.wdot_B[k,0] += self.w_B[k,1] / self.h
            self.wdot_B[k,0] -= self.w_B[k,0] / self.h
            self.wdot_B[k,1:-1] += self.w_B[k,2:] / 2.0 / self.h
            self.wdot_B[k,1:-1] -= self.w_B[k,:-2] / 2.0 / self.h
            self.wdot_B[k,-1] += self.w_B[k,-1] / self.h
            self.wdot_B[k,-1] -= self.w_B[k,-2] / self.h

    def applyDer(self, arg, result):
        if 'w_B' in arg:
            result['wdot_B'] = np.zeros((3, self.n))
            for k in xrange(3):
                result['wdot_B'][k,0] += arg['w_B'][k,1] / self.h
                result['wdot_B'][k,0] -= arg['w_B'][k,0] / self.h
                result['wdot_B'][k,1:-1] += arg['w_B'][k,2:] / 2.0 / self.h
                result['wdot_B'][k,1:-1] -= arg['w_B'][k,:-2] / 2.0 / self.h
                result['wdot_B'][k,-1] += arg['w_B'][k,-1] / self.h
                result['wdot_B'][k,-1] -= arg['w_B'][k,-2] / self.h
        return result

    def applyDerT(self, arg, result):
        if 'wdot_B' in arg:
            result['w_B'] = np.zeros((3, self.n))
            for k in xrange(3):
                result['w_B'][k,1] += arg['wdot_B'][k,0] / self.h
                result['w_B'][k,0] -= arg['wdot_B'][k,0] / self.h
                result['w_B'][k,2:] += arg['wdot_B'][k,1:-1] / 2.0 / self.h
                result['w_B'][k,:-2] -= arg['wdot_B'][k,1:-1] / 2.0 / self.h
                result['w_B'][k,-1] += arg['wdot_B'][k,-1] / self.h
                result['w_B'][k,-2] -= arg['wdot_B'][k,-1] / self.h
        return result


class Attitude_Attitude(Component):

    def __init__(self, n=2):
        super(Attitude_Attitude, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.AttitudeLib').lib.AttitudeLib

        self.add('O_RI', Array(np.zeros((3, 3, n)), iotype='out', 
                               shape=(3, 3, n)))

        self.add('r_e2b_I', Array(np.zeros((6, n)), iotype='in', shape=(6, n)))

    def linearize(self):
        #self.dO_dr = self.lib.computejacobianori(self.n, self.r_e2b_I)
        for i in range(0,self.n):
            r = r_e2b_I[0:3,i]
            v = r_e2b_I[3:6,i]
            
            normr = np.sqrt(np.dot(r,r))
            normv = np.sqrt(np.dot(v,v))
            if normr < 1e-10:
                normr = 1e-10
            if normv < 1e-10:
                normv = 1e-10
            
            r = r / normr
            v = v / normv

            for k in range(0,3):
                dr_dr[k,k] = dr_dr[k,k] + 1.0 / normr
                dv_dv[k,k] = dv_dv[k,k] + 1.0 / normv
                dr_dr[:,k] = dr_dr[:,k] - r_e2b_I[0:3,i] / normr**2 * r_e2b_I[k,i] / normr
                dv_dv[:,k] = dv_dv[:,k] - r_e2b_I[3:6,i] / normv**2 * r_e2b_I[2+k,i]/ normv

            vx[0,:] = (0., -v[2], v[1])
            vx[1,:] = (v[2], 0., -v[0])
            vx[2,:] = (-v[1], v[0], 0.)
    
            dvx_dv[0,:,0] = (0., 0., 0.)
            dvx_dv[1,:,0] = (0., 0., -1.)
            dvx_dv[2,:,0] = (0., 1., 0.)
    
            dvx_dv[0,:,1] = (0., 0., 1.)
            dvx_dv[1,:,1] = (0., 0., 0.)
            dvx_dv[2,:,1] = (-1., 0., 0.)
    
            dvx_dv[0,:,2] = (0., -1., 0.)
            dvx_dv[1,:,2] = (1., 0., 0.)
            dvx_dv[2,:,2] = (0., 0., 0.)
    
            iB =  np.matmul(vx, r)
            jB = -np.matmul(vx, iB)
    
            diB_dr = vx
            diB_dv[:,0] = np.matmul(dvx_dv[:,:,0],r)
            diB_dv[:,1] = np.matmul(dvx_dv[:,:,1],r)
            diB_dv[:,3] = np.matmul(dvx_dv[:,:,2],r)
    
            djB_diB = -vx
            djB_dv[:,0] = -np.matmul(dvx_dv[:,:,0],iB)
            djB_dv[:,1] = -np.matmul(dvx_dv[:,:,1],iB)
            djB_dv[:,2] = -np.matmul(dvx_dv[:,:,2],iB)
            
            dO_dr[i,0,:,0:3] = np.matmul(diB_dr , dr_dr)
            dO_dr[i,0,:,3:6] = np.matmul(diB_dv , dv_dv)
            
            dO_dr[i,1,:,0:3] = np.matmul(matmul(djB_diB, diB_dr) , dr_dr)
            dO_dr[i,1,:,3:6] = np.matmul(matmul(djB_diB, diB_dv) + djB_dv , dv_dv)
            
            dO_dr[i,2,:,3:6] = -dv_dv


    def execute(self):
        #self.O_RI = self.lib.computeori(self.n, self.r_e2b_I)
        for i in range(0,self.n):
            r = self.r_e2b_I[0:3,i]
            v = self.r_e2b_I[3:6,i]
            
            normr = np.sqrt(np.dot(r,r))
            normv = np.sqrt(np.dot(v,v))
            if normr < 1e-10:
                normr = 1e-10
            if normv < 1e-10:
                normv = 1e-10
            
            r = r / normr
            v = v / normv

            vx = np.zeros((3,3))
            vx[0,:] = (0., -v[2], v[1])
            vx[1,:] = (v[2], 0., -v[0])
            vx[2,:] = (-v[1], v[0], 0.)
            
            iB =  np.dot(vx,r)
            jB = -np.dot(vx,iB)

            self.O_RI[0,:,i] = iB
            self.O_RI[1,:,i] = jB
            self.O_RI[2,:,i] = -v

    def applyDer(self, arg, result):
        if 'r_e2b_I' in arg:
            result['O_RI'] = np.zeros((3, 3, self.n))
            for k in xrange(3):
                for j in xrange(3):
                    for i in xrange(6):
                        result['O_RI'][k,j,:] += self.dO_dr[:,k,j,i] * arg['r_e2b_I'][i,:]
        return result

    def applyDerT(self, arg, result):
        if 'O_RI' in arg:
            result['r_e2b_I'] = np.zeros((6, self.n))
            for k in xrange(3):
                for j in xrange(3):
                    for i in xrange(6):
                        result['r_e2b_I'][i,:] += self.dO_dr[:,k,j,i] * arg['O_RI'][k,j,:]
        return result


class Attitude_Roll(Component):

    def __init__(self, n=2):
        super(Attitude_Roll, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.AttitudeLib').lib.AttitudeLib

        self.add('O_BR', Array(np.zeros((3, 3, n)), iotype='out', 
                               shape=(3, 3, n)))

        self.add('Gamma', Array(np.zeros(n), iotype='in', shape=(n,)))

    def linearize(self):
        #self.dO_dg = self.lib.computejacobianobr(self.n, self.Gamma)
        for i in range(0,self.n):
            self.dO_dg[i,0,:] = (-np.sin(Gamma[i]), np.cos(Gamma[i]), 0.)
            self.dO_dg[i,1,:] = (-np.cos(Gamma[i]), -np.sin(Gamma[i]), 0.)
            self.dO_dg[i,2,:] = (0., 0., 0.)

    def execute(self):
        #self.O_BR = self.lib.computeobr(self.n, self.Gamma)
        for i in range(0,self.n):
            self.O_BR[0,:,i] = (np.cos(self.Gamma[i]), np.sin(self.Gamma[i]), 0.)
            self.O_BR[1,:,i] = (-np.sin(self.Gamma[i]), np.cos(self.Gamma[i]), 0.)
            self.O_BR[2,:,i] = (0., 0., 1.)

    def applyDer(self, arg, result):
        if 'Gamma' in arg:
            result['O_BR'] = np.zeros((3, 3, self.n))
            for k in xrange(3):
                for j in xrange(3):
                    result['O_BR'][k,j,:] += self.dO_dg[:,k,j] * arg['Gamma']
        return result

    def applyDerT(self, arg, result):
        if 'O_BR' in arg:
            result['Gamma'] = np.zeros(self.n)
            for k in xrange(3):
                for j in xrange(3):
                    result['Gamma'] += self.dO_dg[:,k,j] * arg['O_BR'][k,j,:]
        return result
    
    
class Attitude_RotationMtx(Component):

    def __init__(self, n=2):
        super(Attitude_RotationMtx, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.AttitudeLib').lib.AttitudeLib

        self.add('O_BI', Array(np.zeros((3, 3, n)), iotype='out', 
                               shape=(3, 3, n)))
        
        self.add('O_BR', Array(np.zeros((3, 3, n)), iotype='in', 
                               shape=(3, 3, n)))
        self.add('O_RI', Array(np.zeros((3, 3, n)), iotype='in', 
                               shape=(3, 3, n)))
    
    def linearize(self):
        return
    
    def execute(self):
        #self.O_BI = self.lib.computeobi(self.n, self.O_BR, self.O_RI)
        for i in range(0,self.n):
            self.O_BI[:,:,i] = np.dot(self.O_BR[:,:,i], self.O_RI[:,:,i])

    def applyDer(self, arg, result):
        if 'O_RI' in arg and 'O_BR' in arg:
            result['O_BI'] = np.zeros((3, 3, self.n))
            for u in xrange(3):
                for v in xrange(3):
                    for k in xrange(3):
                        result['O_BI'][u,v,:] += self.O_BR[u,k,:] * arg['O_RI'][k,v,:]
                        result['O_BI'][u,v,:] += arg['O_BR'][u,k,:] * self.O_RI[k,v,:]
        return result

    def applyDerT(self, arg, result):
        if 'O_BI' in arg:
            result['O_RI'] = np.zeros((3, 3, self.n))
            result['O_BR'] = np.zeros((3, 3, self.n))
            for u in xrange(3):
                for v in xrange(3):
                    for k in xrange(3):
                        result['O_RI'][k,v,:] += self.O_BR[u,k,:] * arg['O_BI'][u,v,:]
                        result['O_BR'][u,k,:] += arg['O_BI'][u,v,:] * self.O_RI[k,v,:]
        return result


class Attitude_RotationMtxRates(Component):

    def __init__(self, n=1500):
        super(Attitude_RotationMtxRates, self).__init__()
        self.n = n
        self.add('h', Float(28.8, iotype='in'))
        self.add('Odot_BI', Array(np.zeros((3, 3, n)), iotype='out', 
                                  shape=(3, 3, n)))
        
        self.add('O_BI', Array(np.zeros((3, 3, n)), iotype='in', 
                               shape=(3, 3, n)))


    def linearize(self):
        return

    def execute(self):
        self.Odot_BI = np.zeros((3, 3, self.n))
        for k in range(3):
            for j in range(3):
                self.Odot_BI[k,j,0] += self.O_BI[k,j,1] / self.h
                self.Odot_BI[k,j,0] -= self.O_BI[k,j,0] / self.h
                self.Odot_BI[k,j,1:-1] += self.O_BI[k,j,2:] / 2.0 / self.h
                self.Odot_BI[k,j,1:-1] -= self.O_BI[k,j,:-2] / 2.0 / self.h
                self.Odot_BI[k,j,-1] += self.O_BI[k,j,-1] / self.h
                self.Odot_BI[k,j,-1] -= self.O_BI[k,j,-2] / self.h

    def applyDer(self, arg, result):
        if 'O_BI' in arg:
            result['Odot_BI'] = np.zeros((3, 3, self.n))
            for k in xrange(3):
                for j in xrange(3):
                    result['Odot_BI'][k,j,0] += arg['O_BI'][k,j,1] / self.h
                    result['Odot_BI'][k,j,0] -= arg['O_BI'][k,j,0] / self.h
                    result['Odot_BI'][k,j,1:-1] += arg['O_BI'][k,j,2:] / 2.0 / self.h
                    result['Odot_BI'][k,j,1:-1] -= arg['O_BI'][k,j,:-2] / 2.0 / self.h
                    result['Odot_BI'][k,j,-1] += arg['O_BI'][k,j,-1] / self.h
                    result['Odot_BI'][k,j,-1] -= arg['O_BI'][k,j,-2] / self.h
        return result

    def applyDerT(self, arg, result):
        if 'Odot_BI' in arg:
            result['O_BI'] = np.zeros((3, 3, self.n))
            for k in xrange(3):
                for j in xrange(3):
                    result['O_BI'][k,j,1] += arg['Odot_BI'][k,j,0] / self.h
                    result['O_BI'][k,j,0] -= arg['Odot_BI'][k,j,0] / self.h
                    result['O_BI'][k,j,2:] += arg['Odot_BI'][k,j,1:-1] / 2.0 / self.h
                    result['O_BI'][k,j,:-2] -= arg['Odot_BI'][k,j,1:-1] / 2.0 / self.h
                    result['O_BI'][k,j,-1] += arg['Odot_BI'][k,j,-1] / self.h
                    result['O_BI'][k,j,-2] -= arg['Odot_BI'][k,j,-1] / self.h
        return result


class Attitude_Sideslip(Component):

    def __init__(self, n=2):
        super(Attitude_Sideslip, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.KinematicsLib').lib.KinematicsLib

        self.add('v_e2b_B', Array(np.zeros((3, n)), iotype='out', shape=(3,n)))
        
        self.add('r_e2b_I', Array(np.zeros((6, n)), iotype='in', shape=(6, n)))
        self.add('O_BI', Array(np.zeros((3, 3, n)), iotype='in', shape=(3, 3, n)))        

    def linearize(self):
        #self.J1, self.J2 = self.lib.computepositionrotdjacobian(self.n,self.r_e2b_I[3:,:], self.O_BI)
        eye = np.eye(3, dtype=double)
        r = self.r_e2b_I[3:,:]
        for i in range(0,n):
            for k in range(0,3):
                for u in range(0,3):
                    for v in range(0,3):
                        self.J1[i,k,u,v] = eye[k,u]*r[v,i]
                for j in range(0,3):
                    self.J2[i,k,j] = self.O_BI[k,j,i]

    def execute(self):
        #self.v_e2b_B = self.lib.computepositionrotd(self.n, self.r_e2b_I[3:,:],self.O_BI)
        for i in range(0,self.n):
            self.v_e2b_B[:,i] = np.dot(self.O_BI[:,:,i],self.r_e2b_I[3:,i])

    def applyDer(self, arg, result):
        if 'O_BI' in arg and 'r_e2b_I' in arg:
            result['v_e2b_B'] = np.zeros((3, self.n))
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['v_e2b_B'][k,:] += self.J1[:,k,u,v] * arg['O_BI'][u,v,:]
                for j in xrange(3):
                    result['v_e2b_B'][k,:] += self.J2[:,k,j] * arg['r_e2b_I'][3+j,:]
        return result

    def applyDerT(self, arg, result):
        if 'v_e2b_B' in arg:
            result['O_BI'] = np.zeros((3, 3, self.n))
            result['r_e2b_I'] = np.zeros((6, self.n))
            for k in xrange(3):
                for u in xrange(3):
                    for v in xrange(3):
                        result['O_BI'][u,v,:] += self.J1[:,k,u,v] * arg['v_e2b_B'][k,:]
                for j in xrange(3):
                    result['r_e2b_I'][3+j,:] += self.J2[:,k,j] * arg['v_e2b_B'][k,:]
        return result


class Attitude_Torque(Component):

    J = np.zeros((3,3))
    J[0,:] = (0.018, 0., 0.)
    J[1,:] = (0., 0.018, 0.)
    J[2,:] = (0., 0., 0.006)

    def __init__(self, n=2):
        super(Attitude_Torque, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.AttitudeLib').lib.AttitudeLib

        self.add('T_tot', Array(np.zeros((3, n)), iotype='out', shape=(3,n)))

        self.add('w_B', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        self.add('wdot_B', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        
    def linearize(self):
        #self.dT_dw, self.dT_dwdot = self.lib.computejacobiant(self.n, self.w_B, self.wdot_B)
        dwx_dw = np.zeros((3,3,3))
        
        dwx_dw[0,:,0] = (0., 0., 0.)
        dwx_dw[1,:,0] = (0., 0., -1.)
        dwx_dw[2,:,0] = (0., 1., 0.)
        
        dwx_dw[0,:,1] = (0., 0., 1.)
        dwx_dw[1,:,1] = (0., 0., 0.)
        dwx_dw[2,:,1] = (-1., 0, 0.)
        
        dwx_dw[0,:,2] = (0., -1., 0)
        dwx_dw[1,:,2] = (1., 0., 0.)
        dwx_dw[2,:,2] = (0., 0., 0.)

        for i in range(0,self.n):
            wx[0,:] = (0., -w_B[2,i], w_B[1,i])
            wx[1,:] = (w_B[2,i], 0., -w_B[0,i])
            wx[2,:] = (-w_B[1,i], w_B[0,i], 0.)
            
            dT_dwdot[i,:,:] = self.J
            dT_dw[i,:,:] = np.dot(wx, self.J)
            for k in range(0,3):
                dT_dw[i,:,k] = dT_dw[i,:,k] + np.dot(dwx_dw[:,:,k], np.dot(self.J,w_B[:,i]))


    def execute(self):
        #self.T_tot = self.lib.computet(self.n, self.w_B, self.wdot_B)
        wx = np.zeros((3,3))
        for i in range(0,self.n):
            wx[0,:] = (0., -self.w_B[2,i], self.w_B[1,i])
            wx[1,:] = (self.w_B[2,i], 0., -self.w_B[0,i])
            wx[2,:] = (-self.w_B[1,i], self.w_B[0,i], 0.)
            self.T_tot[:,i] = np.dot(self.J,self.wdot_B[:,i]) + np.dot(wx, np.dot(self.J,self.w_B[:,i]))

    def applyDer(self, arg, result):
        if 'w_B' in arg and 'wdot_B' in arg:
            result['T_tot'] = np.zeros((3, self.n))
            for k in xrange(3):
                for j in xrange(3):
                    result['T_tot'][k,:] += self.dT_dw[:,k,j] * arg['w_B'][j,:]
                    result['T_tot'][k,:] += self.dT_dwdot[:,k,j] * arg['wdot_B'][j,:]
        return result

    def applyDerT(self, arg, result):
        if 'T_tot' in arg:
            result['w_B'] = np.zeros((3, self.n))
            result['wdot_B'] = np.zeros((3, self.n))
            for k in xrange(3):
                for j in xrange(3):
                    result['w_B'][j,:] += self.dT_dw[:,k,j] * arg['T_tot'][k,:]
                    result['wdot_B'][j,:] += self.dT_dwdot[:,k,j] * arg['T_tot'][k,:]
        return result

