import numpy as np
import scipy.sparse

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array

class Sun_LOS( Component ):

    def __init__(self, n=2): 
        super(Sun_LOS, self).__init__()

        self.n = n 

        self.r1 = 6378.137*0.85 # Earth's radius is 6378 km. 0.85 is the alpha in John Hwang's paper
        self.r2 = 6378.137

        self.add('r_e2b_I', Array(np.zeros((6, n), order='F'), size=(6,n, ), dtype=np.float, iotype="in"))
        self.add('r_e2s_I', Array(np.zeros((3, n), order='F'), size=(3,n, ), dtype=np.float, iotype="in"))

        self.add('LOS', Array(np.zeros((n, ), order='F'), size=(n, ), dtype=np.float, iotype="out", 
            desc="Line of Sight over Time"))

        self.lib = __import__('CADRE.lib.SunLib').lib.SunLib


    def execute(self): 
        for i in range( self.n ):
            r_b = self.r_e2b_I[:3,i]
            r_s = self.r_e2s_I[:3,i]
            dot = np.dot( r_b, r_s )
            cross = np.cross( r_b, r_s )
            dist = np.sqrt( cross.dot(cross) )
            
            if dot >= 0.0 :
                self.LOS[i] = 1.0
            elif dist <= self.r1 :
                self.LOS[i] = 0.0
            elif dist >= self.r2 :
                self.LOS[i] = 1.0
            else :
                x = ( dist - self.r1 ) / ( self.r2 - self.r1 )
                self.LOS[i] = 3 *x ** 2 - 2 * x**3

    def linearize(self): 

        Jab, Jib, Jjb, Jas, Jis, Jjs = \
                  self.lib.computelosjacobian(self.n, 3*self.n,
                                              self.r1, self.r2,
                                              self.r_e2b_I[:],
                                              self.r_e2s_I[:])
        self.Jb = scipy.sparse.csc_matrix((Jab,(Jib,Jjb)),shape=(self.n,6*self.n))
        self.Js = scipy.sparse.csc_matrix((Jas,(Jis,Jjs)),shape=(self.n,3*self.n))
        self.JbT = self.Jb.transpose()
        self.JsT = self.Js.transpose()

        return

        # Python version of Jacobian code - does not work yet

        # nj = self.n
        
        # Jab = np.zeros(shape=(nj,), dtype = np.float)
        # Jib = np.zeros(shape=(nj,), dtype = np.int)
        # Jjb = np.zeros(shape=(nj,), dtype = np.int)
        # Jas = np.zeros(shape=(nj,), dtype = np.float)
        # Jis = np.zeros(shape=(nj,), dtype = np.int)
        # Jjs = np.zeros(shape=(nj,), dtype = np.int)

        # # Working arrays
        # r_b = np.zeros(shape=(3,), dtype = np.int)
        # r_s = np.zeros(shape=(3,), dtype = np.int)
        # Bx = np.zeros(shape=(3,3,), dtype = np.int)
        # Sx = np.zeros(shape=(3,3,), dtype = np.int)
        # cross = np.zeros(shape=(3,), dtype = np.int)
        # ddist_cross = np.zeros(shape=(3,), dtype = np.int)
        # dcross_drb = np.zeros(shape=(3,3,), dtype = np.int)
        # dcross_drs = np.zeros(shape=(3,3,), dtype = np.int)
        # dLOS_dx = np.zeros(shape=(3,), dtype = np.int)
        # dLOS_drs = np.zeros(shape=(3,), dtype = np.int)
        # dLOS_drb = np.zeros(shape=(3,), dtype = np.int)

        # for i in range(self.n) :
        #     r_b = self.r_e2b_I[:3,i]
        #     r_s = self.r_e2s_I[:3,i]
        #     Bx = crossMatrix( r_b )
        #     Sx = crossMatrix(-  r_s )
        #     dot = np.dot( r_b, r_s )
        #     cross = np.cross( r_b, r_s )
        #     dist = np.sqrt( cross.dot(cross) )

        #     if dot >= 0.0 :
        #         dLOS_drb[:] = 0.0
        #         dLOS_drs[:] = 0.0
        #     elif dist <= self.r1 :
        #         dLOS_drb[:] = 0.0
        #         dLOS_drs[:] = 0.0
        #     elif dist >= self.r2 :
        #         dLOS_drb[:] = 0.0
        #         dLOS_drs[:] = 0.0
        #     else:
        #         x = (dist-self.r1)/(self.r2-self.r1)
        #         LOS = 3*x**2 - 2*x**3
        #         ddist_dcross = cross/dist
        #         dcross_drb = Sx
        #         dcross_drs = Bx
        #         dx_ddist = 1.0/(self.r2-self.r1)
        #         dLOS_dx = 6*x - 6*x**2
        #         dLOS_drb = dLOS_dx*dx_ddist*np.dot(ddist_dcross,dcross_drb)
        #         dLOS_drs = dLOS_dx*dx_ddist*np.dot(ddist_dcross,dcross_drs)

        #     for k in range(3) :
        #         iJ = i*3 + k
        #         Jab[iJ] = dLOS_drb[k]
        #         Jib[iJ] = i - 1
        #         Jjb[iJ] = (i-1)*6 + k - 1
        #         Jas[iJ] = dLOS_drs[k]
        #         Jis[iJ] = i - 1
        #         Jjs[iJ] = (i-1)*3 + k - 1

        # self.Jb = scipy.sparse.csc_matrix((Jab,(Jib,Jjb)),shape=(self.n,6*self.n))
        # self.Js = scipy.sparse.csc_matrix((Jas,(Jis,Jjs)),shape=(self.n,3*self.n))
        # self.JbT = self.Jb.transpose()
        # self.JsT = self.Js.transpose()

    def applyDer(self, arg, result):

        if not result['LOS'] :
            result['LOS'] = np.zeros( (self.n,) )
            
        if 'r_e2b_I' in arg: 
            r_e2b_I = arg['r_e2b_I'][:].reshape((6*self.n),order='F')
            result['LOS'] += self.Jb.dot(r_e2b_I) 

        if 'r_e2s_I' in arg: 
            r_e2s_I = arg['r_e2s_I'][:].reshape((3*self.n),order='F')
            result['LOS'] += self.Js.dot(r_e2s_I) 

        #     # from Justin battery code
        # if 'P_bat' in arg: 
        #     result['I_bat'] = self.dI_dP * arg['P_bat'] 
        # if 'temperature' in arg:  
        #     result['I_bat'] += self.dI_dT * arg['temperature'][4,:]
        # if 'SOC' in arg:     
        #     result['I_bat'] += self.dI_dSOC * arg['SOC'][0,:]

        #     # from CADRE 
        # r_e2b_I = arg('r_e2b_I')[:].reshape((6*self.n),order='F')
        # r_e2s_I = arg('r_e2s_I')[:].reshape((3*self.n),order='F')
        # res('LOS')[:] += self.Jb.dot(r_e2b_I) + self.Js.dot(r_e2s_I)

        return result   

    def applyDerT(self, arg, result): 

        if 'LOS' in arg: 
            LOS = arg['LOS']

            if not result['r_e2b_I'] :
                result['r_e2b_I'] = np.zeros( (6,self.n,) )
                
            if not result['r_e2s_I'] :
                result['r_e2s_I'] = np.zeros( (3,self.n,) )
                
            result['r_e2b_I'] += self.JbT.dot(LOS).reshape((6,self.n),order='F')
            result['r_e2s_I'] += self.JsT.dot(LOS).reshape((3,self.n),order='F')

        # # From Justin battery code
        # if 'I_bat' in arg: 
        #     result['P_bat'] = self.dI_dP * arg['I_bat']
        #     result['temperature'] = self.dI_dT * arg['I_bat']
        #     result['SOC'] = self.dI_dSOC * arg['I_bat']

        #     # from CADRE
        #             LOS = arg('LOS')[:]
        # res('r_e2b_I')[:] += self.JbT.dot(LOS).reshape((6,self.n),order='F')
        # res('r_e2s_I')[:] += self.JsT.dot(LOS).reshape((3,self.n),order='F')



        return result

def crossMatrix(v):

        # so m[1,0] is v[2], for example
        m = np.array( [
                         [ 0.0, -v[2], v[1] ],
                         [ v[2],  0.0, -v[0]],
                         [ -v[1], v[0], 0.0 ]
                       ] )
        
        # subroutine crossMatrix(v, M)
        
        #   implicit none
        
        #   !Input
        #   double precision, intent(in) ::  v(3)
        
        #   !Output
        #   double precision, intent(out) ::  M(3,3)
        
        #   M(1,:) = (/ dble(0), -v(3), v(2) /)
        #   M(2,:) = (/ v(3), dble(0), -v(1) /)
        #   M(3,:) = (/ -v(2), v(1), dble(0) /)
        
        # end subroutine crossMatrix


        return m

            
        
        
