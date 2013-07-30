from operator import mul

import numpy as np
import scipy.sparse, scipy.sparse.linalg

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Str


class RK4(Component):
    
    h = Float(.01, units="s", iotype="in",
              desc="time step used for the integration")
    
    state_var = Str("", iotype="in",
                    desc="name of the variable to be used for time integration")
    
    init_state_var = Str("", iotype="in",
                         desc="name of the variable to be used for initial conditions")
    
    external_vars = Array([], iotype="in", dtype=str,
                          desc="list of names of variables that are external to the system, but DO vary with time")
    
    fixed_external_vars = Array([], iotype="in", dtype=str,
                                desc="list of names of variables that are external to the system, but DO NOT vary with time")
    
    def __init__(self):
        super(RK4, self).__init__()
        #self.n = n
        #self.n_states = n_states
        #self.time_step = time_step
        #self.n_external = n_external
        
        #self.ny = self.n_states*self.n
        #self.nJ = self.n_states*self.n+self.n_states**2*(self.n-1)
    
    
    def check_config(self):
        super(RK4, self).check_config()

        y = self.y = self.get(self.state_var)
        y0 = self.y0 = self.get(self.init_state_var)
        
        self.n_states = y.shape[0]
        self.n = y.shape[1]

        ext = []
        self.ext_index_map = {}
        for e in self.external_vars:
            var = self.get(e)
            self.ext_index_map[e] = len(ext)
            
            np.set_printoptions(threshold=np.nan)
            import pprint
            print "printing var for "
            pprint.pprint( var )
        
            #TODO: Check that shape[-1]==self.n
            ext.extend(var.reshape(-1,self.n))
    
        
        for e in self.fixed_external_vars:
            var = self.get(e)
            self.ext_index_map[e] = len(ext)
            
            flat_var = var.flatten()
            #create n copies of the var
            ext.extend(np.tile(flat_var,(self.n, 1)).T)
                
        self.external = np.array(ext)
                
        #TODO
        #check that len(y0) = self.n_states
        
        self.n_external = len(ext)
        self.reverse_name_map = {
            self.state_var:'y',
            self.init_state_var:'y0'
        }
        e_vars = np.hstack((self.external_vars,self.fixed_external_vars))
        for i,var in enumerate(e_vars):
            self.reverse_name_map[var] = i
    
        self.name_map = dict([(v,k) for k,v in self.reverse_name_map.iteritems()])
        
        
        #TODO
        #  check that all ext arrays of of shape (self.n, )
        
        #TODO
        #check that length of state var and external
        # vars are the same length
        
        self.ny = self.n_states*self.n
        self.nJ = self.n_states*self.n+self.n_states**2*(self.n-1)
    
        print "TESTING HERE"
    
    
    def f_dot(self,external,state):
        """time rate of change of state variables
            external: array or external variables for a single time step
            state: array of state variables for a single time step
            """
        raise NotImplementedError

    def df_dy(self,external,state):
        """derivatives of states with respect to states
            external: array or external variables for a single time step
            state: array of state variables for a single time step
            """
        
        raise NotImplementedError
    
    def df_dx(self,external,state):
        """derivatives of states with respect to external vars
            external: array or external variables for a single time step
            state: array of state variables for a single time step
            """
        raise NotImplementedError
    
    def execute(self):
        
        #self.check_config()
        
        self.y = self.y.reshape((self.ny, ))
        self.y[0:self.n_states] = self.y0[:]
        
        size = (self.n_states, self.n)
        self.a = np.zeros(size)
        self.b = np.zeros(size)
        self.c = np.zeros(size)
        self.d = np.zeros(size)
        
        for k in xrange(0,self.n-1):
            k1 = (k)*self.n_states
            k2 = (k+1)*self.n_states
            ex = self.external[:,k] if self.external.shape[0] else np.array([])
            y = self.y[k1:k2]
            
            a = self.a[:,k] = self.f_dot(ex,y)
            b = self.b[:,k] = self.f_dot(ex, y + self.h/2.*a)
            c = self.c[:,k] = self.f_dot(ex, y + self.h/2.*b)
            d = self.d[:,k] = self.f_dot(ex, y + self.h*c)
            
            self.y[self.n_states+k1:self.n_states+k2] = y + self.h/6.*(a + 2*b + 2*c + d)
        
        state_var_name = self.name_map['y']
        setattr(self, state_var_name, self.y.T.reshape((self.n,self.n_states)).T)
    
    
    def linearize(self):
        
        I = np.eye(self.n_states)
        self.Ja = np.zeros((self.nJ, ))
        self.Ji = np.zeros((self.nJ, ))
        self.Jj = np.zeros((self.nJ, ))
        self.Jx = np.zeros((self.n, self.n_external, self.n_states))
        
        self.Ja[:self.ny] = np.ones((self.ny, ))
        self.Ji[:self.ny] = np.arange(self.ny)
        self.Jj[:self.ny] = np.arange(self.ny)
        
        for k in xrange(0,self.n-1):
            k1 = (k)*self.n_states
            k2 = (k+1)*self.n_states
            ex = self.external[:,k] if self.external.shape[0] else np.array([])
            y = self.y[k1:k2]
            
            a = self.a[:,k]
            b = self.b[:,k]
            c = self.c[:,k]
            d = self.d[:,k]
            
            #state vars
            df_dy = self.df_dy(ex,y)
            dg_dy = self.df_dy(ex, y + self.h/2.*a)
            dh_dy = self.df_dy(ex, y + self.h/2.*b)
            di_dy = self.df_dy(ex, y + self.h*c)
            
            da_dy = df_dy
            db_dy = dg_dy + dg_dy.dot(self.h/2.*da_dy)
            dc_dy = dh_dy + dh_dy.dot(self.h/2.*db_dy)
            dd_dy = di_dy + di_dy.dot(self.h*dc_dy)
            
            dR_dy = - I - self.h/6.*(da_dy + 2*db_dy + 2*dc_dy + dd_dy)
            
            for i in xrange(self.n_states):
                for j in xrange(self.n_states):
                    iJ = self.ny + (k)*self.n_states**2 + (j)*self.n_states + i
                    self.Ja[iJ] = dR_dy[i,j]
                    self.Ji[iJ] = (k+1)*self.n_states + i
                    self.Jj[iJ] = (k)*self.n_states + j
            
            #external vars
            df_dx = self.df_dx(ex,y)
            dg_dx = self.df_dx(ex, y + self.h/2.*a)
            dh_dx = self.df_dx(ex, y + self.h/2.*b)
            di_dx = self.df_dx(ex, y + self.h*c)
            
            da_dx = df_dx
            db_dx = dg_dx + dg_dy.dot(self.h/2*da_dx)
            dc_dx = dh_dx + dh_dy.dot(self.h/2*db_dx)
            dd_dx = di_dx + di_dy.dot(self.h*dc_dx)
            
            self.Jx[k+1,:,:] = -self.h/6*(da_dx + 2*db_dx + 2*dc_dx + dd_dx).T
        
        self.J = scipy.sparse.csc_matrix((self.Ja,(self.Ji,self.Jj)),shape=(self.ny,self.ny))
        self.JT = self.J.transpose()
        self.Minv = scipy.sparse.linalg.splu(self.J).solve
    
    
    
    def applyJ(self, arg, result):
        
        r1 = self.applyJint(arg, result)
        #r2 = self.applyJext(arg, result)
        r2 = self._applyJext(arg, result)
        
        r3 = dict(r1)
        for k,v in r2.iteritems():
            if k in r3 and r3[k] is not None:
                r3[k] += v
            else:
                r3[k] = v
        return r3

    
    def applyJint(self, arg, result):
        arg = dict([(self.reverse_name_map[k],v) for k,v in arg.iteritems()])
        result = dict([(self.reverse_name_map[k],v) for k,v in result.iteritems()])
        
        if "y" in arg:
            flat_y = arg['y'].reshape((self.n_states*self.n),order='F')
            result["y"] = self.J.dot(flat_y).reshape((self.n_states,self.n),order='F')
        
        result =  dict([(self.name_map[k],v) for k,v in result.iteritems()])
        return result

    def applyJext(self, arg, result):
        raise NotImplementedError
    
    def _applyJext(self, arg, result):
        #Jx --> (n_times,n_external,n_states)
        state_var = getattr(self,self.state_var)
        result[self.state_var] = np.zeros(state_var.shape)
        
        for ext_var_name in self.external_vars:
            if ext_var_name in arg:
                ext_var = getattr(self,ext_var_name)
                i_ext = self.ext_index_map[ext_var_name]
                ext_length = np.prod(ext_var.shape)/self.n
                for k in xrange(self.n_states):
                    for j in xrange(ext_length):
                        result[self.state_var][k,1:] += self.Jx[1:,j+i_ext,k] * self.external[j+i_ext,:-1]
        
        for ext_var_name in self.fixed_external_vars:
            if ext_var_name in arg:
                ext_var = getattr(self,ext_var_name)
                i_ext = self.ext_index_map[ext_var_name]
                ext_length = np.prod(ext_var.shape)
                for k in xrange(self.n_states):
                    result[self.state_var][k,1:] += self.Jx[1,i_ext:i_ext+ext_length,k].dot(self.external[i_ext:i_ext+ext_length,0])
        if self.init_state_var in arg:
            result[self.state_var][:,0] -= arg[self.init_state_var][0]
        
        
        return result

    def applyJT(self, arg, result):
        
        r1 = self.applyJintT(arg, result)
        #r2 = self.applyJextT(arg, result)
        r2 = self._applyJextT(arg, result)
        
        
        r3 = dict(r1)
        for k,v in r2.iteritems():
            if k in r3 and r3[k] is not None:
                r3[k] += v
            else:
                r3[k] = v
        return r3

    def applyJintT(self, arg, result):
        
        arg = dict([(self.reverse_name_map[k],v) for k,v in arg.iteritems()])
        result = dict([(self.reverse_name_map[k],v) for k,v in result.iteritems()])
        
        if 'y' in arg:
            flat_y = arg['y'].flatten()
            result['y'] = self.JT.dot(flat_y).reshape((self.n_states,self.n))
        
        result =  dict([(self.name_map[k],v) for k,v in result.iteritems()])
        return result

    def applyJextT(self, arg, result):
        raise NotImplementedError
    
    def _applyJextT(self, arg, result):
        if self.state_var in arg:
            
            for ext_var_name in self.external_vars:
                ext_var = getattr(self,ext_var_name)
                i_ext = self.ext_index_map[ext_var_name]
                ext_length = np.prod(ext_var.shape)/self.n
                result[ext_var_name] = np.zeros((ext_length,self.n))
                for k in xrange(self.n_states):
                    for j in xrange(ext_length):
                        result[ext_var_name][j,:-1] += self.Jx[1:,j+i_ext,k] * arg[self.state_var][k,1:]
            
            
            for ext_var_name in self.fixed_external_vars:
                ext_var = getattr(self,ext_var_name)
                i_ext = self.ext_index_map[ext_var_name]
                ext_length = np.prod(ext_var.shape)
                result[ext_var_name] = np.zeros((ext_length, ))
                for k in xrange(self.n_states): 
                    result[ext_var_name] += self.Jx[1:,i_ext:i_ext+ext_length,k].T.dot(arg[self.state_var][k,1:])
            
            
            result[self.init_state_var] = -arg[self.state_var][:,0]
        
        for k,v in result.iteritems(): 
            ext_var = getattr(self,k)
            result[k] = v.reshape(ext_var.shape, order='F')
        
        return result





