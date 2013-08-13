""" RK4 time integration component """

import numpy as np
import scipy.sparse, scipy.sparse.linalg

from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Str


class RK4(Component):
    """Inherit from this component to use.
    
    State variable dimension: (num_states, num_time_points)
    External input dimension: (input width, num_time_points)
    """
    
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
    
    def initialize(self):
        """Set up dimensions and other data structures."""
        
        self.y = self.get(self.state_var)
        self.y0 = self.get(self.init_state_var)
        
        self.n_states, self.n = self.y.shape
        self.ny = self.n_states*self.n
        self.nJ = self.n_states*(self.n + self.n_states*(self.n-1))
    
        ext = []
        self.ext_index_map = {}
        for e in self.external_vars:
            var = self.get(e)
            self.ext_index_map[e] = len(ext)
            
            #TODO: Check that shape[-1]==self.n
            ext.extend(var.reshape(-1, self.n))
    
        
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
        e_vars = np.hstack((self.external_vars, self.fixed_external_vars))
        for i,var in enumerate(e_vars):
            self.reverse_name_map[var] = i
    
        self.name_map = dict([(v,k) for k,v in self.reverse_name_map.iteritems()])
        
        
        #TODO
        #  check that all ext arrays of of shape (self.n, )
        
        #TODO
        #check that length of state var and external
        # vars are the same length
    
    def f_dot(self, external, state):
        """time rate of change of state variables
            external: array or external variables for a single time step
            state: array of state variables for a single time step
            """
        raise NotImplementedError

    def df_dy(self, external, state):
        """derivatives of states with respect to states
            external: array or external variables for a single time step
            state: array of state variables for a single time step
            """
        
        raise NotImplementedError
    
    def df_dx(self, external, state):
        """derivatives of states with respect to external vars
            external: array or external variables for a single time step
            state: array of state variables for a single time step
            """
        raise NotImplementedError
    
    def execute(self):
        """Solve for the states."""
        
        self.initialize()
        
        self.y = self.y.reshape((self.ny, ))
        
        # Copy initial state into state array for t=0
        self.y[0:self.n_states] = self.y0[:]
        
        size = (self.n_states, self.n)
        self.a = np.zeros(size)
        self.b = np.zeros(size)
        self.c = np.zeros(size)
        self.d = np.zeros(size)
        
        for k in xrange(0, self.n-1):
            k1 = (k)*self.n_states
            k2 = (k+1)*self.n_states
            
            # Next state a function of current input
            ex = self.external[:, k] if self.external.shape[0] else np.array([])
            
            # Next state a function of current state
            y = self.y[k1:k2]
            
            a = self.a[:,k] = self.f_dot(ex, y)
            b = self.b[:,k] = self.f_dot(ex, y + self.h/2.*a)
            c = self.c[:,k] = self.f_dot(ex, y + self.h/2.*b)
            d = self.d[:,k] = self.f_dot(ex, y + self.h*c)
            
            self.y[self.n_states+k1:self.n_states+k2] = y + self.h/6.*(a + 2*b + 2*c + d)
        
        state_var_name = self.name_map['y']
        setattr(self, state_var_name, 
                self.y.T.reshape((self.n, self.n_states)).T)
    
        #print "executed", self.name
    
    def linearize(self):
        """Linearize about current point."""
        
        n_state = self.n_states
        I = np.eye(n_state)
        
        # Sparse Jacobian with respect to states
        self.Ja = np.zeros((self.nJ, ))
        self.Ji = np.zeros((self.nJ, ))
        self.Jj = np.zeros((self.nJ, ))
        
        # Full Jacobian with respect to inputs
        self.Jx = np.zeros((self.n, self.n_external, self.n_states))
        
        self.Ja[:self.ny] = np.ones((self.ny, ))
        self.Ji[:self.ny] = np.arange(self.ny)
        self.Jj[:self.ny] = np.arange(self.ny)
        
        for k in xrange(0, self.n-1):
            k1 = k*n_state
            k2 = k1 + n_state
            ex = self.external[:, k] if self.external.shape[0] else np.array([])
            y = self.y[k1:k2]
            
            a = self.a[:, k]
            b = self.b[:, k]
            c = self.c[:, k]
            d = self.d[:, k]
            
            # State vars
            df_dy = self.df_dy(ex, y)
            dg_dy = self.df_dy(ex, y + self.h/2.*a)
            dh_dy = self.df_dy(ex, y + self.h/2.*b)
            di_dy = self.df_dy(ex, y + self.h*c)
            
            da_dy = df_dy
            db_dy = dg_dy + dg_dy.dot(self.h/2.*da_dy)
            dc_dy = dh_dy + dh_dy.dot(self.h/2.*db_dy)
            dd_dy = di_dy + di_dy.dot(self.h*dc_dy)
            
            dR_dy = -I - self.h/6.*(da_dy + 2*db_dy + 2*dc_dy + dd_dy)
            
            for i in xrange(n_state):
                for j in xrange(n_state):
                    iJ = self.ny + i + n_state*(j + n_state*k)
                    self.Ja[iJ] = dR_dy[i, j]
                    self.Ji[iJ] = (k+1)*n_state + i
                    self.Jj[iJ] = (k)*n_state + j
            
            # External vars (Inputs)
            df_dx = self.df_dx(ex, y)
            dg_dx = self.df_dx(ex, y + self.h/2.*a)
            dh_dx = self.df_dx(ex, y + self.h/2.*b)
            di_dx = self.df_dx(ex, y + self.h*c)
            
            da_dx = df_dx
            db_dx = dg_dx + dg_dy.dot(self.h/2*da_dx)
            dc_dx = dh_dx + dh_dy.dot(self.h/2*db_dx)
            dd_dx = di_dx + di_dy.dot(self.h*dc_dx)
            
            # Input-State Jacobian at each time point.
            # No Jacobian with respect to previous time points.
            self.Jx[k+1,:,:] = self.h/6*(da_dx + 2*db_dx + 2*dc_dx + dd_dx).T
        
        self.J = scipy.sparse.csc_matrix((self.Ja, (self.Ji, self.Jj)),
                                         shape=(self.ny, self.ny))
        self.JT = self.J.transpose()
        self.Minv = scipy.sparse.linalg.splu(self.J).solve
    
    
    def apply_deriv(self, arg, result):
        """Matrix-vector product between Jacobian and arg. Result placed in
        result.
        """
        #result = self._applyJint(arg, result)
        result_ext = self._applyJext(arg)
        
        svar = self.state_var
        if svar in result:
            result[svar] += result_ext
        else:
            result[svar] = result_ext

    
    def _applyJint(self, arg, result):
        """Apply derivatives with respect to state variables."""
        
        arg = dict([(self.reverse_name_map[k],v) for k,v in arg.iteritems()])
        res1 = dict([(self.reverse_name_map[k],v) for k,v in result.iteritems()])
        
        if "y" in arg:
            flat_y = arg['y'].reshape((self.n_states*self.n),order='F')
            result["y"] = self.J.dot(flat_y).reshape((self.n_states,self.n),order='F')
        
        res1 = dict([(self.name_map[k],v) for k,v in res1.iteritems()])
        return res1

    def _applyJext(self, arg):
        """Apply derivatives with respect to inputs"""
        
        #Jx --> (n_times, n_external, n_states)
        n_state = self.n_states
        n_time = self.n
        result = np.zeros((n_state, n_time))
        
        # Collapse incoming a*b*...*c*n down to (ab...c)*n
        for name in self.external_vars:
            if name in arg:
                var = self.get(name)
                shape = var.shape
                arg[name] = arg[name].reshape((np.prod(shape[:-1]), 
                                               shape[-1]))
            
        # Time-varying inputs
        for k in xrange(n_time):
            for name in self.external_vars:
                if name in arg:
                    i_ext = self.ext_index_map[name]
                    ext_length = np.prod(arg[name][:, 0].shape)
                    for j in xrange(k):
                        Jsub = self.Jx[j+1, i_ext:i_ext+ext_length, :]
                        result[:, k] += Jsub.T.dot(arg[name][:, j])

        # Time-invariant inputs
        for name in self.fixed_external_vars:
            if name in arg:
                ext_var = getattr(self, name)
                if len(ext_var) > 1:
                    arg[name] = arg[name].flatten()
                i_ext = self.ext_index_map[name]
                ext_length = np.prod(ext_var.shape)
                for k in xrange(self.n_states):
                    result[k, 1:] += self.Jx[1, i_ext:i_ext+ext_length, k].dot(arg[name])
        
        return result

    def apply_derivT(self, arg, result):
        
        r1 = self.applyJintT(arg, result)
        r2 = self._applyJextT(arg, result)
        
        for k,v in r2.iteritems():
            if k in r1 and r1[k] is not None:
                r1[k] += v
            else:
                r1[k] = v

    def applyJintT(self, arg, result):
        
        res1 = dict([(self.reverse_name_map[k],v) for k,v in result.iteritems()])
        
        if self.state_var in arg:
            flat_y = arg[self.state_var].flatten()
            res1['y'] = self.JT.dot(flat_y).reshape((self.n_states, self.n))
        
        res1 =  dict([(self.name_map[k],v) for k,v in res1.iteritems()])
        return res1

    def _applyJextT(self, arg, required_results):
        
        #Jx --> (n_times, n_external, n_states)
        n_state = self.n_states
        n_time = self.n
        result = {}
        
        if self.state_var in arg:
            
            argsv = arg[self.state_var]
            
            # Time-varying inputs
            for name in self.external_vars:
                if name in required_results:
                    ext_var = getattr(self, name)
                    i_ext = self.ext_index_map[name]
                    ext_length = np.prod(ext_var.shape)
                    result[name] = np.zeros((n_state, n_time))
                    for k in xrange(n_time):
                        for j in xrange(k+1, n_time):
                            Jsub = self.Jx[j, i_ext:i_ext+ext_length, :].T
                            result[name][:, k] += Jsub.dot(argsv[:, j])
            
            # Time-invariant inputs
            for name in self.fixed_external_vars:
                ext_var = getattr(self, name)
                i_ext = self.ext_index_map[name]
                ext_length = np.prod(ext_var.shape)
                result[name] = np.zeros((ext_length, ))
                for k in xrange(n_state): 
                    result[name] += self.Jx[1:, i_ext:i_ext+ext_length, k].T.dot(argsv[k,1:])
            
            
            #result[self.init_state_var] = -arg[self.state_var][:,0]
        for k, v in result.iteritems(): 
            ext_var = getattr(self, k)
            result[k] = v.reshape(ext_var.shape, order='F')
        
        return result





