from openmdao.lib.datatypes.api import Float, Dict, Array, List
from openmdao.main.api import Component, Assembly
import numpy as np

class Attitude(Assembly):
    """
    CADRE Attitude assembly
    """

    def __init__(self, n=2):
        super(Attitude, self).__init__()
        
        self.add("angular", Attitude_Angular(n=n))
        self.driver.workflow.add("angular")
        
        self.make_connections()
        
    def make_connections(self):
        """
        Collects the names of all input and output variables for all
        components within the assembly (drivers excluded). 
        
        Then establishes connections between
        any output variable and input variable that has the same name, so 
        long as the variable name does not exist as an output to more than
        a single component (so excludes default outputs)
        """
        inputs, outputs = {}, {}
        for compname in self.list_components():
            if compname == "driver":
                continue
            
            comp_inputs = self.get(compname).list_inputs()
            
            for input_name in comp_inputs:
                if input_name not in inputs:
                    inputs[input_name] = [compname]
                else:
                    inputs[input_name].append(compname)
                
            comp_outputs = self.get(compname).list_outputs()
            
            for output_name in comp_outputs:
                if output_name not in outputs:
                    outputs[output_name] = [compname]
                else:
                    outputs[output_name].append(compname)
            
        print
        for varname in outputs.keys():
            comps = outputs[varname]
            if len(comps) > 1:
                continue
            frompath = '.'.join([comps[0], varname])
            
            if varname in inputs:
                for compname in inputs[varname]:
                    topath = '.'.join([compname, varname])
                    self.connect(frompath, topath)
                    print "Connecting", frompath, "to", topath, "..."

            
            
            
        
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
        self.dw_dO, self.dw_dOdot = self.lib.computejacobianw(self.n, self.O_BI,
                                                              self.Odot_BI)

    def execute(self):
        self.w_B = self.lib.computew(self.n, self.O_BI, self.Odot_BI)

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
        self.h = h

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
        self.dO_dr = self.lib.computejacobianori(self.n, self.r_e2b_I)

    def execute(self):
        self.O_RI = self.lib.computeori(self.n, self.r_e2b_I)

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
        self.dO_dg = self.lib.computejacobianobr(self.n, self.Gamma)

    def execute(self):
        self.O_BR = self.lib.computeobr(self.n, self.Gamma)

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
        self.O_BI = self.lib.computeobi(self.n, self.O_BR, self.O_RI)

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

    def __init__(self, n=2, h=0.01):
        super(Attitude_RotationMtxRates, self).__init__()
        self.n = n
        self.h = h

        self.add('Odot_BI', Array(np.zeros((3, 3, n)), iotype='out', 
                                  shape=(3, 3, n)))
        
        self.add('O_BI', Array(np.zeros((3, 3, n)), iotype='in', 
                               shape=(3, 3, n)))

    def linearize(self):
        return

    def execute(self):
        for k in xrange(3):
            for j in xrange(3):
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
        self.J1, self.J2 = self.lib.computepositionrotdjacobian(self.n, 
                                                                self.r_e2b_I[3:,:], 
                                                                self.O_BI)

    def execute(self):
        self.v_e2b_B = self.lib.computepositionrotd(self.n, self.r_e2b_I[3:,:], 
                                                    self.O_BI)

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

    def __init__(self, n=2):
        super(Attitude_Torque, self).__init__()
        self.n = n
        self.lib = __import__('CADRE.lib.AttitudeLib').lib.AttitudeLib

        self.add('T_tot', Array(np.zeros((3, n)), iotype='out', shape=(3,n)))

        self.add('w_B', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        self.add('wdot_B', Array(np.zeros((3, n)), iotype='in', shape=(3, n)))
        
    def linearize(self):
        self.dT_dw, self.dT_dwdot = self.lib.computejacobiant(self.n, self.w_B, 
                                                              self.wdot_B)

    def execute(self):
        self.T_tot = self.lib.computet(self.n, self.w_B, self.wdot_B)

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

