from openmdao.main.api import Assembly
from openmdao.main.datatypes.api import Float, Array, Int
import numpy as np
from CADRE_assembly import CADRE
from pyopt_driver import pyopt_driver

class CADRE_Optimization(Assembly):
    
    def __init__(self, n=3, m=300):
        super(CADRE_Optimization, self).__init__()
        
        npts = 6
        #add SNOPT driver
        self.add("driver", pyopt_driver.pyOptDriver())
        self.driver.optimizer = "SNOPT"
        
        # Parameters with values common to all analysis points
        self.add("cellInstd", Array(np.ones((7,12)), size=(7,12), dtype=np.float, 
            iotype="in", desc="Cell/Radiator indication", low=0, high=1)
        )
        self.add("finAngle", Float(0., iotype="in", copy=None))
        self.add("antAngle", Float(0., iotype="in", copy=None))
        common_parameters = ["cellInstd", "finAngle", "antAngle"]
        
        # Raw data to load
        solar_raw1 = np.genfromtxt('CADRE/data/Solar/Area10.txt')
        solar_raw2 = np.loadtxt("CADRE/data/Solar/Area_all.txt")
        comm_rawGdata = np.genfromtxt('CADRE/data/Comm/Gain.txt')
        comm_raw = (10**(comm_rawGdata/10.0)).reshape((361,361),order='F')
        power_raw = np.genfromtxt('CADRE/data/Power/curve.dat')
        
        # Initialize analysis points
        for i in xrange(npts):
            print "pt",i
            aname = ''.join(["pt", str(i)])
            self.add(aname, CADRE(n, m, solar_raw1, solar_raw2, 
                                  comm_raw, power_raw))
            for param in common_parameters:
                self.connect(param, '.'.join([aname, param]))
            
            # add parameters to driver
            for k in xrange(12):
                print "CP_Isetpt",k
                for j in xrange(m):
                    param = ''.join(["pt", str(i), ".CP_Isetpt[", str(k), "][",
                                     str(j), "]"])
                    self.driver.add_parameter(param, low=0, high=1)
            for k in xrange(m):
                print "CP_gamma",k
                param = ''.join(["pt",str(i),".CP_gamma[",str(k),"]"])
                self.driver.add_parameter(param, low=0, high=1)
            for k in xrange(m):
                print "CP_comm",k
                param = ''.join(["pt",str(i),".CP_P_comm[",str(k),"]"])
                self.driver.add_parameter(param, low=0, high=1)
                
            # add battery constraints
            constr = ''.join(["pt",str(i),".ConCh >= 0"])
            self.driver.add_constraint(constr)  
            
            constr = ''.join(["pt",str(i),".ConDs >= 0"])
            self.driver.add_constraint(constr)  
            
            constr = ''.join(["pt",str(i),".ConS0 >= 0"])
            self.driver.add_constraint(constr)  
            
            constr = ''.join(["pt",str(i),".ConS1 >= 0"])
            self.driver.add_constraint(constr)  
        
        #add rest of parameters to driver
        for i in xrange(7):
            print "Cellinstd",i
            for k in xrange(12):
                param = ''.join(["cellInstd[",str(i),"][",str(k),"]"])
                self.driver.add_parameter(param, low=0, high=1)
        self.driver.add_parameter("finAngle", low=0, high=np.pi/2.)
        self.driver.add_parameter("antAngle", low=0, high=np.pi)
        
        #add objective
        obj = ''.join([''.join(["-np.sum(pt",str(i),".Data)"]) for i in xrange(npts)])
        self.driver.add_objective(obj)
        
if __name__ == "__main__":
    a = CADRE_Optimization(1500)
    