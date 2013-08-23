from openmdao.main.api import Assembly
from openmdao.main.datatypes.api import Float, Array, Int
import numpy as np
from CADRE_assembly import CADRE
from pyopt_driver import pyopt_driver

class CADRE_Optimization(Assembly):
    
    def __init__(self, n=3):
        super(CADRE_Optimization, self).__init__()
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
        for i in xrange(6):
            aname = ''.join(["pt", str(i)])
            self.add(aname, CADRE(n, solar_raw1, solar_raw2, 
                                  comm_raw, power_raw))
            for param in common_parameters:
                self.connect(param, '.'.join([aname, param]))
        
        #add parameters
        #self.driver.add_parameter()
        
        #add objective
        #self.driver.add_objective()
        
        #add constraints
        #self.driver.add_constraint()
        
if __name__ == "__main__":
    a = CADRE_Optimization(1500)
    