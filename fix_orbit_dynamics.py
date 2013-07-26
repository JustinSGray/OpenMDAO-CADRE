from CADRE.power import Power_CellVoltage
from pprint import pprint
import numpy as np
import pickle
import pylab

data = pickle.load(open("data1346.pkl", 'rb'))
pkloutput = data["0:V_sol"]
n=1500
comp = Power_CellVoltage(n)
comp.time_step = 12*3600. / (n - 1)

comp.LOS = data["0:LOS"]
comp.temperature = data["0:temperature"]
comp.exposedArea = data["0:exposedArea"]
comp.Isetpt = data["0:Isetpt"]
comp.run()

print "loaded:",pkloutput
print 
print "calc:",comp.V_sol
print "error:", np.linalg.norm(comp.V_sol - pkloutput) / np.linalg.norm(pkloutput)




