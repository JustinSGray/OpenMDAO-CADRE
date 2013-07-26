from CADRE.orbit import Orbit_Dynamics
from pprint import pprint
import numpy as np
import pickle
import pylab

data = pickle.load(open("data1346.pkl", 'rb'))
pkloutput = data["0:r_e2b_I"]
n=1500
comp = Orbit_Dynamics(n)
comp.time_step = 12*3600. / (n - 1)
comp.r_e2b_I0[:3] = data["0:r_e2b_I0"]
comp.r_e2b_I0[3:] = data["0:v_e2b_I0"]
comp.run()

print "loaded:",pkloutput
print 
print "calc:",comp.r_e2b_I
print "error:", np.linalg.norm(comp.r_e2b_I - pkloutput) / np.linalg.norm(pkloutput)




