from CADRE.orbit import Orbit_Dynamics
from pprint import pprint
import numpy as np
import pickle
import pylab

r_e2b_I0 = np.loadtxt("r_e2b_I0.dat")
r_e2b_I = np.loadtxt("r_e2b_I.dat")

od = Orbit_Dynamics(1500)
od.r_e2b_I0 = r_e2b_I0

od.run()

print np.linalg.norm(od.r_e2b_I - r_e2b_I)/np.linalg.norm(r_e2b_I)


