<<<<<<< HEAD
from CADRE.orbit import Orbit_Dynamics
=======
from CADRE.reactionwheel import ReactionWheel_Dynamics
>>>>>>> 763a137ee915413eb6e851d79486a3b524f32667
from pprint import pprint
import numpy as np
import pickle
import pylab

data = pickle.load(open("data1346.pkl", 'rb'))
<<<<<<< HEAD
pkloutput = data["0:r_e2b_I"]
n=1500
comp = Orbit_Dynamics(n)
comp.time_step = 12*3600. / (n - 1)
comp.r_e2b_I0[:3] = data["0:r_e2b_I0"]
comp.r_e2b_I0[3:] = data["0:v_e2b_I0"]
=======
pkloutput = data["0:w_RW"]
n=1500
comp = ReactionWheel_Dynamics(n)
comp.h = 12*3600. / (n - 1)

comp.w_B = data["0:w_B"]
comp.T_RW = data["0:T_RW"]
#comp.w_RW0 = data["0:w_RW0"]
>>>>>>> 763a137ee915413eb6e851d79486a3b524f32667
comp.run()

print "loaded:",pkloutput
print 
<<<<<<< HEAD
print "calc:",comp.r_e2b_I
print "error:", np.linalg.norm(comp.r_e2b_I - pkloutput) / np.linalg.norm(pkloutput)
=======
print "calc:",comp.w_RW
print "error:", np.linalg.norm(comp.w_RW - pkloutput) / np.linalg.norm(pkloutput)
>>>>>>> 763a137ee915413eb6e851d79486a3b524f32667




