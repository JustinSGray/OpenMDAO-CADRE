from CADRE.reactionwheel import ReactionWheel_Dynamics
from pprint import pprint
import numpy as np
import pickle
import pylab

data = pickle.load(open("data1346.pkl", 'rb'))
pkloutput = data["0:w_RW"]
n=1500
comp = ReactionWheel_Dynamics(n)
comp.h = 12*3600. / (n - 1)

comp.w_B = data["0:w_B"]
comp.T_RW = data["0:T_RW"]
#comp.w_RW0 = data["0:w_RW0"]
comp.run()

print "loaded:",pkloutput
print 
print "calc:",comp.w_RW
print "error:", np.linalg.norm(comp.w_RW - pkloutput) / np.linalg.norm(pkloutput)




