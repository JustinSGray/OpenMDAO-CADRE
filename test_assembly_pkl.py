from openmdao.main.api import set_as_top
from CADRE.CADRE_assembly import CADRE
from pprint import pprint
import numpy as np
import pickle
#import pylab

idx = '5'

setd = {}
data = pickle.load(open("data1346.pkl", 'rb'))

for key in data.keys():
    if key[0] == idx or not key[0].isdigit():
        if not key[0].isdigit():
            shortkey = key
        else:
            shortkey = key[2:]
        if data[key].shape == (1,) and shortkey != "iSOC": #set floats correctly
            setd[shortkey] = data[key][0]
        else:
            setd[shortkey] = data[key]

n = setd['P_comm'].size
assembly = set_as_top(CADRE(n=n))

setd['r_e2b_I0'] = np.zeros(6)
setd['r_e2b_I0'][:3] = data[idx+":r_e2b_I0"]
setd['r_e2b_I0'][3:] = data[idx+":v_e2b_I0"]
setd['Gamma'] = data[idx+":gamma"]


assembly.print_set_vals(setvals=setd, printvals="none")
assembly.run()

for key in setd.keys():
        print "checking",key
        #print setd[key]
        assembly.print_set_vals(printvals=key, tval=setd[key])
        print

print assembly.ThermalTemperature.temperature
print data['3:temperature']



