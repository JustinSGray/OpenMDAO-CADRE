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
m = setd['CP_P_comm'].size

raw1 = np.genfromtxt('CADRE/data/Solar/Area10.txt')
raw2 = np.loadtxt("CADRE/data/Solar/Area_all.txt")

comm_rawGdata = np.genfromtxt('CADRE/data/Comm/Gain.txt')
comm_raw = (10**(comm_rawGdata/10.0)).reshape((361,361),order='F')

power_raw = np.genfromtxt('CADRE/data/Power/curve.dat')

assembly = set_as_top(CADRE(n, m, raw1, raw2, comm_raw, power_raw))

setd['r_e2b_I0'] = np.zeros(6)
setd['r_e2b_I0'][:3] = data[idx+":r_e2b_I0"]
setd['r_e2b_I0'][3:] = data[idx+":v_e2b_I0"]
setd['Gamma'] = data[idx+":gamma"]


assembly.print_set_vals(setvals=setd, printvals="none")
from time import time
tzero = time()
assembly.run()
print "Execution time", time()-tzero

for key in setd.keys():
        print "checking",key
        #print setd[key]
        assembly.print_set_vals(printvals=key, tval=setd[key])
        print

#print assembly.ThermalTemperature.temperature
#print data['3:temperature']
print "done"
# Temp testing stuff
#inputs = ['Power_CellVoltage.LOS', 'Power_CellVoltage.temperature', 
#          'Power_CellVoltage.exposedArea', 'Power_CellVoltage.Isetpt']
#inputs = ['Power_CellVoltage.LOS']
#outputs = ['Power_CellVoltage.V_sol']
#assembly.driver.workflow.check_gradient(inputs=inputs, outputs=outputs)

