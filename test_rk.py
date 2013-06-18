import numpy as np
from numpy import array

import timeit

from CADRE.battery import BatteryConstraints
from CADRE.battery import BatteryPower
from CADRE.battery import BatterySOC

SIZE = 5

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

#battery SOC
io_spec = [
    ('temperature', (SIZE, )),
    ('P_bat', (SIZE,)),
    ('iSOC', (1,)),
    ('SOC', (1,SIZE)),
]

io_map = {
    'SOC':'y'
}


baseline = eval(open('comp_check_baseline.out','rb').read())

comp = BatterySOC(n_times=SIZE, time_step=1)
inputs = ['iSOC', 'P_bat', 'temperature']#comp.list_inputs()
outputs = ['SOC'] #comp.list_outputs()

for name,size in io_spec: 
    if name in inputs: 
        value = baseline['execute'][name]
        comp.set(name,value)

comp.run()        

print 5*"#######" 
print "Testing execute"
print 5*"#######" 

for name,size in io_spec: 
    if name in outputs: 
        baseline_value = baseline['execute'][name]
        error = np.linalg.norm(baseline_value - comp.get(name))
        print (name+": ").ljust(10), 'OK' if error < 1e-9 else 'Wrong -- error: %.20f'%error
        #print baseline_value, comp.get(name)
comp.linearize()

arg = {}
result = {}
for name, size in io_spec: 
    arg[name] = comp.get(name)
    result[name] = None


print 5*"#######" 
print "Testing ApplyJ"
print 5*"#######" 
result = comp.applyJ(arg, result)
for name, baseline_value in baseline['applyDer'].iteritems():
    #print name, ": ", baseline_value, result[name]
    if name in outputs: 
        error = np.linalg.norm(baseline_value - result[name])
        print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error 
        #print (name+": ").ljust(10), baseline_value, result[name]


print 5*"#######" 
print "Testing ApplyJT"
print 5*"#######" 
result = comp.applyJT(arg, result)
for name, baseline_value in baseline['applyDerT'].iteritems():
    if name in inputs: 
        error = np.linalg.norm(baseline_value - result[name])
        print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error         
        #print (name+": ").ljust(10), baseline_value, result[name]


