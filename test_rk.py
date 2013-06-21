import numpy as np
from numpy import array

import timeit

from CADRE.battery import BatteryConstraints
from CADRE.battery import BatteryPower
from CADRE.battery import BatterySOC

from CADRE.thermal_temperature import ThermalTemperature

SIZE = 5

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

#battery SOC
io_spec = [
    ('temperature', (5, SIZE)),
    ('P_bat', (SIZE,)),
    ('iSOC', (1,)),
    ('SOC', (1,SIZE)),
]

#thermal_temperature
io_spec = [
    ('temperature', (5, SIZE)),
    ('exposedArea', (7,12,SIZE)), 
    ('cellInstd',(7,12)), 
    ('LOS', (SIZE, )), 
    ('P_comm', (SIZE, ))
]


baseline = eval(open('comp_check_baseline.out','rb').read())

comp = BatterySOC(n_times=SIZE, time_step=1)
comp = ThermalTemperature(n_times=SIZE, time_step=1)
inputs = comp.list_inputs()
outputs = comp.list_outputs()

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
        if (not(baseline_value.shape==(1,) and isinstance(comp.get(name), float))) and baseline_value.shape != comp.get(name).shape: 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, result[name].shape
            continue
        is_error = np.allclose(baseline_value, comp.get(name))
        print (name+": ").ljust(10), 'OK' if is_error else 'Wrong' 
        #print baseline_value, "\n\n" ,comp.get(name)

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
        if (baseline_value.shape != result[name].shape) and (not(baseline_value.shape==(1,) and isinstance(result[name], float))): 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, result[name].shape
            continue
        is_error = np.allclose(result[name],baseline_value,rtol=1e-1,atol=5)
        print (name+": ").ljust(10), 'OK' if is_error else 'Wrong' 
        #print (name+": ").ljust(10), "\n" , baseline_value, "\n\n", result[name]

print 5*"#######" 
print "Testing ApplyJT"
print 5*"#######" 
result = comp.applyJT(arg, result)
for name, baseline_value in baseline['applyDerT'].iteritems():
    if name in inputs: 
        if (baseline_value.shape != result[name].shape) and (not(baseline_value.shape==(1,) and isinstance(result[name], float))): 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, result[name].shape
            continue
        #error = np.linalg.norm(baseline_value - result[name])
        is_error = np.allclose(baseline_value, result[name],rtol=1e-1,atol=5)
        print (name+": ").ljust(10), 'OK' if is_error else 'Wrong'    
        #print (name+": ").ljust(10), baseline_value, "\n\n", result[name]


