import numpy as np
from numpy import array

from CADRE.battery import BatteryConstraints
from CADRE.battery import BatteryPower
from CADRE.battery import BatterySOC

SIZE = 5

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

#battery constraints
io_spec = [
    ('ConCh', (1,)),
    ('ConDs', (1,)),
    ('ConS0', (1,)),
    ('ConS1', (1,)),
    ('I_bat', (SIZE,)),
    ('SOC', (1,SIZE)),
]

"""
#battery Power
io_spec = [
    ('temperature', (5,SIZE)),
    ('P_bat', (SIZE,)),
    ('I_bat', (SIZE,)),
    ('SOC', (1,SIZE)),
]"""




baseline = eval(open('comp_check_baseline.out','rb').read())

comp = BatteryConstraints(n=SIZE)
#comp = BatteryPower(n=SIZE)

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
        if (baseline_value.shape != comp.get(name).shape) and (not(baseline_value.shape==(1,) and isinstance(comp.get(name), float))): 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, comp.get(name).shape
            continue
        error = np.linalg.norm(baseline_value - comp.get(name))
        print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong'

comp.linearize()

arg = {}
result = {}
for name, size in io_spec: 
    arg[name] = comp.get(name)
    result[name] = None


print 5*"#######" 
print "Testing ApplyDer"
print 5*"#######" 
result = comp.applyDer(arg, result)
for name, baseline_value in baseline['applyDer'].iteritems():
    #print name, ": ", baseline_value, result[name]
    if name in outputs: 
        if (baseline_value.shape != result[name].shape) and (not(baseline_value.shape==(1,) and isinstance(result[name], float))): 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, result[name].shape, 
            continue
        error = np.linalg.norm(baseline_value - result[name])
        print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error 
        #print (name+": ").ljust(10), baseline_value, result[name]

print 5*"#######" 
print "Testing ApplyDerT"
print 5*"#######" 
result = comp.applyDerT(arg, result)
for name, baseline_value in baseline['applyDerT'].iteritems():
    if name in inputs: 
        if  baseline_value.shape != result[name].shape and (not(baseline_value.shape==(1,) and isinstance(result[name], float))): 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, result[name].shape
            continue
        #print name, ": ", baseline_value, result[name]
        error = np.linalg.norm(baseline_value - result[name])
        print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error         


