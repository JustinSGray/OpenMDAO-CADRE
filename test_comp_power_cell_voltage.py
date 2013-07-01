import numpy as np
from numpy import array

from CADRE.power import Power_CellVoltage

SIZE = 10

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

# LOS, temperature and exposedArea and Isetpt are inputs. V_sol is output
io_spec = [
    ('LOS', (SIZE,)), # all zeros
    ('temperature', (5,SIZE,)), # set to 273
    ('exposedArea', (7,12,SIZE,)), # all zeros
    ('Isetpt', (12,SIZE,)), # all 0.2
    ('V_sol', (12,SIZE,)),
]

baseline = eval(open('comp_check_baseline_power_cell_voltage.out','rb').read())

rawP = np.genfromtxt('../CADRE/CADRE/data/Power/curve.dat')
comp = Power_CellVoltage( n=SIZE, rawP=rawP)
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
        is_error = np.allclose( baseline_value, comp.get(name) )
        print (name+": ").ljust(10), 'OK' if is_error else 'Wrong'


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
    if name in outputs: 
        if (baseline_value.shape != result[name].shape) and (not(baseline_value.shape==(1,) and isinstance(result[name], float))): 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, result[name].shape, 
            continue
        is_error = np.allclose( baseline_value, result[name] )
        print (name+": ").ljust(10), 'OK' if is_error else 'Wrong'



print 5*"#######" 
print "Testing ApplyDerT"
print 5*"#######" 

result = comp.applyDerT(arg, result)
for name, baseline_value in baseline['applyDerT'].iteritems():
    if name in inputs: 
        if  baseline_value.shape != result[name].shape and (not(baseline_value.shape==(1,) and isinstance(result[name], float))): 
            print (name+": ").ljust(10), 'wrong shaped result: ', baseline_value.shape, result[name].shape
            continue
        is_error = np.allclose( baseline_value, result[name] )
        print (name+": ").ljust(10), 'OK' if is_error else 'Wrong'
