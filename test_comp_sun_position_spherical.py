import numpy as np
from numpy import array

from CADRE.sun import Sun_PositionSpherical

SIZE = 10

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

# azimuth and elevation are outputs
io_spec = [
    ('azimuth', (SIZE,)),
    ('elevation', (SIZE,)),
    ('r_e2s_B', (3,SIZE,)),
]

baseline = eval(open('comp_check_baseline_position_spherical.out','rb').read())

comp = Sun_PositionSpherical(n=SIZE)
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
