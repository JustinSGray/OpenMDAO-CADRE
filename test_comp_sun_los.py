import numpy as np
from numpy import array

from CADRE.sun_los import Sun_LOS

SIZE = 10

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

io_spec = [
    ('LOS', (SIZE,)),
    ('r_e2b_I', (6,SIZE,)),
    ('r_e2s_I', (3,SIZE,)),
]

baseline = eval(open('comp_check_baseline.out','rb').read())

comp = Sun_LOS(n=SIZE)
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
        error = np.linalg.norm(baseline_value - result[name])
        print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error 


print 5*"#######" 
print "Testing ApplyDerT"
print 5*"#######" 
result = comp.applyDerT(arg, result)
for name, baseline_value in baseline['applyDerT'].iteritems():
    if name in inputs: 
        #print name, ": ", baseline_value, result[name]
        error = np.linalg.norm(baseline_value - result[name])
        print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error         


