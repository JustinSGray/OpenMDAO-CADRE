import numpy as np
from numpy import array

from CADRE.ReactionWheel_Motor import ReactionWheel_Motor
from CADRE.ReactionWheel_Power import ReactionWheel_Power
from CADRE.ReactionWheel_Torque import ReactionWheel_Torque

SIZE = 5

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

#ReactionWheel_Motor
'''
io_spec = [
    ('T_RW', (3,SIZE)),
    ('w_B', (3,SIZE)),
    ('w_RW', (3,SIZE)),
    ('T_m', (3,SIZE)),
]
'''
#ReactionWheel_Power
io_spec = [
    ('w_RW', (3,SIZE)),
    ('T_RW', (3,SIZE)),
    ('P_RW', (3,SIZE)),
]
'''
#ReactionWheel_Torque
io_spec = [
    ('T_tot', (3,SIZE)),
    ('T_RW', (3,SIZE))
]
'''
baseline = eval(open('comp_check_baseline.out','rb').read())

#comp = ReactionWheel_Motor(n=SIZE)
comp = ReactionWheel_Power(n=SIZE)
#comp = ReactionWheel_Torque(n=SIZE)
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
        #error = np.linalg.norm(baseline_value - comp.get(name))
        #print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong'
        print baseline_value - comp.get(name) #TEST
        error2 = np.nanmax( np.absolute( baseline_value - comp.get(name) ) / np.absolute( baseline_value ) )
        print (name+": ").ljust(10), 'OK' if error2 < 1e-5 else 'Wrong: %f'%ferror2

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
        #error = np.linalg.norm(baseline_value - result[name])
        #print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error 
        print baseline_value - result[name] #TEST        
        error2 = np.nanmax( np.absolute( baseline_value - result[name] ) / np.absolute( baseline_value ) )
        print (name+": ").ljust(10), 'OK' if error2 < 1e-5 else 'Wrong: %f'%error2

print 5*"#######" 
print "Testing ApplyDerT"
print 5*"#######" 
result = comp.applyDerT(arg, result)
for name, baseline_value in baseline['applyDerT'].iteritems():
    if name in inputs: 
        #print name, ": ", baseline_value, result[name]
        #error = np.linalg.norm(baseline_value - result[name])
        #print (name+": ").ljust(10), 'OK' if error < 1e-5 else 'Wrong: %f'%error
        print baseline_value, result[name] #TEST
        error2 = np.nanmax( np.absolute( baseline_value - result[name] ) / np.absolute( baseline_value ) )
        print (name+": ").ljust(10), 'OK' if error2 < 1e-5 else 'Wrong: %f'%error2        
print "COMPLETE"


