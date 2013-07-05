import numpy as np
from numpy import array

import timeit

from CADRE.battery import BatteryConstraints
from CADRE.battery import BatteryPower
from CADRE.battery import BatterySOC

from CADRE.thermal_temperature import ThermalTemperature

from CADRE.reactionwheel import ReactionWheel_Dynamics

from CADRE.comm import Comm_DataDownloaded

from CADRE.orbit import Orbit_Dynamics

SIZE = 5

############################################################################
# Edit the io_spec to match your component -- same as from other test file
############################################################################

#battery SOC

io_specs = []

io_spec = [
    ('temperature', (5, SIZE)),
    ('P_bat', (SIZE,)),
    ('iSOC', (1,)),
    ('SOC', (1,SIZE)),
]

io_specs.append(io_spec)

#thermal_temperature
io_spec = [
    ('temperature', (5, SIZE)),
    ('exposedArea', (7,12,SIZE)), 
    ('cellInstd',(7,12)), 
    ('LOS', (SIZE, )), 
    ('P_comm', (SIZE, ))
]

io_specs.append(io_spec)

#reactionwheel_dynamics

io_spec = [
    ('w_B', (3,SIZE)),
    ('T_RW', (3,SIZE)),
    ('w_RW', (3,SIZE)),
]

io_specs.append(io_spec)

io_spec = [
    ('Dr',(SIZE,)),
    ('Data',(1,SIZE)),
]

io_specs.append(io_spec)

io_spec = [
    ('r_e2b_I',(6,SIZE)),
    ('r_e2b_I0',(6,)), 
]

io_specs.append(io_spec)


baselines = []
baselines.append('comp_check_SOC')
baselines.append('comp_check_thermal_temp')
baselines.append('comp_check_reactionwheel_dynamics')
baselines.append('comp_check_comm_datadownloaded')
baselines.append('comp_check_orbit_dynamics')



comps = []

comp = BatterySOC(n_times=SIZE, time_step=1)
comps.append(comp)

comp = ThermalTemperature(n_times=SIZE, time_step=1)
comps.append(comp)

comp = ReactionWheel_Dynamics(n_times=SIZE, time_step=1)
comps.append(comp)

comp = Comm_DataDownloaded(n_times=SIZE, time_step=1)
comps.append(comp)

comp = Orbit_Dynamics(n_times=SIZE, time_step=1)
comps.append(comp)




for comp,io_spec,baseline_file in zip(comps,io_specs,baselines): 

    print 5*"#######" 
    print 5*"#######" 
    print baseline_file
    print 5*"#######" 
    print 5*"#######" 

    baseline = eval(open(baseline_file,'rb').read())
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
            is_ok = np.allclose(baseline_value, comp.get(name))
            print (name+": ").ljust(10), 'OK' if is_ok else 'Wrong' 
            if not is_ok: 
                print "HERE"
                print baseline_value, "\n\n" ,comp.get(name)
                exit()

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
            is_ok = np.allclose(result[name],baseline_value,rtol=1e-1,atol=5)
            print (name+": ").ljust(10), 'OK' if is_ok else 'Wrong' 
            if not is_ok: 
                print (name+": ").ljust(10), "\n" , baseline_value, "\n\n", result[name]
                exit()

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
            is_ok = np.allclose(baseline_value, result[name],rtol=1e-1,atol=5)
            print (name+": ").ljust(10), 'OK' if is_ok else 'Wrong'    
            
            if not is_ok: 
                print (name+": ").ljust(10), "\n" , baseline_value, "\n\n", result[name]
                exit()
    print 5*"\n"
