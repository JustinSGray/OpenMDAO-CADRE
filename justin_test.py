import numpy as np
import pickle

from openmdao.main.api import Assembly, set_as_top

from CADRE.attitude import Attitude_Angular, Attitude_AngularRates, \
     Attitude_Attitude, Attitude_Roll, Attitude_RotationMtx, \
     Attitude_RotationMtxRates, Attitude_Sideslip, Attitude_Torque
from CADRE.battery import BatteryConstraints, BatteryPower, BatterySOC
from CADRE.parameters import BsplineParameters
from CADRE.comm import Comm_AntRotation, Comm_AntRotationMtx, Comm_BitRate, \
     Comm_DataDownloaded, Comm_Distance, Comm_EarthsSpin, Comm_EarthsSpinMtx, \
     Comm_GainPattern, Comm_GSposEarth, Comm_GSposECI, Comm_LOS, Comm_VectorAnt, \
     Comm_VectorBody, Comm_VectorECI, Comm_VectorSpherical
#from MultiPtParameters import MultiPtParameters ??
from CADRE.orbit import Orbit_Initial, Orbit_Dynamics
from CADRE.reactionwheel import ReactionWheel_Motor, ReactionWheel_Power, \
     ReactionWheel_Torque, ReactionWheel_Dynamics
from CADRE.solar import Solar_ExposedArea
from CADRE.sun import Sun_LOS, Sun_PositionBody, Sun_PositionECI,\
     Sun_PositionSpherical
from CADRE.thermal_temperature import ThermalTemperature
from CADRE.power import Power_CellVoltage, Power_SolarPower, Power_Total


NTIME = 4

cadre = set_as_top(Assembly())

cadre.add('comp', Attitude_Attitude(NTIME))

inputs = ['comp.r_e2b_I']
outputs = ['comp.O_RI']

for i_name in inputs: 
    var = cadre.get(i_name)
    val = np.random.random(var.shape)
    cadre.set(i_name, val)


cadre.driver.workflow.add('comp')
cadre.comp.h = .01
cadre.run()
#cadre.driver.workflow.check_gradient(inputs=inputs, outputs=outputs)
cadre.driver.workflow.check_gradient(inputs=inputs, outputs=outputs, adjoint=True)



cadre.driver.update_parameters()
cadre.driver.workflow.config_changed()        
Jn = cadre.driver.workflow.calc_gradient(inputs=inputs,
                                         outputs=outputs,
                                         fd=True)
cadre.driver.update_parameters()
cadre.driver.workflow.config_changed()        
Jf = cadre.driver.workflow.calc_gradient(inputs=inputs,
                                         outputs=outputs)
diff = abs(Jf - Jn)
print diff.max()


diff = np.nan_to_num(abs(Jf - Jn)/Jn)
print diff.max()
