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

cadre.add('comp', ReactionWheel_Dynamics(NTIME))
shape = cadre.comp.w_B.shape
cadre.comp.w_B = np.random.random(shape)
shape = cadre.comp.T_RW.shape
cadre.comp.T_RW = np.random.random(shape)*200
shape = cadre.comp.w_RW0.shape
cadre.comp.w_RW0 = np.random.random(shape)
inputs = ['comp.T_RW']
outputs = ['comp.w_RW']

#cadre.add('comp', BatteryPower(NTIME))
#shape = cadre.comp.SOC.shape
#cadre.comp.SOC = np.random.random(shape)
#shape = cadre.comp.temperature.shape
#cadre.comp.temperature = np.random.random(shape)
#shape = cadre.comp.P_bat.shape
#cadre.comp.P_bat = np.random.random(shape)
#inputs = ['comp.SOC', 'comp.temperature', 'comp.P_bat']
#outputs = ['comp.I_bat']

#cadre.add('comp', ReactionWheel_Power(NTIME))
#shape = cadre.comp.w_RW.shape
#cadre.comp.w_RW = np.random.random(shape)
#shape = cadre.comp.T_RW.shape
#cadre.comp.T_RW = np.random.random(shape)
#inputs = ['comp.T_RW']
#outputs = ['comp.P_RW']

#cadre.add('comp', BatteryConstraints(NTIME))
#shape = cadre.comp.I_bat.shape
#cadre.comp.I_bat = np.random.random(shape)
#shape = cadre.comp.SOC.shape
#cadre.comp.SOC = np.random.random(shape)
#inputs = ['comp.I_bat', 'comp.SOC']
#outputs = ['comp.ConCh', 'comp.ConDs', 'comp.ConS0', 'comp.ConS1']

#cadre.add('comp', Attitude_AngularRates(NTIME))
#shape = cadre.comp.w_B.shape
#cadre.comp.w_B = np.random.random(shape)
#inputs = ['comp.w_B']
#outputs = ['comp.wdot_B']

#cadre.add('comp', ThermalTemperature(NTIME))
#shape = cadre.comp.exposedArea.shape
#cadre.comp.exposedArea = np.random.random(shape)
#shape = cadre.comp.cellInstd.shape
#cadre.comp.cellInstd = np.random.random(shape)
#shape = cadre.comp.LOS.shape
#cadre.comp.LOS = np.random.random(shape)
#shape = cadre.comp.P_comm.shape
#cadre.comp.P_comm = np.random.random(shape)
#inputs = ['comp.exposedArea', 'comp.cellInstd', 
          #'comp.LOS', 'comp.P_comm']
##inputs = ['ThermalTemperature.exposedArea', 'ThermalTemperature.LOS', 'ThermalTemperature.P_comm']
##inputs = ['ThermalTemperature.cellInstd']
#outputs = ['comp.temperature']

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
