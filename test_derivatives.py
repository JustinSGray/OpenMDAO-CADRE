
import unittest
from numpy.random import random

from openmdao.main.api import Assembly, set_as_top
from openmdao.util.testutil import assert_rel_error

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


NTIME = 5

class Testcase_provideJ(unittest.TestCase):
    """ Test run/step/stop aspects of a simple workflow. """

    def test_Comm_DataDownloaded(self):
        
        model = set_as_top(Assembly())
        model.add('comp', Comm_DataDownloaded(NTIME))
        model.driver.workflow.add('comp')
        
        for item in ['Dr', 'Data0']:
            val = model.comp.get(item)
            shape1 = val.shape
            model.comp.set(item, random(shape1))
        
        model.comp.h = 0.01
        model.run()
        wflow = model.driver.workflow
        
        # Numeric
        model.driver.update_parameters()
        wflow.config_changed()        
        Jn = wflow.calc_gradient(inputs=['comp.Dr'],
                                 outputs=['comp.Data'],
                                 fd=True)
        
        print Jn
        
        # Analytic forward
        model.driver.update_parameters()
        wflow.config_changed()        
        Jf = wflow.calc_gradient(inputs=['comp.Dr'],
                                 outputs=['comp.Data'])
        
        print Jf
        
        diff = abs(Jf - Jn)
        assert_rel_error(self, diff.max(), 0.0, 1e-6)
        
        # Analytic adjoint
        model.driver.update_parameters()
        wflow.config_changed()        
        Ja = wflow.calc_gradient(inputs=['comp.Dr'],
                                 outputs=['comp.Data'],
                                 mode='adjoint')
        
        print Ja
        
        diff = abs(Ja - Jn)
        assert_rel_error(self, diff.max(), 0.0, 1e-6)
        
    def test_Comm_Comm_AntRotation(self):
        
        model = set_as_top(Assembly())
        model.add('comp', Comm_AntRotation(NTIME))
        model.driver.workflow.add('comp')
        
        for item in ['antAngle']:
            val = model.comp.get(item)
            shape1 = val.shape
            model.comp.set(item, random(shape1))
        
        model.run()
        wflow = model.driver.workflow
        
        # Numeric
        model.driver.update_parameters()
        wflow.config_changed()        
        Jn = wflow.calc_gradient(inputs=['comp.antAngle'],
                                 outputs=['comp.q_A'],
                                 fd=True)
        
        print Jn
        
        # Analytic forward
        model.driver.update_parameters()
        wflow.config_changed()        
        Jf = wflow.calc_gradient(inputs=['comp.Dr'],
                                 outputs=['comp.Data'])
        
        print Jf
        
        diff = abs(Jf - Jn)
        assert_rel_error(self, diff.max(), 0.0, 1e-6)
        
        # Analytic adjoint
        model.driver.update_parameters()
        wflow.config_changed()        
        Ja = wflow.calc_gradient(inputs=['comp.Dr'],
                                 outputs=['comp.Data'],
                                 mode='adjoint')
        
        print Ja
        
        diff = abs(Ja - Jn)
        assert_rel_error(self, diff.max(), 0.0, 1e-6)
        
if __name__ == "__main__":
    unittest.main()        