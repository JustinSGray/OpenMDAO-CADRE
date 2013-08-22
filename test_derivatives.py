
import unittest
import numpy as np
import pickle
import random

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

class Testcase_CADRE(unittest.TestCase):
    """ Test run/step/stop aspects of a simple workflow. """

    def setUp(self):
        """ Called before each test. """
        self.model = set_as_top(Assembly())

    def tearDown(self):
        """ Called after each test. """
        self.model = None
        
    def setup(self, compname, inputs, state0):
        
        self.model.add('comp', eval('%s(NTIME)' % compname))
        self.model.driver.workflow.add('comp')
        
        #data = pickle.load(open("data1346.pkl", 'rb'))
        
        for item in inputs+state0:
            val = self.model.comp.get(item)
            #key = "%d:%s" % (NTIME, item)
            if hasattr(val, 'shape'):
                shape1 = val.shape
                #self.model.comp.set(item, data[key])
                self.model.comp.set(item, np.random.random(shape1))
            else:
                #self.model.comp.set(item, data[key][0])
                self.model.comp.set(item, random.random())
        
    def run_model(self):
    
        self.model.comp.h = 0.01
        self.model.run()
        
    def compare_derivatives(self, var_in, var_out):
        
        wflow = self.model.driver.workflow
        inputs = ['comp.%s' % v for v in var_in]
        outputs = ['comp.%s' % v for v in var_out]
        
        # Numeric
        self.model.driver.update_parameters()
        wflow.config_changed()        
        Jn = wflow.calc_gradient(inputs=inputs,
                                 outputs=outputs,
                                 fd=True)
        #print Jn
        
        # Analytic forward
        self.model.driver.update_parameters()
        wflow.config_changed()        
        Jf = wflow.calc_gradient(inputs=inputs,
                                 outputs=outputs)
        
        #print Jf
        
        diff = abs(Jf - Jn)
        assert_rel_error(self, diff.max(), 0.0, 1e-5)
        
        # Analytic adjoint
        self.model.driver.update_parameters()
        wflow.config_changed()        
        Ja = wflow.calc_gradient(inputs=inputs,
                                 outputs=outputs,
                                 mode='adjoint')
        
        #print Ja
        
        diff = abs(Ja - Jn)
        assert_rel_error(self, diff.max(), 0.0, 1e-5)
    
    def test_Comm_DataDownloaded(self):
        
        compname = 'Comm_DataDownloaded'
        inputs = ['Dr']
        outputs = ['Data']
        state0 = ['Data0']
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_AntRotation(self):
        
        compname = 'Comm_AntRotation'
        inputs = ['antAngle']
        outputs = ['q_A']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_AntRotationMtx(self):
        
        compname = 'Comm_AntRotationMtx'
        inputs = ['q_A']
        outputs = ['O_AB']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_BitRate(self):
        
        compname = 'Comm_BitRate'
        inputs = ['P_comm', 'gain', 'GSdist', 'CommLOS']
        outputs = ['Dr']
        state0 = []
        
        self.setup(compname, inputs, state0)
        
        # These need to be a certain magnitude so it doesn't blow up
        shape = self.model.comp.P_comm.shape
        self.model.comp.P_comm = np.ones(shape)
        shape = self.model.comp.GSdist.shape
        self.model.comp.GSdist = np.random.random(shape)*1e3
        
        self.run_model()
        
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_Distance(self):
        
        compname = 'Comm_Distance'
        inputs = ['r_b2g_A']
        outputs = ['GSdist']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_EarthsSpin(self):
        
        compname = 'Comm_EarthsSpin'
        inputs = ['t']
        outputs = ['q_E']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_EarthsSpinMtx(self):
        
        compname = 'Comm_EarthsSpinMtx'
        inputs = ['q_E']
        outputs = ['O_IE']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_GainPattern(self):
        
        compname = 'Comm_GainPattern'
        inputs = ['azimuthGS', 'elevationGS']
        outputs = ['gain']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_GSposEarth(self):
        
        compname = 'Comm_GSposEarth'
        inputs = ['lon', 'lat', 'alt']
        outputs = ['r_e2g_E']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_GSposECI(self):
        
        compname = 'Comm_GSposECI'
        inputs = ['O_IE', 'r_e2g_E']
        outputs = ['r_e2g_I']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_LOS(self):
        
        compname = 'Comm_LOS'
        inputs = ['r_b2g_I', 'r_e2g_I']
        outputs = ['CommLOS']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_VectorAnt(self):
        
        compname = 'Comm_VectorAnt'
        inputs = ['r_b2g_B', 'O_AB']
        outputs = ['r_b2g_A']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_VectorBody(self):
        
        compname = 'Comm_VectorBody'
        inputs = ['r_b2g_I', 'O_BI']
        outputs = ['r_b2g_B']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_VectorECI(self):
        
        compname = 'Comm_VectorECI'
        inputs = ['r_e2g_I', 'r_e2b_I']
        outputs = ['r_b2g_I']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_Comm_VectorSpherical(self):
        
        compname = 'Comm_VectorSpherical'
        inputs = ['r_b2g_A']
        outputs = ['azimuthGS', 'elevationGS']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)
        
    def test_ThermalTemperature(self):
        
        compname = 'ThermalTemperature'
        inputs = ['exposedArea', 'cellInstd', 'LOS', 'P_comm']
        outputs = ['temperature']
        state0 = ['T0']
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)

    def test_Attitude_Angular(self):
        
        compname = 'Attitude_Angular'
        inputs = ['O_BI', 'Odot_BI']
        outputs = ['w_B']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)     

    def test_Attitude_AngularRates(self):

        compname = 'Attitude_AngularRates'
        inputs = ['h', 'w_B']
        outputs = ['wdot_B']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)        

    def test_Attitude_Attitude(self):

        compname = 'Attitude_Attitude'
        inputs = ['r_e2b_I']
        outputs = ['O_RI']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)

    def test_Attitude_Roll(self):
        
        compname = 'Attitude_Roll'
        inputs = ['Gamma']
        outputs = ['O_BR']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)       

    def test_Attitude_RotationMtx(self):

        compname = 'Attitude_RotationMtx'
        inputs = ['O_BR', 'O_RI']
        outputs = ['O_BI']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_Attitude_RotationMtxRates(self):
        
        compname = 'Attitude_RotationMtxRates'
        inputs = ['h', 'O_BI']
        outputs = ['Odot_BI']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_Attitude_Sideslip(self):
        
        compname = 'Attitude_Sideslip'
        inputs = ['r_e2b_I', 'O_BI']
        outputs = ['v_e2b_B']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_Attitude_Torque(self):
        
        compname = 'Attitude_Torque'
        inputs = ['w_B', 'wdot_B']
        outputs = ['T_tot']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)     

    def test_Sun_LOS(self):
        
        compname = 'Sun_LOS'
        inputs = ['r_e2b_I', 'r_e2s_I']
        outputs = ['LOS']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)      

    def test_Sun_PositionBody(self):
        
        compname = 'Sun_PositionBody'
        inputs = ['O_BI', 'r_e2s_I']
        outputs = ['r_e2s_B']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)   

    def test_Sun_PositionECI(self):
        
        compname = 'Sun_PositionECI'
        inputs = ['t', 'LD']
        outputs = ['r_e2s_I']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)   

    def test_Solar_ExposedArea(self):

        compname = 'Solar_ExposedArea'
        inputs = ['finAngle', 'azimuth', 'elevation']
        outputs = ['exposedArea']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs) 

    #def test_Power_CellVoltage(self):
    #    #fix
    #    compname = 'Power_CellVoltage'
    #    inputs = ['LOS', 'temperature', 'exposedArea', 'Isetpt']
    #    outputs = ['V_sol']
    #    state0 = []
    #    
    #    self.setup(compname, inputs, state0)
    #    self.run_model()
    #    self.compare_derivatives(inputs, outputs)  

    def test_Power_SolarPower(self):

        compname = 'Power_SolarPower'
        inputs = ['V_sol', 'Isetpt']
        outputs = ['P_sol']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_Power_Total(self):

        compname = 'Power_Total'
        inputs = ['P_sol', 'P_comm', 'P_RW']
        outputs = ['P_bat']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)  

    def test_ReactionWheel_Motor(self):

        compname = 'ReactionWheel_Motor'
        inputs = ['T_RW', 'w_B', 'w_RW']
        outputs = ['T_m']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)   

    def test_ReactionWheel_Dynamics(self):

        compname = 'ReactionWheel_Dynamics'
        inputs = ['w_B', 'T_RW']
        outputs = ['w_RW']
        state0 = ['w_RW0']
        
        self.setup(compname, inputs, state0)
        
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_ReactionWheel_Power(self):

        compname = 'ReactionWheel_Power'
        inputs = ['w_RW', 'T_RW']
        outputs = ['P_RW']
        state0 = []
        
        self.setup(compname, inputs, state0)
        
        self.run_model()
        self.compare_derivatives(inputs, outputs)   

    def test_ReactionWheel_Torque(self):

        compname = 'ReactionWheel_Torque'
        inputs = ['T_tot']
        outputs = ['T_RW']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_BatterySOC(self):

        compname = 'BatterySOC'
        inputs = ['P_bat', 'temperature']
        outputs = ['SOC']
        state0 = ['iSOC']
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_BatteryPower(self):

        compname = 'BatteryPower'
        inputs = ['SOC', 'temperature', 'P_bat']
        outputs = ['I_bat']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_BatteryConstraints(self):

        compname = 'BatteryConstraints'
        inputs = ['I_bat', 'SOC']
        outputs = ['ConCh', 'ConDs', 'ConS0', 'ConS1']
        state0 = []
        
        self.setup(compname, inputs, state0)
        self.run_model()
        self.compare_derivatives(inputs, outputs)    

    def test_Orbit_Dynamics(self):
        # This comp has no inputs, so no derivs.
        pass

        #compname = 'Orbit_Dynamics'
        #inputs = ['']
        #outputs = ['r_e2b_I']
        #state0 = ['r_e2b_I0']
        
        #self.setup(compname, inputs, state0)
        #self.run_model()
        #self.compare_derivatives(inputs, outputs)    

        
if __name__ == "__main__":
    
    unittest.main()
