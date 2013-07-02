from openmdao.main.api import Assembly

from attitude import Attitude_Angular, Attitude_AngularRates, Attitude_Attitude, \
     Attitude_Roll, Attitude_RotationMtx, \
     Attitude_RotationMtxRates, Attitude_Sideslip, Attitude_Torque
from battery import BatteryConstraints, BatteryPower, BatterySOC
#from bsplineparameters import BsplineParameters ??
from comm import Comm_AntRotation, Comm_AntRotationMtx, Comm_BitRate, \
     Comm_DataDownloaded, Comm_Distance, Comm_EarthsSpin, Comm_EarthsSpinMtx, \
     Comm_GainPattern, Comm_GSposEarth, Comm_GSposECI, Comm_LOS, Comm_VectorAnt, \
     Comm_VectorBody, Comm_VectorECI, Comm_VectorSpherical
#from MultiPtParameters import MultiPtParameters ??
from orbit import Orbit_Initial
#from parameters import Parameters ??
from reactionwheel import ReactionWheel_Motor, ReactionWheel_Power, \
     ReactionWheel_Torque
#from solar import Solar
from sun import Sun_LOS, Sun_PositionBody, Sun_PositionECI, Sun_PositionSpherical
from thermal_temperature import ThermalTemperature

#rk4 components:
#Comm_DataDownloaded, BatterySOC, ThermalTemperature, Orbit_Dynamics

class CADRE(Assembly):
    """
    OpenMDAO implementation of the CADRE model
    """
    def __init__(self, n=2):
        super(CADRE, self).__init__()
        
        # Additude components
        self.add("Attitude_Angular", Attitude_Angular(n))
        self.driver.workflow.add("Attitude_Angular")
        
        self.add("Attitude_AngularRates", Attitude_AngularRates(n))
        self.driver.workflow.add("Attitude_AngularRates")
        
        self.add("Attitude_Attitude", Attitude_Attitude(n))
        self.driver.workflow.add("Attitude_Attitude")
        
        self.add("Attitude_Roll", Attitude_Roll(n))
        self.driver.workflow.add("Attitude_Roll")
        
        self.add("Attitude_RotationMtx", Attitude_RotationMtx(n))
        self.driver.workflow.add("Attitude_RotationMtx")
        
        self.add("Attitude_RotationMtxRates", Attitude_RotationMtxRates(n))
        self.driver.workflow.add("Attitude_RotationMtxRates")
        
        self.add("Attitude_Sideslip", Attitude_Sideslip(n))
        self.driver.workflow.add("Attitude_Sideslip")
        
        self.add("Attitude_Torque", Attitude_Torque(n))
        self.driver.workflow.add("Attitude_Torque")
        
        # Battery components
        self.add("BatteryConstraints", BatteryConstraints(n))
        self.driver.workflow.add("BatteryConstraints")
        
        self.add("BatteryPower", BatteryPower(n))
        self.driver.workflow.add("BatteryPower")
        
        self.add("BatterySOC", BatterySOC(n))
        self.driver.workflow.add("BatterySOC")
        
        # Comm components
        self.add("Comm_AntRotation", Comm_AntRotation(n))
        self.driver.workflow.add("Comm_AntRotation")
        
        self.add("Comm_AntRotationMtx", Comm_AntRotationMtx(n))
        self.driver.workflow.add("Comm_AntRotationMtx")

        self.add("Comm_BitRate", Comm_BitRate(n))
        self.driver.workflow.add("Comm_BitRate")
        
        self.add("Comm_DataDownloaded", Comm_DataDownloaded(n))
        self.driver.workflow.add("Comm_DataDownloaded")
        
        self.add("Comm_Distance", Comm_Distance(n))
        self.driver.workflow.add("Comm_Distance")
        
        self.add("Comm_EarthsSpin", Comm_EarthsSpin(n))
        self.driver.workflow.add("Comm_EarthsSpin")

        self.add("Comm_EarthsSpinMtx", Comm_EarthsSpinMtx(n))
        self.driver.workflow.add("Comm_EarthsSpinMtx")
        
        self.add("Comm_GainPattern", Comm_GainPattern(n))
        self.driver.workflow.add("Comm_GainPattern")
        
        self.add("Comm_GSposEarth", Comm_GSposEarth(n))
        self.driver.workflow.add("Comm_GSposEarth")
        
        self.add("Comm_GSposECI", Comm_GSposECI(n))
        self.driver.workflow.add("Comm_GSposECI")
        
        self.add("Comm_LOS", Comm_LOS(n))
        self.driver.workflow.add("Comm_LOS")
        
        self.add("Comm_VectorAnt", Comm_VectorAnt(n))
        self.driver.workflow.add("Comm_VectorAnt")
        
        self.add("Comm_VectorBody", Comm_VectorBody(n))
        self.driver.workflow.add("Comm_VectorBody")
        
        self.add("Comm_VectorECI", Comm_VectorECI(n))
        self.driver.workflow.add("Comm_VectorECI")
        
        self.add("Comm_VectorSpherical", Comm_VectorSpherical(n))
        self.driver.workflow.add("Comm_VectorSpherical")
        
		# Orbit components
        self.add("Orbit_Initial", Orbit_Initial())
        self.driver.workflow.add("Orbit_Initial")
        
        # Reaction wheel components
        self.add("ReactionWheel_Motor", ReactionWheel_Motor(n))
        self.driver.workflow.add("ReactionWheel_Motor")
        
        self.add("ReactionWheel_Power", ReactionWheel_Power(n))
        self.driver.workflow.add("ReactionWheel_Power")
        
        self.add("ReactionWheel_Torque", ReactionWheel_Torque(n))
        self.driver.workflow.add("ReactionWheel_Torque")

        # Sun components
        self.add("Sun_LOS", Sun_LOS(n))
        self.driver.workflow.add("Sun_LOS")
        
        self.add("Sun_PositionBody", Sun_PositionBody(n))
        self.driver.workflow.add("Sun_PositionBody")

        self.add("Sun_PositionECI", Sun_PositionECI(n))
        self.driver.workflow.add("Sun_PositionECI")
        
        self.add("Sun_PositionSpherical", Sun_PositionSpherical(n))
        self.driver.workflow.add("Sun_PositionSpherical")

        # Thermal temp components
        self.add("ThermalTemperature", ThermalTemperature(n))
        self.driver.workflow.add("ThermalTemperature")
        
        self.make_connections()
        
    def make_connections(self):
        """
        Collects the names of all input and output variables for all
        components within the assembly (drivers excluded). 
        
        Then establishes connections between
        any output variable and input variable that has the same name, so 
        long as the variable name does not exist as an output to more than
        a single component (so excludes default outputs)
        """
        inputs, outputs = {}, {}
        for compname in self.list_components():
            if compname == "driver":
                continue
            
            comp_inputs = self.get(compname).list_inputs()
            
            for input_name in comp_inputs:
                if input_name not in inputs:
                    inputs[input_name] = [compname]
                else:
                    inputs[input_name].append(compname)
                
            comp_outputs = self.get(compname).list_outputs()
            
            for output_name in comp_outputs:
                if output_name not in outputs:
                    outputs[output_name] = [compname]
                else:
                    outputs[output_name].append(compname)
            
        print
        for varname in outputs.keys():
            comps = outputs[varname]
            if len(comps) > 1:
                continue
            frompath = '.'.join([comps[0], varname])
            
            if varname in inputs:
                for compname in inputs[varname]:
                    topath = '.'.join([compname, varname])
                    self.connect(frompath, topath)
                    print "Connecting", frompath, "to", topath, "..."