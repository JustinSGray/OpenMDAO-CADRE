import numpy as np
import pickle

from CADRE.CADRE_assembly import CADRE
from openmdao.main.api import Dataflow

cadre = CADRE()
#cadre.Comm_AntRotation.antAngle = 0.1
cadre.run()

#cadre.driver.workflow.check_gradient(inputs=['Orbit_Dynamics.r_e2b_I0'],
#                                     outputs=['Orbit_Dynamics.r_e2b_I'])

#cadre.driver.workflow.check_gradient(inputs=['Comm_DataDownloaded.Dr'],
#                                     outputs=['Comm_DataDownloaded.Data'])

#cadre.driver.workflow.check_gradient(inputs=['Comm_AntRotation.antAngle'],
#                                     outputs=['Comm_AntRotation.q_A'])

#shape = cadre.Comm_AntRotationMtx.q_A.shape
#cadre.Comm_AntRotationMtx.q_A = np.random.random(shape)
#cadre.run()
#cadre.driver.workflow.check_gradient(inputs=['Comm_AntRotationMtx.q_A'],
                                     #outputs=['Comm_AntRotationMtx.O_AB'])

#shape = cadre.Comm_BitRate.P_comm.shape
#cadre.Comm_BitRate.P_comm = np.ones(shape)
#shape = cadre.Comm_BitRate.gain.shape
#cadre.Comm_BitRate.gain = np.random.random(shape)
#shape = cadre.Comm_BitRate.GSdist.shape
#cadre.Comm_BitRate.GSdist = np.random.random(shape)*1e3
#shape = cadre.Comm_BitRate.CommLOS.shape
#cadre.Comm_BitRate.CommLOS = np.random.random(shape)
#cadre.run()
#cadre.driver.workflow.check_gradient(inputs=['Comm_BitRate.P_comm',
                                             #'Comm_BitRate.gain',
                                             #'Comm_BitRate.GSdist',
                                             #'Comm_BitRate.CommLOS'],
                                     #outputs=['Comm_BitRate.Dr'])

shape = cadre.ThermalTemperature.exposedArea.shape
cadre.ThermalTemperature.exposedArea = np.random.random(shape)
shape = cadre.ThermalTemperature.cellInstd.shape
cadre.ThermalTemperature.cellInstd = np.random.random(shape)
shape = cadre.ThermalTemperature.LOS.shape
cadre.ThermalTemperature.LOS = np.random.random(shape)
shape = cadre.ThermalTemperature.P_comm.shape
cadre.ThermalTemperature.P_comm = np.random.random(shape)
#data = pickle.load(open("data1346.pkl", 'rb'))
#cadre.ThermalTemperature.set('exposedArea', data['5:exposedArea'])
#cadre.ThermalTemperature.set('cellInstd', data['cellInstd'])
#cadre.ThermalTemperature.set('LOS', data['5:LOS'])
#cadre.ThermalTemperature.set('P_comm', data['5:P_comm'])
cadre.run()
inputs = ['ThermalTemperature.exposedArea', 'ThermalTemperature.cellInstd', 
          'ThermalTemperature.LOS', 'ThermalTemperature.P_comm']
#inputs = ['ThermalTemperature.exposedArea', 'ThermalTemperature.LOS', 'ThermalTemperature.P_comm']
#inputs = ['ThermalTemperature.cellInstd']
outputs = ['ThermalTemperature.temperature']
cadre.driver.workflow.check_gradient(inputs=inputs,
                                     outputs=outputs)
