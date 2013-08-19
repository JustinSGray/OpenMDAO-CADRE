import numpy as np

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

shape = cadre.Comm_AntRotationMtx.q_A.shape
cadre.Comm_AntRotationMtx.q_A = np.random.random(shape)
cadre.run()
cadre.driver.workflow.check_gradient(inputs=['Comm_AntRotationMtx.q_A'],
                                     outputs=['Comm_AntRotationMtx.O_AB'])

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
