from CADRE.CADRE_assembly import CADRE
from openmdao.main.api import Dataflow

cadre = CADRE()
cadre.driver.workflow = Dataflow()
#cadre.driver.workflow.add(['Orbit_Dynamics'])
cadre.driver.workflow.add(['Comm_DataDownloaded'])
cadre.run()

#cadre.driver.workflow.check_gradient(inputs=['Orbit_Dynamics.r_e2b_I0'],
#                                     outputs=['Orbit_Dynamics.r_e2b_I'])

cadre.driver.workflow.check_gradient(inputs=['Comm_DataDownloaded.Dr'],
                                     outputs=['Comm_DataDownloaded.Data'])