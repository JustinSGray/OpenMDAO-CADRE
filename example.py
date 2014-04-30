from openmdao.lib.casehandlers.api import CSVCaseRecorder
from CADRE import CADRE_Optimization

print "setting up"
top = CADRE_Optimization(n=1500, m=300)
top.recorders = [CSVCaseRecorder(filename='CADRE.csv')]
printvars = []
for var in ['Data', 'ConCh', 'ConDs', 'ConS0', 'ConS1', 'SOC']:
    printvars += ["pt" + str(i) + ".Data" for i in xrange(6)]
top.printvars = printvars
print "running"
top.run()
