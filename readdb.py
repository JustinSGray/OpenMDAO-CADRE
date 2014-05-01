import csv
import numpy as np
import pylab

import pprint 

from openmdao.lib.casehandlers.dbcase import list_db_vars, DBCaseIterator


#print list_db_vars('CADRE.db')

dataIter = DBCaseIterator('CADRE.db')


X, Y, Z = [], [], []

pcom = []

start_time = None
case_count = 0
seconds_per_case = []


for case in dataIter: 
    if ('driver.workflow.itername' in case.keys()): 
        #print case['driver.workflow.itername']
        
        if start_time is None: 
            start_time = case.timestamp

        case_count+=1
        seconds_per_case.append((case.timestamp-start_time)/case_count)

        #pprint.pprint(case.keys())
        data = [case["pt" + str(i) + ".Data"][0][1499] for i in xrange(6)]
        sumdata = sum([float(i) for i in data if i])
        c1 = [case["Constraint ( pt" + str(i) + ".ConCh<=0 )"][0] for i in xrange(6)]
        c2 = [case["Constraint ( pt" + str(i) + ".ConDs<=0 )"][0] for i in xrange(6)]
        c3 = [case["Constraint ( pt" + str(i) + ".ConS0<=0 )"][0] for i in xrange(6)]
        c4 = [case["Constraint ( pt" + str(i) + ".ConS1<=0 )"][0] for i in xrange(6)]
        c5 = [case["Constraint ( pt" + str(i) + ".SOC[0][0]=pt" + str(i) + ".SOC[0][-1] )"][0]
              for i in xrange(6)]
        # c1_f = np.all([float(i) < 0 for i in c1 if i])
        # c2_f = np.all([float(i) < 0 for i in c2 if i])
        # c3_f = np.all([float(i) < 0 for i in c3 if i])
        # c4_f = np.all([float(i) < 0 for i in c4 if i])
        # c5_f = np.all([float(i) < 0 for i in c4 if i])

        c1_f = sum([float(i) for i in c1 if i])
        c2_f = sum([float(i) for i in c2 if i])
        c3_f = sum([float(i) for i in c3 if i])
        c4_f = sum([float(i) for i in c4 if i])
        c5_f = sum([float(i) for i in c5 if i])

        feasible = [c1_f, c2_f,  c3_f, c4_f, c5_f]

        X.append(sumdata), Y.append(sum(feasible)), Z.append(feasible)

        pcom.append([float(case["pt5.CP_gamma"][i])
                    for i in xrange(300)])

        # print sumdata, sum(feasible), max(feasible) #,[ '%.1f' % i for i in
        # feasible]
        #print sumdata

end_time = case.timestamp

tot_time = end_time - start_time
print "total time: %f hours (%f sec)"%(tot_time/3600., tot_time)

#pylab.figure()
#pylab.plot(pcom[-1])


pylab.figure()
pylab.plot(seconds_per_case)

Z = np.array(Z)
if not len(Z):
    print "no data yet..."
    quit()
pylab.figure()
pylab.subplot(311)
pylab.title("total data")
pylab.plot(X, 'b')
pylab.plot([0, len(X)], [3e4, 3e4], 'k--', marker="o")
pylab.subplot(312)
pylab.title("Sum of Constraints")
pylab.plot([0, len(Y)], [0, 0], 'k--', marker="o")
pylab.plot(Y, 'k')
pylab.subplot(313)
pylab.title("Max of Constraints")
pylab.plot([0, len(Z)], [0, 0], 'k--')
pylab.plot(Z[:, 0], marker="o", label="c1")
pylab.plot(Z[:, 1], marker="o", label="c2")
pylab.plot(Z[:, 2], marker="o", label="c3")
pylab.plot(Z[:, 3], marker="o", label="c4")
pylab.plot(Z[:, 4], marker="o", label="c5")
pylab.legend(loc="best")
pylab.show()