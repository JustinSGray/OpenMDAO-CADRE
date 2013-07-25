from CADRE.CADRE_assembly import CADRE
from pprint import pprint
import pickle

idx = '5'

setd = {}
data = pickle.load(open("data1346.pkl", 'rb'))
for key in data.keys():
    if key[0] == idx or not key[0].isdigit():
        if not key[0].isdigit():
            shortkey = key
        else:
            shortkey = key[2:]
        if data[key].shape == (1,) and shortkey != "iSOC": #set floats correctly
            setd[shortkey] = data[key][0]
        else:
            setd[shortkey] = data[key]
n = setd['P_comm'].size

assembly = CADRE(n=n)

assembly.print_set_vals(setvals=setd, printvals="none")

assembly.run()

print

for key in data.keys():
    if key[0] == idx or not key[0].isdigit():
        shortkey = key[2:]
        if not key[0].isdigit():
            shortkey = key
        
        print "checking",key
        print data[key]
        assembly.print_set_vals(printvals=shortkey)
        print



