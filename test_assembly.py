from CADRE.CADRE_assembly import CADRE
import networkx as nx
import numpy as np
import pprint

raw1 = np.genfromtxt('CADRE/data/Solar/Area10.txt')
raw2 = np.loadtxt("CADRE/data/Solar/Area_all.txt")

comm_rawGdata = np.genfromtxt('CADRE/data/Comm/Gain.txt')
comm_raw = (10**(comm_rawGdata/10.0)).reshape((361,361),order='F')

power_raw = np.genfromtxt('CADRE/data/Power/curve.dat')
assembly = CADRE(1500, raw1, raw2, comm_raw, power_raw)

pprint.pprint(assembly.varnames)

#assembly.get(assembly.varnames[var]).get(var)

graph = assembly._depgraph._graph

graph.remove_node("@xin")
graph.remove_node("@xout")
graph.remove_node("@bin")
graph.remove_node("@bout")
graph.remove_node("driver")

ag = nx.to_agraph(graph)
ag.layout('dot')
ag.draw('design.png')

