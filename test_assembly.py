from CADRE.CADRE_assembly import CADRE
import networkx as nx

assembly = CADRE()
graph = assembly._depgraph._graph

ag = nx.to_agraph(graph)
ag.layout('dot')
ag.draw('design.png')