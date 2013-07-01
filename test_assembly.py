from CADRE.CADRE_assembly import CADRE
import networkx as nx

assembly = CADRE()
graph = assembly._depgraph._graph

graph.remove_node("@xin")
graph.remove_node("@xout")
graph.remove_node("@bin")
graph.remove_node("@bout")
graph.remove_node("driver")

ag = nx.to_agraph(graph)
ag.layout('dot')
ag.draw('design.png')