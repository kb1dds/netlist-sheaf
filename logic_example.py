import json
import numpy as np
import networkx as nx
import netlist_sheaf
import pysheaf as ps
import matplotlib.pyplot as plt

parts=json.load(open('parts.json'))
nets=json.load(open('nets.json'))

shf = netlist_sheaf.NetlistSheaf(parts,nets)

#pos=nx.layout.spring_layout(shf)
#nx.draw_networkx_labels(shf,pos)
#nx.draw_networkx_edges(shf,pos)
#plt.show()

for nd in shf.nodes():
    shf.GetCell(nd).mOptimizationCell = True

for a in [0,1]:
    for b in [0,1]:
        for c in [0,1]:
            q = min((1-a*b)+c,1)

            shf.GetCell('A').SetDataAssignment(ps.Assignment('net',np.array((a,))))
            shf.GetCell('A').mOptimizationCell = False
            shf.GetCell('B').SetDataAssignment(ps.Assignment('net',np.array((b,))))
            shf.GetCell('B').mOptimizationCell = False
            shf.GetCell('C').SetDataAssignment(ps.Assignment('net',np.array((c,))))
            shf.GetCell('C').mOptimizationCell = False

            shf.FuseAssignment()

            qp = shf.GetCell('Q').mDataAssignment
            
            print('A={}, B={}, C={},  Q={} =?= {} '.format(a,b,c,q,qp))
            
