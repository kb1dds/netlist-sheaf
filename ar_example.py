# Autoregressive model example

import csv
import json
import numpy as np
import scipy.linalg
import networkx as nx
import pysheaf as ps
import matplotlib.pyplot as plt
import netlist_sheaf
import numpy.random
from collections import defaultdict

# Load the sheaf structure
parts=json.load(open('ar_parts.json'))
nets=json.load(open('ar_nets.json'))

npts = 6
ar = 3
shf = netlist_sheaf.NetlistSheaf(parts,nets,npts=npts,ar=ar)

#pos=nx.layout.spring_layout(shf)
#nx.draw_networkx_labels(shf,pos)
#nx.draw_networkx_edges(shf,pos)
#plt.show()

# Construct an autoregressive signal
ar_coefs = np.random.randn(ar)
roots = np.roots(np.concatenate(([1],-ar_coefs)))
restriction = scipy.linalg.toeplitz(np.concatenate(([0],ar_coefs,np.zeros((npts-ar-1,)))),np.zeros((npts,)))
x=roots[0]**np.arange(0,npts)

print(restriction)
print(ar_coefs)
print(roots)
print('Signal is {}'.format(x))
print('Applying restriction {}'.format(np.dot(restriction,x)))
print('Autoregression equation LHS is {} '.format(x-np.dot(restriction,x)))
print('Expected error is {}'.format(np.linalg.norm((x-np.dot(restriction,x))[ar:])))

shf.GetCell('signal').SetDataAssignment(ps.Assignment("net",x))
shf.GetCell('signal-lagvar').SetDataAssignment(ps.Assignment("net",x[ar:]))

print('Optimizing...')
shf.FuseAssignment()

print('Consistency radius is {}'.format(shf.ComputeConsistencyRadius()))

for c in shf.GetCellIndexList():
    print('Value at {} is {}'.format(c,shf.GetCell(c).mDataAssignment.mValue))

print('Designed coefficients {}'.format(ar_coefs))
print('Estimated coefficients {}'.format(shf.GetCell('signal-lag').mDataAssignment.mValue[0:ar]))

restriction2 = scipy.linalg.toeplitz(np.concatenate(([0],shf.GetCell('signal-lag').mDataAssignment.mValue[0:ar],np.zeros((npts-ar-1,)))),np.zeros((npts,)))

print('Residual is {}'.format(np.linalg.norm((x-np.dot(restriction2,x))[ar:])))
