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
import numpy.polynomial
from collections import defaultdict

# Basic parameters of the signal
npts = 12
ar = 3

# Construct an autoregressive signal
roots = np.random.randn(ar) # ... roots of the characteristic polynomial
ar_coefs = -numpy.polynomial.Polynomial.fromroots(roots).coef[:-1][::-1]
#roots = np.roots(np.concatenate(([1],-ar_coefs)))
restriction = scipy.linalg.toeplitz(np.concatenate(([0],ar_coefs,np.zeros((npts-ar-1,)))),np.zeros((npts,)))
x=roots[0]**np.arange(0,npts)

#print(restriction)
print('Roots of characteristic equation {}'.format(roots))
print('Signal is {}'.format(x))
#print('Applying restriction {}'.format(np.dot(restriction,x)))
#print('Autoregression equation LHS is {} '.format(x-np.dot(restriction,x)))
print('Expected error is {}'.format(np.linalg.norm((x-np.dot(restriction,x))[ar:]))) # Should be small!

print('---')
 
# Construct the sheaf structure
parts=json.load(open('ar_parts.json'))
nets=json.load(open('ar_nets.json'))

shf = netlist_sheaf.NetlistSheaf(parts,nets,npts=npts,ar=ar)

# The timeseries, as variables
shf.GetCell('signal').SetDataAssignment(ps.Assignment("net",x))
shf.GetCell('signal-lagvar').SetDataAssignment(ps.Assignment("net",x[ar:]))

# Extra copies of the timeseries in other stalks; can improve convergence of the optimizer
#shf.GetCell('signal-block').SetDataAssignment(ps.Assignment("part",x))
#shf.GetCell('signal-block').mOptimizationCell = False
#shf.GetCell('signal-lag').SetDataAssignment(ps.Assignment("part",np.concatenate((np.zeros((ar,)),x))))

print('Optimizing...')
shf.FuseAssignment()


for c in shf.GetCellIndexList():
    print('Value at {} is {}'.format(c,shf.GetCell(c).mDataAssignment.mValue))
    print('  local consistency radius is {}'.format(shf.ComputeStarLocalConsistencyRadius(c)))

print('---')

# These should be similar, but are often not... because the optimizer tends to get "lost"
print('Designed coefficients {}'.format(ar_coefs))
print('Estimated coefficients {}'.format(shf.GetCell('signal-lag').mDataAssignment.mValue[0:ar]))

restriction2 = scipy.linalg.toeplitz(np.concatenate(([0],shf.GetCell('signal-lag').mDataAssignment.mValue[0:ar],np.zeros((npts-ar-1,)))),np.zeros((npts,)))

# These should be similar
print('Consistency radius is {}'.format(shf.ComputeConsistencyRadius()))
print('Residual is {}'.format(np.linalg.norm((x-np.dot(restriction2,x))[ar:],ord=2)))

# Now let's verify the sheaf was actually correct...
# Putting in the correct values for the AR coefficients should yield a consistency radius
# that is the same as the expected error above
shf.GetCell('signal').SetDataAssignment(ps.Assignment("net",x))
shf.GetCell('signal-lagvar').SetDataAssignment(ps.Assignment("net",x[ar:]))
shf.GetCell('signal-block').SetDataAssignment(ps.Assignment("part",x))
shf.GetCell('signal-lag').SetDataAssignment(ps.Assignment("part",np.concatenate((ar_coefs,x))))

for c in shf.GetCellIndexList():
    shf.MaximallyExtendCell(c)

print('Consistency radius of global section {}'.format(shf.ComputeConsistencyRadius()))
if np.abs(shf.ComputeConsistencyRadius()-np.linalg.norm((x-np.dot(restriction,x))[ar:]))/np.linalg.norm((x-np.dot(restriction,x))[ar:])<0.01:
    print('Correct to within 1 percent')
else:
    print('Global section is in error!')
