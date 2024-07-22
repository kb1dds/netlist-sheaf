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
npts = 30
ar = 3

# Construct an autoregressive signal

# We will start with the roots of the characteristic polynomial to ensure they're real
roots = 2/np.pi*np.arctan(np.random.randn(ar))

# From the roots, we recover the autoregressive coefficients
ar_coefs = -numpy.polynomial.Polynomial.fromroots(roots).coef[:-1][::-1]

# This is a copy of the lagging map; it's the correct non-identity restriction in the sheaf
restriction = scipy.linalg.toeplitz(np.concatenate(([0],ar_coefs,np.zeros((npts-ar-1,)))),np.zeros((npts,)))

# Build the signal as a random linear combination of the roots
x=np.zeros((npts,))
for i in range(ar):
    x = x + np.random.randn(1)*roots[i]**np.arange(0,npts)

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

print('Drift in signal-lag stalk: {}'.format(np.linalg.norm(shf.GetCell('signal-lag').mDataAssignment.mValue[ar:]-x)))

print('---')

ar_est = shf.GetCell('signal-lag').mDataAssignment.mValue[0:ar]
roots_est = numpy.polynomial.Polynomial(np.concatenate(([1],-ar_est))[::-1]).roots()

# The estimated coefficients and should be similar, though may be in a different order
print('Designed coefficients {}'.format(ar_coefs))
print('Estimated coefficients {}'.format(ar_est))
print('Designed roots {}'.format(roots))
print('Estimated roots {}'.format(roots_est))

# Constructing the only non-identity restriction map for use outside of the Sheaf instance
restriction2 = scipy.linalg.toeplitz(np.concatenate(([0],ar_est,np.zeros((npts-ar-1,)))),np.zeros((npts,)))

# These should be similar
print('Consistency radius is {}'.format(shf.ComputeConsistencyRadius()))
print('Residual is {}'.format(np.linalg.norm((x-np.dot(restriction2,x))[ar:],ord=2)))
print(np.dot(restriction,x)[ar:])
print(np.dot(restriction2,x)[ar:])
print('Norm of difference is {}'.format(np.linalg.norm(np.dot(restriction-restriction2,x)[ar:])))

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
if np.abs(shf.ComputeConsistencyRadius()-np.linalg.norm((x-np.dot(restriction,x))[ar:]))/np.linalg.norm((x-np.dot(restriction,x))[ar:])>0.01:
    print('Global section is in error!')
