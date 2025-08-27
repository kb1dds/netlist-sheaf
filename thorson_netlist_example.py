# Example based upon Figure 4 in
# Thorson et al., "Dynamic structural equation models synthesize ecosystem dynamics constrained by ecological mechanisms", DOI: 10.1111/2041-210X.14289
#
# This model uses the netlist sheaf construction.

import csv
import json
import numpy as np
import scipy.linalg
import pysheaf as ps
import matplotlib.pyplot as plt
import netlist_sheaf
from collections import defaultdict

# Names in the Thorson data disagree with their paper
names_dict={
    "log_CP" : "ColdPool",
    "log_seaice" : "SeaIce",
    "log_Cfall" : "Copepods",
    "log_Esummer" : "Krill",
    "log_RperS" : "Survival",
    "SSB" : "Spawners",    
    "log_PercentCop" : "Diet_Cop",
    "log_PercentEuph" : "Diet_Krill"
    }

def lag_fcn(npts,ar,x):
    '''Convenience function for building lag matrices'''
    if ar:
        return scipy.linalg.toeplitz(np.concatenate(([0],x[0:ar],np.zeros((npts-ar-1,)))),np.zeros((npts,)))
    else: # Disable autoregression
        return np.eye(npts)

# Load the sheaf structure
parts=json.load(open('thorson_parts.json'))
nets=json.load(open('thorson_nets.json'))

startyear = 1963
npts = 61
ar = 0
shf = netlist_sheaf.NetlistSheaf(parts,nets,npts=npts,ar=ar,lag_fcn=lag_fcn)

#pos=nx.layout.spring_layout(shf)
#nx.draw_networkx_labels(shf,pos)
#nx.draw_networkx_edges(shf,pos)
#plt.show()

# Load the data; adding in new cells to house the known values
measurements = defaultdict(list)
with open('bering_sea.csv') as fp:
    reader = csv.DictReader(fp, delimiter = ',')
    idx = 0
    for row in reader:
        for data_name,sheaf_name in names_dict.items():
            if row[data_name] != 'NA':
                year = startyear+idx
                datum = str(year) + "-" + sheaf_name
                value = float(row[data_name])

                measurements[sheaf_name].append((year,value))
                
                shf.AddCell(datum, ps.Cell('datum',dataDimension=1))
                shf.GetCell(datum).SetDataAssignment(ps.Assignment('datum', np.array([value])))
                shf.AddCoface(sheaf_name,datum,
                              ps.Coface('net','datum', lambda x,idx=idx : x[idx] ))
        
        idx = idx + 1

# Load in DSEM predictions
dsem_values = defaultdict(list)
with open('dsem_bering_sea.csv') as fp:
    reader = csv.DictReader(fp, delimiter = ',')
    idx = 0
    for row in reader:
        for data_name,sheaf_name in names_dict.items():
            if row[data_name] != 'NA':
                year = startyear+idx
                datum = str(year) + "-" + sheaf_name
                value = float(row[data_name])

                dsem_values[sheaf_name].append((year,value))        
        idx = idx + 1

for c in shf.GetCellIndexList():
    shf.MaximallyExtendCell(c)

print('===')
print("Consistency radius before optimization: {}".format(shf.ComputeConsistencyRadius()))

print('Optimizing...')
shf.FuseAssignment()

print("Consistency radius after optimization: {}".format(shf.ComputeConsistencyRadius()))

print('===')
print("Path coefficients:")
print('---')

for c in shf.GetCellIndexList():
    if "-pc" in c:
        print("{} : {}".format(c,shf.GetCell(c).mDataAssignment.mValue))

print('===')
print("Offsets:")
print('---')

for c in shf.GetCellIndexList():
    if c == 'ColdPool-block':
        print("{} : {}".format(c,shf.GetCell(c).mDataAssignment.mValue[0:2]))
    elif "-block" in c:
        print("{} : {}".format(c,shf.GetCell(c).mDataAssignment.mValue[0]))

print('===')
print("Autoregressive cofficients:")
print('---')

for data_name,sheaf_name in names_dict.items():
    try:
        print("{} : {}".format(sheaf_name,shf.GetCell(sheaf_name + "-lag").mDataAssignment.mValue[0:ar]))
    except KeyError:
        pass


print('')
print('==========')
print('All stalks')
print('---')
for c in shf.GetCellIndexList():
    print("{} : {}".format(c,shf.GetCell(c).mDataAssignment.mValue))

years=[startyear + year for year in range(npts)]

fig = plt.figure()
axes = fig.subplots(3,3, sharex='all')
row = 0
col = 0
for data_name,sheaf_name in names_dict.items():
    axes[row,col].plot([x for x,y in dsem_values[sheaf_name]],[y for x,y in dsem_values[sheaf_name]],'k:',label='DSEM')
    axes[row,col].plot(years,shf.GetCell(sheaf_name).mDataAssignment.mValue,'b',label='Sheaf')
    axes[row,col].plot([x for x,y in measurements[sheaf_name]],[y for x,y in measurements[sheaf_name]],'*',label='Measurement')
    axes[row,col].set_title(sheaf_name)
    axes[row,col].legend()

    row += 1
    if row > 2:
        row = 0
        col += 1

plt.show()
