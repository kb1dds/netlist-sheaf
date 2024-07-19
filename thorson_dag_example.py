# Example based upon Figure 4 in
# Thorson et al., "Dynamic structural equation models synthesize ecosystem dynamics constrained by ecological mechanisms", DOI: 10.1111/2041-210X.14289
#
# This model uses the graph as shown in the figure as the Hasse diagram of the
# base space poset (with some modification to permit path coefficients to vary;
# see below), even though naive diagram does not commute. This can be
# remedied by the product trick at the "Survival" cell.
# The point of the netlist sheaf (not this script; see the other one!)
# is that it handles this a little more uniformly.
#
# Unlike what is shown in the figure, each of the upstream cells (everthing
# but "survival") has an additional cell that carries the corresponding
# path coefficient.  This allows path coefficients to be adjusted easily by
# the researcher, including allowing them to participate in optimizations

import numpy as np
import networkx as nx
import pysheaf as ps
import matplotlib.pyplot as plt

# Timeseries length
npts = 10

# Default path coefficients
seaIceColdPool = 0.6
coldPoolCopepods = 1.79
coldPoolKrill = 0.18
copepodsDietCop = 0.29
dietCopSurvival = 0.15
krillDietKrill = 0.06
dietKrillSurvival = 0.13
spawnersSurvival = -0.59

shf = ps.Sheaf()

# Construct all cells and their path coefficient connections

# SeaIce holds in its stalk path coefficients for: 
# 0: SeaIce->ColdPool,
# 1: ColdPool->Copepods,
# 2: ColdPool->Krill,
# 3: Copepods->Diet_Cop,
# 4: Diet_Cop->Survival,
# 5: Krill->Diet_Krill,
# 6: Diet_Krill->Survival (NB: if we had commutativity, this would be redundant)
shf.AddCell("SeaIce", ps.Cell(None,
                              dataDimension=npts+7))
shf.GetCell("SeaIce").mOptimizationCell=True
shf.GetCell("SeaIce").SetDataAssignment(ps.Assignment(None,np.concatenate(([seaIceColdPool,coldPoolCopepods,coldPoolKrill,copepodsDietCop,dietCopSurvival,krillDietKrill,dietKrillSurvival],np.zeros((npts,))))))
shf.AddCell("SeaIce-pathcoeff",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("SeaIce-pathcoeff").mOptimizationCell=False
shf.GetCell("SeaIce-pathcoeff").SetDataAssignment(ps.Assignment(None,np.array([seaIceColdPool])))
shf.AddCoface("SeaIce","SeaIce-pathcoeff",
              ps.Coface(None,None,lambda x: x[0]))

# ColdPool holds in its stalk path coefficients for: 
# 0: ColdPool->Copepods,
# 1: ColdPool->Krill,
# 2: Copepods->Diet_Cop,
# 3: Diet_Cop->Survival,
# 4: Krill->Diet_Krill,
# 5: Diet_Krill->Survival
shf.AddCell("ColdPool", ps.Cell(None,
                                dataDimension=npts+6))
shf.GetCell("ColdPool").mOptimizationCell=True
shf.GetCell("ColdPool").SetDataAssignment(ps.Assignment(None,np.concatenate(([coldPoolCopepods,coldPoolKrill,copepodsDietCop,dietCopSurvival,krillDietKrill,dietKrillSurvival],np.zeros((npts,))))))
shf.AddCell("ColdPool-pathcoeff1",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("ColdPool-pathcoeff1").mOptimizationCell=False
shf.GetCell("ColdPool-pathcoeff1").SetDataAssignment(ps.Assignment(None,np.array([coldPoolCopepods])))
shf.AddCoface("ColdPool","ColdPool-pathcoeff1",
              ps.Coface(None,None,lambda x: x[0]))
shf.AddCell("ColdPool-pathcoeff2",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("ColdPool-pathcoeff2").mOptimizationCell=False
shf.GetCell("ColdPool-pathcoeff2").SetDataAssignment(ps.Assignment(None,np.array([coldPoolKrill])))
shf.AddCoface("ColdPool","ColdPool-pathcoeff2",
              ps.Coface(None,None,lambda x: x[1]))

# Copepods holds in its stalk path coefficients for: 
# 0: Copepods->Diet_Cop,
# 1: Diet_Cop->Survival,
shf.AddCell("Copepods", ps.Cell(None,
                                dataDimension=npts+2))
shf.GetCell("Copepods").mOptimizationCell=True
shf.GetCell("Copepods").SetDataAssignment(ps.Assignment(None,np.concatenate(([copepodsDietCop,dietCopSurvival],np.zeros((npts,))))))
shf.AddCell("Copepods-pathcoeff",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("Copepods-pathcoeff").mOptimizationCell=False
shf.GetCell("Copepods-pathcoeff").SetDataAssignment(ps.Assignment(None,np.array([copepodsDietCop])))
shf.AddCoface("Copepods","Copepods-pathcoeff",
              ps.Coface(None,None,lambda x: x[0]))

# Krill holds in its stalk path coefficients for: 
# 0: Krill->Diet_Krill,
# 1: Diet_Krill->Survival
shf.AddCell("Krill", ps.Cell(None,
                              dataDimension=npts+2))
shf.GetCell("Krill").mOptimizationCell=True
shf.GetCell("Krill").SetDataAssignment(ps.Assignment(None,np.concatenate(([krillDietKrill,dietKrillSurvival],np.zeros((npts,))))))
shf.AddCell("Krill-pathcoeff",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("Krill-pathcoeff").mOptimizationCell=False
shf.GetCell("Krill-pathcoeff").SetDataAssignment(ps.Assignment(None,np.array([krillDietKrill])))
shf.AddCoface("Krill","Krill-pathcoeff",
              ps.Coface(None,None,lambda x: x[0]))

# Diet_Cop holds in its stalk path coefficients for: 
# 0: Diet_Cop->Survival,
shf.AddCell("Diet_Cop", ps.Cell(None,
                              dataDimension=npts+1))
shf.GetCell("Diet_Cop").mOptimizationCell=True
shf.GetCell("Diet_Cop").SetDataAssignment(ps.Assignment(None,np.concatenate(([dietCopSurvival],np.zeros((npts,))))))
shf.AddCell("Diet_Cop-pathcoeff",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("Diet_Cop-pathcoeff").mOptimizationCell=False
shf.GetCell("Diet_Cop-pathcoeff").SetDataAssignment(ps.Assignment(None,np.array([dietCopSurvival])))
shf.AddCoface("Diet_Cop","Diet_Cop-pathcoeff",
              ps.Coface(None,None,lambda x: x[0]))

# Diet_Krill holds in its stalk path coefficients for: 
# 0: Diet_Krill->Survival
shf.AddCell("Diet_Krill", ps.Cell(None,
                              dataDimension=npts+1))
shf.GetCell("Diet_Krill").mOptimizationCell=True
shf.GetCell("Diet_Krill").SetDataAssignment(ps.Assignment(None,np.concatenate(([dietKrillSurvival],np.zeros((npts,))))))
shf.AddCell("Diet_Krill-pathcoeff",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("Diet_Krill-pathcoeff").mOptimizationCell=False
shf.GetCell("Diet_Krill-pathcoeff").SetDataAssignment(ps.Assignment(None,np.array([dietKrillSurvival])))
shf.AddCoface("Diet_Krill","Diet_Krill-pathcoeff",
              ps.Coface(None,None,lambda x: x[0]))

# Spawners holds in its stalk path coefficients for: 
# 0: Spawners->Survival
shf.AddCell("Spawners", ps.Cell(None,
                              dataDimension=npts+1))
shf.GetCell("Spawners").mOptimizationCell=True
shf.GetCell("Spawners").SetDataAssignment(ps.Assignment(None,np.concatenate(([spawnersSurvival],np.zeros((npts,))))))
shf.AddCell("Spawners-pathcoeff",ps.Cell(None,
                                       dataDimension=1))
shf.GetCell("Spawners-pathcoeff").mOptimizationCell=False
shf.GetCell("Spawners-pathcoeff").SetDataAssignment(ps.Assignment(None,np.array([spawnersSurvival])))
shf.AddCoface("Spawners","Spawners-pathcoeff",
              ps.Coface(None,None,lambda x: x[0]))

# Product trick applied to the "Survival" node
shf.AddCell("SurvivalDC", ps.Cell(None,
                              dataDimension=npts))
shf.GetCell("SurvivalDC").mOptimizationCell=True
shf.GetCell("SurvivalDC").SetDataAssignment(ps.Assignment(None,np.zeros((npts,))))
shf.AddCell("SurvivalDK", ps.Cell(None,
                              dataDimension=npts))
shf.GetCell("SurvivalDK").mOptimizationCell=True
shf.GetCell("SurvivalDK").SetDataAssignment(ps.Assignment(None,np.zeros((npts,))))
shf.AddCell("SurvivalS", ps.Cell(None,
                              dataDimension=npts))
shf.GetCell("SurvivalS").mOptimizationCell=True
shf.GetCell("SurvivalS").SetDataAssignment(ps.Assignment(None,np.zeros((npts,))))
shf.AddCell("Survival-prod", ps.Cell(None,
                                     dataDimension=npts*3))
shf.GetCell("Survival-prod").mOptimizationCell=True
shf.GetCell("Survival-prod").SetDataAssignment(ps.Assignment(None,np.zeros((3*npts,))))
shf.AddCell("Survival", ps.Cell(None,
                                     dataDimension=npts))
shf.GetCell("Survival").mOptimizationCell=True
shf.GetCell("Survival").SetDataAssignment(ps.Assignment(None,np.zeros((npts,))))
shf.AddCoface("Survival-prod","SurvivalDC",
              ps.Coface(None,None,lambda x: x[0:npts]))
shf.AddCoface("Survival-prod","SurvivalDK",
              ps.Coface(None,None,lambda x: x[npts:2*npts]))
shf.AddCoface("Survival-prod","SurvivalS",
              ps.Coface(None,None,lambda x: x[2*npts:]))
shf.AddCoface("Survival-prod","Survival",
              ps.Coface(None,None,lambda x: x[0:npts]+x[npts:2*npts]+x[2*npts:]))

# Construct causal connections between cells
shf.AddCoface("SeaIce","ColdPool",
              ps.Coface(None,None,
                        lambda x: np.concatenate((x[1:7],x[0]*x[7:]))))

shf.AddCoface("ColdPool","Copepods",
              ps.Coface(None,None,
                        lambda x: np.concatenate((x[2:4],x[0]*x[6:]))))

shf.AddCoface("ColdPool","Krill",
              ps.Coface(None,None,
                        lambda x: np.concatenate((x[4:6],x[1]*x[6:]))))

shf.AddCoface("Copepods","Diet_Cop",
              ps.Coface(None,None,
                        lambda x: np.concatenate((x[1:2],x[0]*x[2:]))))

shf.AddCoface("Diet_Cop","SurvivalDC",
              ps.Coface(None,None,
                        lambda x: x[0]*x[1:]))

shf.AddCoface("Krill","Diet_Krill",
              ps.Coface(None,None,
                        lambda x: np.concatenate((x[1:2],x[0]*x[2:]))))

shf.AddCoface("Diet_Krill","SurvivalDK",
              ps.Coface(None,None,
                        lambda x: x[0]*x[1:]))

shf.AddCoface("Spawners","SurvivalS",
              ps.Coface(None,None,
                        lambda x: x[0]*x[1:]))


for c in shf.GetCellIndexList():
    shf.MaximallyExtendCell(c)

print(shf.ComputeConsistencyRadius())
