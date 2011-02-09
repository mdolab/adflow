#! /usr/bin/env python

# Import standard modules
import sys
import os

# Import class representing SUmb
sys.path.append(os.path.abspath('..'))
from sumbsolution import SUmbSolution

# Get an instance of the SUmb class
myflow = SUmbSolution('startfile_unsteady')

# Run time steps / iterations.
#   i) With no arguments the values are based on the input file;
#  ii) With only one argument, it overrides the number of time steps;
# iii) A second argument overrrides the number of MG cycles.
myflow.RunIterations()
#myflow.RunIterations(3)
#myflow.RunIterations(3,2)

#print 'RSDMassRMS = ',myflow.GetConvergenceHistory('RSDMassRMS')

# Get the family names from the CGNS file
#families = myflow.GetMesh().GetFamilyNames()
#print 'Family names =',families

# Get the surface coordinates
#xyz_all = myflow.GetMesh().GetSurfaceCoordinates(1)
#xyz_wing = myflow.GetMesh().GetSurfaceCoordinates(1,'Wing')

# Get the surface loads
#loads_all = myflow.GetSurfaceLoads(1)
#loads_wing = myflow.GetSurfaceLoads(1,'Wing')

# Get the surface patche dimensions
#patch_dim_all = myflow.GetMesh().GetSurfacePatchDimensions()
#patch_dim_wing = myflow.GetMesh().GetSurfacePatchDimensions('Wing')

# Get the surface Cp
#cp_all = myflow.GetSurfaceCp(1)
#cp_wing = myflow.GetSurfaceCp(1,'Wing')

# Write CGNS files
# (filenames are optional, by default uses names in startfile)
#myflow.GetMesh().WriteMeshFile()
#myflow.GetMesh().WriteMeshFile('mesh.cgns')
myflow.WriteVolumeSolutionFile()
#myflow.WriteVolumeSolutionFile('vol_sol.cgns')
#myflow.WriteSurfaceSolutionFile()
#myflow.WriteSurfaceSolutionFile('surf_sol.cgns')

# Insert new coordinates

# Run more time steps / iterations
#   i) With no arguments the values are based on the previous call;
#  ii) With only one argument, it overrides the number of time steps;
# iii) A second argument overrrides the number of MG cycles.
#myflow.RunIterations()
myflow.RunIterations(2,3)

# Output convergence history
#print 'RSDMassRMS = ',myflow.GetConvergenceHistory('RSDMassRMS')

print '...end of python script.'
