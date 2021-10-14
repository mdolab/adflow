# Author: Sicheng He
# Date: 2019/10/02
# test adflow with flexible mesh time spectral Euler code
# Classic Davis pitching airfoil case (CT6)
# Sicheng He made this 2019/02/24
# Borrow heavily from  Eirikur Jonsson's file unsteady_run4.py (unsteady time accurate)
from __future__ import print_function
############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI
from baseclasses import AeroProblem
from mdo_regression_helper import *
from commonUtils import *
from idwarp import USMesh
import numpy

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from adflow import ADFLOW' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pyADflow import ADFLOW
# ###################################################################

# util class
class Transfer():

    # simplified transfer class
    # converting csd displacement to cfd surface nodes

    def __init__(self, alpha, xRot, aeroSolver):

        # takes in displacement history

        self.alpha = alpha
        self.ntimeintervalsspectral = len(alpha)

        self.aeroSolver = aeroSolver # a shallow copy of CFD solver

        self.xRot = xRot

    def getUndeformedSurfaceNodes(self):
        
        self.MDGroup = self.aeroSolver.allWallsGroup

        self.cfdPts0 = self.aeroSolver.getSurfaceCoordinates(self.MDGroup,includeZipper=False)

    def setDisplacements(self):

        xRot = self.xRot
        ntimeintervalsspectral = self.ntimeintervalsspectral
        alpha = self.alpha # notice a shallow copy introduced here; dont change the underlying obj!
        cfdPoints_init = self.cfdPts0 # notice a shallow copy introduced here; dont change the underlying obj!

        N_pts = cfdPoints_init.shape[0]

        self.cfdPts = []

        for sps in xrange(ntimeintervalsspectral):

            cfdPoints_deformed = numpy.zeros((N_pts, 3))

            ptch_loc = alpha[sps]

            cc = numpy.cos(ptch_loc)
            ss = numpy.sin(ptch_loc)

            for j in xrange(N_pts):

                cfdPoints_deformed[j, 0] = (  cc*(cfdPoints_init[j, 0] - xRot) + ss*cfdPoints_init[j, 1] + xRot)
                cfdPoints_deformed[j, 1] = (- ss*(cfdPoints_init[j, 0] - xRot) + cc*cfdPoints_init[j, 1])
                cfdPoints_deformed[j, 2] = cfdPoints_init[j, 2]

            self.cfdPts.append(cfdPoints_deformed)

    def setVolumeMesh(self):

        ntimeintervalsspectral = self.ntimeintervalsspectral

        for sps in xrange(ntimeintervalsspectral):

            self.aeroSolver.mesh.setSurfaceCoordinates(self.cfdPts[sps])
            self.aeroSolver.mesh.warpMesh()
            m = self.aeroSolver.mesh.getSolverGrid()
            self.aeroSolver.adflow.warping.setgridforoneinstance(m,sps=sps+1)

        self.aeroSolver._updateGeomInfo = True
        self.aeroSolver.updateGeometryInfo()


# parameters
gridFile = '../inputFiles/naca64A010_euler-L2.cgns'
outputDirectory = 'OUTPUT/'

ntimeintervalsspectral = 3

# Geometry params
areaRef = 1.0
chordRef = 1.0
xRef = 0.248
xRot = 0.248
name = '64A010pitchingTS'

# Called DI 55 in AGARD 702 report
k = 0.202
M = 0.796    
alpha_m = 0
alpha_0 = 1.01

reynolds=12.56e6
T = 305.0 # This is the average of the tunnel temp range

gamma = 1.4
R = 287.085


# Generate motion history
deltaAlpha = -alpha_0*numpy.pi/180.0
alpha = numpy.linspace(0.0, 2.0*numpy.pi, ntimeintervalsspectral + 1)[:-1]
alpha = -numpy.sin(alpha)
alpha *= deltaAlpha

# Unsteady params
# Use the reduced frequency to find the angular velocity
#omega = 2*V*k/c = 2*M*a*k/c = 2*M*sqrt(gamma*R*T)*k/c
#print k, M, R, T, chordRef, gamma

omega = 2*M*numpy.sqrt(gamma*R*T)*k/chordRef 

options = {
    'outputDirectory':outputDirectory,
    'gridfile':gridFile,
    'mgcycle': 'sg',
    'blocksplitting': True, 
    'equationmode':'time spectral',
    'l2convergence':1e-15,
    'ncycles':200000,
    'monitorvariables':['resrho', 'cl'],  
    'timeintervals': ntimeintervalsspectral,
    'alphafollowing': False,
    'usenksolver': True,
    'nkswitchtol': 1e-4, 
    'NKSubSpaceSize':400,
    'useblockettes': False,
    "applypcsubspacesize":400,
    'useanksolver': True,
    'ankswitchtol':1e-2,
    'anksubspacesize':50,
    'useexternaldynamicmesh': True,
    'usetsinterpolatedgridvelocity': True,
    'writeVolumeSolution': False,
    'writeSurfaceSolution': False,
}


ap = AeroProblem(name=name, alpha=alpha_m, mach=M, 
                 reynolds=reynolds, reynoldsLength=chordRef, T=T, R=R,
                 areaRef=areaRef, chordRef=chordRef, xRef=xRef, xRot=xRot,
                 omegaFourier=omega)

CFDSolver = ADFLOW(options=options, comm=MPI.COMM_WORLD)


meshOptions = {
      'gridFile':gridFile,
}

# mesh
mesh = USMesh(options=meshOptions)
CFDSolver.setMesh(mesh)

# deformation
TSTransfer = Transfer(alpha, xRot, CFDSolver)
TSTransfer.getUndeformedSurfaceNodes()
TSTransfer.setDisplacements()
TSTransfer.setVolumeMesh()

 
CFDSolver(ap)