############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI
from baseclasses import AeroProblem
from pyspline import Curve
from pygeo import DVGeometry
from idwarp import USMesh
from mdo_regression_helper import *
from commonUtils import *
from adflow import ADFLOW

# ****************************************************************************
printHeader('Test13: MDO tutorial -- Euler -- Adjoint Test')

# ****************************************************************************
gridFile = '../inputFiles/mdo_tutorial_euler_scalar_jst.cgns'
ffdFile = '../inputFiles/mdo_tutorial_ffd.fmt'

options = copy.copy(adflowDefOpts)

options.update(
    {'gridfile': gridFile,
     'restartfile': gridFile,
     'mgcycle':'2w',
     'ncyclescoarse':250,
     'ncycles':500,
     'monitorvariables':['cpu', 'resrho','cl','cd','cmz','totalr'],
     'usenksolver':True,
     'l2convergence':1e-14,
     'l2convergencecoarse':1e-2,
     'nkswitchtol':1e-2,
     'adjointl2convergence': 1e-14,
     'solutionprecision':'double',
     'gridprecision':'double',
 })

meshOptions = copy.copy(IDWarpDefOpts)
meshOptions.update({'gridFile':gridFile})

# Setup aeroproblem, cfdsolver, mesh and geometry.
ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, P=20000.0, T=220.0,
                 areaRef=45.5, chordRef=3.25, beta=0.0, R=287.87,
                 xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=['cl', 'cd'])

def setup_cb(comm):

    # Create the solver
    CFDSolver = ADFLOW(options=options, comm=comm, debug=False)

    # Setup geometry/mesh
    DVGeo = DVGeometry(ffdFile)
    nTwist = 6
    DVGeo.addRefAxis('wing', Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist),
                                   y=numpy.zeros(nTwist),
                                   z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in range(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = USMesh(options=meshOptions, comm=comm)
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)

    return CFDSolver, mesh, DVGeo, None


if __name__ == "__main__":


    CFDSolver, mesh, DVGeo, _ = setup_cb(MPI.COMM_WORLD)

    adjointTest(CFDSolver, ap)

