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
from adflow import ADFLOW

# ****************************************************************************
printHeader('Test 3: MDO tutorial -- Euler -- Matrix')

# ****************************************************************************
gridFile = '../inputFiles/mdo_tutorial_euler_upwind.cgns'

options = copy.copy(adflowDefOpts)
options.update(
    {'gridfile': gridFile,
     'mgcycle':'2w',
     'ncyclescoarse':250,
     'discretization':'upwind',
     'coarsediscretization':'upwind',
     'ncycles':500,
     'monitorvariables':['cpu', 'resrho','cl','cd','cmz','totalr'],
     'usenksolver':True,
     'vis4':0.1,
     'l2convergence':1e-14,
     'l2convergencecoarse':1e-2,
     'nkswitchtol':1e-2,
     'adjointl2convergence': 1e-14,
     'solutionprecision':'double',
     'gridprecision':'double',
 })


# Setup aeroproblem, cfdsolver, mesh and geometry.
ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, P=20000.0, T=220.0,
                 areaRef=45.5, chordRef=3.25, beta=0.0, R=287.87,
                 xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=defaultFuncList)

def setup_cb(comm):

    # Create the solver
    CFDSolver = ADFLOW(options=options, debug=True)

    return CFDSolver, None, None, None

if __name__ == "__main__":

    solve = True
    if 'solve' not in sys.argv:
        options['restartfile'] = gridFile
        solve = False

    CFDSolver, _, _, _ = setup_cb(MPI.COMM_WORLD)

    standardTest(CFDSolver, ap, solve)
