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

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from adflow import ADFLOW' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pyADflow import ADFLOW
# ###################################################################

# ****************************************************************************
printHeader('Test 4: MDO tutorial -- Random Mesh -- Euler -- Scalar JST')

# ****************************************************************************
gridFile = '../inputFiles/mdo_tutorial_random_euler_scalar_jst.cgns'

options = copy.copy(adflowDefOpts)
options.update(
    {'gridfile': gridFile,
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

# Setup aeroproblem, cfdsolver, mesh and geometry.
ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, P=20000.0, T=220.0,
                 areaRef=45.5, chordRef=3.25, beta=0.0, 
                 xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=defaultFuncList)


if __name__ == "__main__": 

    solve = True
    if 'solve' not in sys.argv:
        options['restartfile'] = gridFile
        solve = False
    # Create the solver
    CFDSolver = ADFLOW(options=options, debug=False)

    standardTest(CFDSolver, ap, solve)
