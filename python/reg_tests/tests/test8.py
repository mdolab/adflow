############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE SUMB
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
# SCRIPTS ONLY. Use 'from sumb import SUMB' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pySUmb import SUMB
# ###################################################################

# ****************************************************************************
printHeader('Test 8: CRM WBT -- Euler -- Scalar JST')

# ****************************************************************************
gridFile = '../inputFiles/CRM_wbt_scalar_jst.cgns'

options = copy.copy(sumbDefOpts)
options.update(
    {'gridfile': gridFile,
     'mgcycle':'sg',
     'cfl':1.5,
     'cflcoarse':1.25,
     'resaveraging':'noresaveraging',
     'ncycles':1000,
     'monitorvariables':['resrho','cl','cd','cmy','yplus','totalr'],
     'usenksolver':True,
     'l2convergence':1e-14,
     'l2convergencecoarse':1e-4,
     'nkswitchtol':1e-1,
     'adjointl2convergence': 1e-14,
     'liftindex':3,
     'solutionprecision':'double',
     'gridprecision':'double',
 })

solve = True
if 'solve' not in sys.argv:
    options['restartfile'] = gridFile
    solve = False
# Create the solver
CFDSolver = SUMB(options=options, debug=False)

# Setup aeroproblem, cfdsolver, mesh and geometry.
ap = AeroProblem(name='CRM', alpha=1.8, mach=0.80, P=20000.0, T=220.0,
                 areaRef=45.5, chordRef=3.25, beta=0.0, 
                 xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=defaultFuncList)

standardTest(CFDSolver, ap, solve)

