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
printHeader('Test10: MDO tutorial -- Solve-CL Test')


gridFile = '../inputFiles/mdo_tutorial_euler_scalar_jst.cgns'
options = copy.copy(sumbDefOpts)
options.update(
    {'gridfile': gridFile,
     'mgcycle':'2w',
     'ncyclescoarse':250,
     'ncycles':500,
     'monitorvariables':['resrho','cl','cd','cmz','totalr'],
     'usenksolver':True,
     'l2convergence':1e-14,
     'l2convergencecoarse':1e-2,
     'nkswitchtol':1e-2,
     'adjointl2convergence': 1e-14,
     'solutionprecision':'single',
     'gridprecision':'double',
 })
ap = AeroProblem(name='mdo_tutorial', alpha=1.20, mach=0.80, altitude=10000.0,
                 areaRef=45.5, chordRef=3.25)

CFDSolver = SUMB(options=options, debug=True)
CFDSolver.solveCL(ap, 0.475, alpha0=1.20, delta=0.025, tol=1e-4, autoReset=False)
funcs = {}
CFDSolver.evalFunctions(ap, funcs, evalFuncs=['cl'])
if MPI.COMM_WORLD.rank == 0:
    print 'CL-CL*'
    reg_write(funcs['mdo_tutorial_cl'] - 0.475, 1e-4, 1e-4)
