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
printHeader('Test9: MDO tutorial -- Euler -- Solution Test')

# ****************************************************************************
gridFile = '../inputFiles/mdo_tutorial_euler_scalar_jst.cgns'

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

def setup_cb(comm): 

    # Create the solver
    CFDSolver = ADFLOW(options=options, debug=False)
    
    return CFDSolver, None, None, None

if __name__ == "__main__": 

        
    CFDSolver, _, _, _ = setup_cb(MPI.COMM_WORLD)

    solutionTest(CFDSolver, ap)



