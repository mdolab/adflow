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
printHeader('Test 16: Euler Convergenent Nozzle -- Flow Properties Integration ')

# ****************************************************************************
gridFile = '../inputFiles/euler_conv_nozzle.cgns'

options = copy.copy(adflowDefOpts)

options = {
    'gridFile':gridFile,
    'equationType':'euler',
    'smoother':'dadi',
    'liftIndex':2,
    'CFL':3.,
    'CFLCoarse':1.5,
    'MGCycle':'2w',
    'MGStartLevel':2,
    'nCyclesCoarse':500,
    'nCycles':2500,
    'monitorvariables':['resrho','cl','cd', 'yplus'],
    'nsubiterturb':3,
    'useNKSolver':True,
    'NKSubSpaceSize':60,
    'L2Convergence':1e-14,
    'L2ConvergenceCoarse':1e-2,
    'NKSwitchTol':1e-2,
    'nkadpc': False, 
    'vis4':0.006,
    'vis2': 0.0, 
    'blocksplitting': True, 
    'solutionPrecision':'double', 
    'flowtype':'internal'
}

ap = AeroProblem(name='conv_nozzle', alpha=30.0,  mach=0.25, T=500, P=79326.7,
                 areaRef=1., chordRef=2.,
                 evalFuncs=['mdot', 'mdot_up', 'mdot_down',
                            'mavgptot_up', 'mavgptot_down',
                            'mavgttot_up', 'mavgttot_down',
                            'mavgps_up', 'mavgps_down', 
                            'mavgmn_up', 'mavgmn_down', 
                            'thrust', 
                            'pk_up', 'pk_down', 
                            'distortionmn', 'distortionptot'
                            ], )

# Creat the solver

CFDSolver = ADFLOW(options=options, debug=False)

CFDSolver.addFamilyGroup('upstream',['INFLOW'])
CFDSolver.addFamilyGroup('downstream',['OUTFLOW'])
CFDSolver.addFamilyGroup('all_flow',['INFLOW', 'OUTFLOW'])
CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")

CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")

CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")

CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")

CFDSolver.addFunction('mavgmn', 'downstream', name="mavgmn_down")
CFDSolver.addFunction('mavgmn', 'upstream', name="mavgmn_up")

CFDSolver.addFunction('pk', 'downstream', name="pk_down")
CFDSolver.addFunction('pk', 'upstream', name="pk_up")

CFDSolver.addFunction('sigmamn', 'upstream', name="distortionmn")
CFDSolver.addFunction('sigmaptot', 'upstream', name="distortionptot")

CFDSolver.addFunction('drag', 'all_flow', name="thrust") # this naming makes it seem like wishful thinking

CFDSolver(ap)


# Check the residual
res = CFDSolver.getResidual(ap)
totalR0, totalRStart, totalRFinal = CFDSolver.getResNorms()
res /= totalR0

parPrint('Norm of residual')
reg_par_write_norm(res, 1e-10, 1e-10)

# Get and check the states
parPrint('Norm of state vector')
reg_par_write_norm(CFDSolver.getStates(), 1e-10, 1e-10)


funcs = {}
CFDSolver.evalFunctions(ap, funcs)
if MPI.COMM_WORLD.rank == 0:
    print 'Eval Functions:'
    reg_write_dict(funcs, 1e-10, 1e-10)



