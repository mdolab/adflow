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
printHeader('Test 16: Euler Convergenent Nozzle -- Flow Properties Integration ')

# ****************************************************************************
gridFile = '../inputFiles/euler_conv_nozzle.cgns'

options = copy.copy(sumbDefOpts)

options = {
    # Common Parameters
    'gridFile':gridFile,
    # Physics Parameters
    'equationType':'euler',
    'smoother':'dadi',
    'liftIndex':2,
    # Common Parameters
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
    # Convergence Parameters
    'L2Convergence':1e-12,
    'L2ConvergenceCoarse':1e-2,
    #'miniterationnum':100,
    'NKSwitchTol':1e-2,
    'nkadpc': False, 
    #artifical viscosity
    'vis4':0.006,
    'vis2': 0.0
}

ap = AeroProblem(name='conv_nozzle', alpha=0.0,  mach=0.25, T=500, P=100000,
                 areaRef=1., chordRef=2.,
                 evalFuncs=['mdot', 'mdot_up', 'mdot_down',
                            'mavgptot_up', 'mavgptot_down',
                            'mavgttot_up', 'mavgttot_down',
                            'mavgps_up', 'mavgps_down'])

solve = True
if 'solve' not in sys.argv:
    options['restartfile'] = gridFile
    solve = False
# Creat the solver

CFDSolver = SUMB(options=options, debug=False)

CFDSolver.addFamilyGroup('upstream',['INFLOW'])
CFDSolver.addFamilyGroup('downstream',['OUTFLOW'])
CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")

if solve:
    # We are told that we must first solve the problem, most likely
    # for a training run. 
    CFDSolver(ap)

funcs = {}
CFDSolver.evalFunctions(ap, funcs)
if MPI.COMM_WORLD.rank == 0:
    reg_write_dict(funcs, 1e-10, 1e-10)
    # print 'Functions Sens:', funcsSens
