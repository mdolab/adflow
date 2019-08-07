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
    'flowtype':'internal',
}

ap = AeroProblem(name='conv_nozzle', alpha=00.0,  mach=0.25, T=500, P=79326.7,
                 areaRef=1., chordRef=2., R=287.87,
                 evalFuncs=['mdot', 'mdot_up', 'mdot_down',
                            'mavgptot_up', 'mavgptot_down',
                            'aavgptot_up', 'aavgptot_down',
                            'mavgttot_up', 'mavgttot_down',
                            'mavgps_up', 'mavgps_down',
                            'aavgps_up', 'aavgps_down',
                            'mavgmn_up', 'mavgmn_down',
                            'thrust',
                            'thrust_pressure', 'thrust_viscous', 'thrust_momentum'
                            ], )


def setup_cb(comm):
    solver = ADFLOW(options=options, comm=comm, debug=False)

    solver.addFamilyGroup('upstream',['INFLOW'])
    solver.addFamilyGroup('downstream',['OUTFLOW'])
    solver.addFamilyGroup('all_flow',['INFLOW', 'OUTFLOW'])
    solver.addFunction('mdot', 'upstream', name="mdot_up")
    solver.addFunction('mdot', 'downstream', name="mdot_down")

    solver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
    solver.addFunction('mavgptot', 'upstream', name="mavgptot_up")

    solver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
    solver.addFunction('aavgptot', 'upstream', name="aavgptot_up")

    solver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
    solver.addFunction('mavgttot', 'upstream', name="mavgttot_up")

    solver.addFunction('mavgps', 'downstream', name="mavgps_down")
    solver.addFunction('mavgps', 'upstream', name="mavgps_up")

    solver.addFunction('aavgps', 'downstream', name="aavgps_down")
    solver.addFunction('aavgps', 'upstream', name="aavgps_up")

    solver.addFunction('mavgmn', 'downstream', name="mavgmn_down")
    solver.addFunction('mavgmn', 'upstream', name="mavgmn_up")

    solver.addFunction('drag', 'all_flow', name="thrust") # this naming makes it seem like wishful thinking

    solver.addFunction('dragpressure', 'all_flow', name="thrust_pressure")
    solver.addFunction('dragviscous', 'all_flow', name="thrust_viscous")
    solver.addFunction('dragmomentum', 'all_flow', name="thrust_momentum")

    return solver, None, None, None

if __name__ == "__main__":

    # Creat the solver

    CFDSolver, _, _, _ = setup_cb(MPI.COMM_WORLD)

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
        print('Eval Functions:')
        reg_write_dict(funcs, 1e-10, 1e-10)
