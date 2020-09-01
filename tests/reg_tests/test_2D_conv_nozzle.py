# built-ins
import unittest
import numpy
import os
import sys
import copy

# MACH classes
from pygeo import DVGeometry
from pyspline import Curve
from idwarp import USMesh


# MACH testing class
from baseclasses import BaseRegTest

# get the directories that we will use to import some packages in this repo 
baseDir = os.path.dirname(os.path.abspath(__file__))

BaseRegTest.setLocalPaths(baseDir, sys.path)

# import the ADFLOW module this test lives in 
from adflow import ADFLOW

# import the testing utilities that live a few directories up 
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts

from reg_aeroproblems import ap_2D_conv_nozzle 
from reg_test_classes import test_objects

refDir, inputDir, outputDir = BaseRegTest.getLocalDirPaths(baseDir)




class TestSolve(test_objects.RegTest):
    '''
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 

    Test 9: MDO tutorial -- Euler -- Solution Test
    '''
    N_PROCS = 4
    ref_file = 'solve_2D_conv_nozzle.json'
    options = {
        'gridfile': os.path.join(inputDir, 'euler_conv_nozzle.cgns'),
        'outputdirectory':outputDir,
        'solutionprecision':'double',
        'gridprecision':'double',
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
    ap = ap_2D_conv_nozzle
    def setUp(self):
        super().setUp()

        # self.ref_file = os.path.join(refDir, 'ref16.json')      




        options = copy.copy(adflowDefOpts)
        options.update(self.options)


        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        

        self.CFDSolver.addFamilyGroup('upstream',['INFLOW'])
        self.CFDSolver.addFamilyGroup('downstream',['OUTFLOW'])
        self.CFDSolver.addFamilyGroup('all_flow',['INFLOW', 'OUTFLOW'])
        self.CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
        self.CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")

        self.CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
        self.CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")

        self.CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
        self.CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")

        self.CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
        self.CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")

        self.CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
        self.CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")

        self.CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
        self.CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")

        self.CFDSolver.addFunction('mavgmn', 'downstream', name="mavgmn_down")
        self.CFDSolver.addFunction('mavgmn', 'upstream', name="mavgmn_up")

        self.CFDSolver.addFunction('drag', 'all_flow', name="thrust") # this naming makes it seem like wishful thinking

        self.CFDSolver.addFunction('dragpressure', 'all_flow', name="thrust_pressure")
        self.CFDSolver.addFunction('dragviscous', 'all_flow', name="thrust_viscous")
        self.CFDSolver.addFunction('dragmomentum', 'all_flow', name="thrust_momentum")


    
    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap)
        utils.assert_states_allclose(self.handler, self.CFDSolver)
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap)



if __name__ == '__main__':
    unittest.main()