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

from reg_aeroproblems import ap_tutorial_wing 
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
    ref_file = 'solve_cl.json'
    def setUp(self):
        super().setUp()

        # self.ref_file = os.path.join(refDir, 'ref10.json')      

        
        gridFile = os.path.join(inputDir, 'mdo_tutorial_euler_scalar_jst.cgns')


        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'outputdirectory':outputDir,

            'solutionprecision':'double',
            'gridprecision':'double',

            'mgcycle':'2w',
            'ncyclescoarse':250,
            'ncycles':500,
            'usenksolver':True,
            'nkswitchtol':1e-2,

            'l2convergence':1e-14,
            'l2convergencecoarse':1e-2,
        })

        # Setup aeroproblem
        self.ap = copy.copy(ap_tutorial_wing)

        # add the default dvs to the problem 
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)


    def test_solve(self):


        self.CFDSolver.solveCL(self.ap, 0.475, alpha0=1.20, delta=0.025, tol=1e-4, autoReset=False)
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=['cl'])

        self.handler.root_add_val(funcs['mdo_tutorial_cl'] - 0.475, 'CL-CL*', rel_tol=1e-4, abs_tol=1e-4)

if __name__ == '__main__':
    unittest.main()