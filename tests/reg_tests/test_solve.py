# built-ins 
import unittest
import numpy
import os
import sys
import copy
from parameterized import parameterized_class
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
# from pyADflow import ADFLOW
from adflow import ADFLOW

# import the testing utilities that live a few directories up 
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts

from reg_aeroproblems import ap_tutorial_wing 
# from reg_test_classes import SolveRegTest, FunctionalsRegTest, AdjointRegTest


refDir, inputDir, outputDir = BaseRegTest.getLocalDirPaths(baseDir)



class test_objects():
    class RegTest(unittest.TestCase):
        def train(self):
            self.handler.setMode(train=True)
            self.handler.setRef({})

            # get all of the testing methods
            # all of the tests written in this framework must start with "test_"
            tests = [x for x in dir(self) if x.startswith('test_')]
            for test in tests:
                test_func = getattr(self, test)
                test_func()

            trained_ref = self.handler.getRef()
    
            utils.writeRefToJson(self.ref_file, trained_ref)



refDir = '/home/josh/repos/MACH/adflow_clean/tests/reg_tests/refs'
@parameterized_class([
    # scalar JST
    { "options": {
        'gridfile':os.path.join(inputDir, 'mdo_tutorial_euler_scalar_jst.cgns'),
        # 'restartfile':os.path.join(inputDir, 'mdo_tutorial_euler_scalar_jst.cgns'),
        'l2convergence': 1e-14,
        'mgcycle':'2w',
        'ncyclescoarse':250,
        'usenksolver':True,
        'useblockettes': False,
        'nkswitchtol':1e-2,

          },
        "ref_file": 'solve_euler_scalar_jst.json', 
        'ap':ap_tutorial_wing,
    }, 
    # # Matrix JST
    # { "options": {
    #     'gridfile':os.path.join(inputDir, 'mdo_tutorial_euler_matrix.cgns'),
    #     'restartfile':os.path.join(inputDir, 'mdo_tutorial_euler_matrix.cgns'),

    #     'mgcycle':'2w',
    #     'ncyclescoarse':250,
    #     'usenksolver':True,
    #     'nkswitchtol':1e-2,

    #     'vis4':0.1,
    #     'discretization':'central plus matrix dissipation',
    #     'coarsediscretization':'central plus matrix dissipation',
    #     'l2convergence':1e-14,
    #     'useblockettes': False,

    #     },
    #     "ref_file":'euler_matrix_jst.json', 
    #     'ap':ap_tutorial_wing,
    # },  
    # # Upwind
    # { "options": {
    #     'gridfile': os.path.join(inputDir, 'mdo_tutorial_euler_upwind.cgns'),
    #     'restartfile': os.path.join(inputDir, 'mdo_tutorial_euler_upwind.cgns'),
       
    #     'outputdirectory':outputDir,
    #     'mgcycle':'2w',
    #     'ncyclescoarse':250,
    #     'usenksolver':True,
    #     'nkswitchtol':1e-2,
    #     'vis4':0.1,
    #     'discretization':'upwind',
    #     'useblockettes': False,
    #     'l2convergence':1e-14
    #     },
    #     "ref_file":'euler_upwind.json',
    #     'ap':ap_tutorial_wing
    # }
])
class TestSolve(test_objects.RegTest):
# class TestSolve(SolveRegTest):
    '''
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 

    Test 9: MDO tutorial -- Euler -- Solution Test
    '''
    N_PROCS = 4

    def setUp(self):
        
        print(self.ref_file)
        ref_file = os.path.join(refDir, self.ref_file)      

        # create the object used to compare the values to the references 
        ref = utils.readJSONRef(ref_file)
        self.handler = BaseRegTest(ref, train=False)
        

            
        options = copy.copy(adflowDefOpts)
        options.update(self.options)


        # Setup aeroproblem
        # self.ap = copy.copy(ap_tutorial_wing)

        # # add the default dvs to the problem 
        # for dv in defaultAeroDVs:
        #     self.ap.addDV(dv)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
    
    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap)
        utils.assert_states_allclose(self.handler, self.CFDSolver)
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap)


if __name__ == '__main__':
    unittest.main()