# built-ins
import unittest
import numpy
import os
import sys
import copy

# MACH testing class
from baseclasses import BaseRegTest

# get the directories that we will use to import some packages in this repo 
baseDir = os.path.dirname(os.path.abspath(__file__))

BaseRegTest.setLocalPaths(baseDir, sys.path)

# import the ADFLOW module this test lives in 
from python.pyADflow import ADFLOW

# import the testing utilities that live a few directories up 
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs

from reg_aeroproblems import ap_tutorial_wing_laminar 
from reg_test_classes import SolveRegTest, FunctionalsRegTest


refDir, inputDir, outputDir = BaseRegTest.getLocalDirPaths(baseDir)


@unittest.skip('This test is broken')
class TestSolve(SolveRegTest, unittest.TestCase):
    '''
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 

    Test 9: MDO tutorial -- Euler -- Solution Test
    '''
    N_PROCS = 4

    def setUp(self):

        self.ref_file = os.path.join(refDir, 'ref11.json')      

        # create the object used to compare the values to the references 
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)
        
        gridFile = os.path.join(inputDir, 'mdo_tutorial_viscous_scalar_jst.cgns')


        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'outputdirectory':outputDir,

            'solutionprecision':'double',
            'gridprecision':'double',

            'equationtype':'Laminar NS',

            'mgcycle':'2w',
            'cfl':1.5,
            'cflcoarse':1.25,
            'ncyclescoarse':250,
            'ncycles':10000,
            'usenksolver':True,
            'nkswitchtol':1e-2,

            'l2convergence':1e-14,
            'l2convergencecoarse':1e-2,
        })

        # Setup aeroproblem
        self.ap = copy.copy(ap_tutorial_wing_laminar)

        # add the default dvs to the problem 
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        
@unittest.skip('This test is broken')
class TestFunctionals(FunctionalsRegTest, unittest.TestCase):
    '''
    Tests that given a flow state the residuals, function, forces/tractions, 
    and jacobian vector products are accurate.

    Test 1: MDO tutorial -- Euler -- Scalar JST
    '''
    N_PROCS = 1


    def setUp(self, train=False):

        self.ref_file = os.path.join(refDir, 'ref5.json')
        
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)


        gridFile = os.path.join(inputDir, 'mdo_tutorial_viscous_scalar_jst.cgns')

        options = copy.copy(adflowDefOpts)
        options.update({
            'outputdirectory':outputDir,
            'gridfile': gridFile,
            'restartfile': gridFile,
            'solutionprecision':'double',
            'gridprecision':'double',
            
            'mgcycle':'2w',

            'equationtype':'Laminar NS',

            'useblockettes': False,
        })


        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=True)


        # Setup aeroproblem
        # self.ap = copy.copy(ap_tutorial_wing_laminar)
        from baseclasses import AeroProblem
        self.ap =    copy.copy(ap_tutorial_wing_laminar)
        # add the default dvs to the problem 
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code 
        self.CFDSolver.getResidual(self.ap)


if __name__ == '__main__':
    unittest.main()