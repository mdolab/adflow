from __future__ import print_function
import unittest
import numpy
from baseclasses import BaseRegTest
import os
import sys
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from python.pyADflow import ADFLOW
import copy
from baseclasses import AeroProblem
# from reg_tests.commonUtils import standard_test, adflowDefOpts, defaultFuncList
from commonUtils import *    
from runpy import run_path


class TestSolution(unittest.TestCase):
    '''
    Test 9: MDO tutorial -- Euler -- Solution Test
    '''
    N_PROCS = 4

    def setUp(self):

        self.ref_file = os.path.join(baseDir, 'ref/ref9.py')
        ref = run_path(self.ref_file)['ref']

        self.handler = BaseRegTest(ref, train=False)




        gridFile = os.path.join(baseDir, '../input_files/mdo_tutorial_euler_scalar_jst.cgns')


        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
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
        ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, P=20000.0,
                        T=220.0, areaRef=45.5, chordRef=3.25, beta=0.0,
                        xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=defaultFuncList)
        # Standard test for solving the problem.
        for dv in defaultAeroDVs:
            ap.addDV(dv)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        
        self.ap = ap


    def test_solve(self):
        self.CFDSolver(self.ap)

        assert_residuals_allclose(self.handler, self.CFDSolver, self.ap)
        assert_functions_allclose(self.handler, self.CFDSolver, self.ap)
        assert_states_allclose(self.handler, self.CFDSolver)

    def train(self):
        self.handler.train = True
        self.handler.setRef({})

        
        self.test_solve()


        trained_ref = self.handler.getRef()
        writeRef(self.ref_file, trained_ref)


class TestFunctionals(unittest.TestCase):
    '''
    Test 1: MDO tutorial -- Euler -- Scalar JST
    '''
    N_PROCS = 4


    def setUp(self, train=False):

        self.ref_file = os.path.join(baseDir, 'ref/ref1.py')
        ref = run_path(self.ref_file)['ref']

        self.handler = BaseRegTest(ref, train=False)


        gridFile = os.path.join(baseDir, '../input_files/mdo_tutorial_euler_scalar_jst.cgns')

        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'restartfile': gridFile,
            'mgcycle':'2w',
            'l2convergence':1e-14,
            'ncyclescoarse':250,
            'ncycles':500,
            'monitorvariables':['cpu', 'resrho','cl','cd','cmz','totalr'],
            'usenksolver':True,
            'l2convergencecoarse':1e-2,
            'nkswitchtol':1e-2,
            'adjointl2convergence': 1e-14,
            'solutionprecision':'double',
            'gridprecision':'double',
        })

        # Setup aeroproblem, cfdsolver, mesh and geometry.
        ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, P=20000.0, T=220.0,
                        areaRef=45.5, chordRef=3.25, beta=0.0, R=287.87,
                        xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=defaultFuncList)

        # Create the solver
        CFDSolver = ADFLOW(options=options, debug=True)


        for dv in defaultAeroDVs:
            ap.addDV(dv)

        self.ap = ap


        # propagates the values from the restart file throughout the code 
        CFDSolver.getResidual(ap)
        self.CFDSolver = CFDSolver

    def test_restart_read(self):

        assert_problem_size_equal(self.handler, self.CFDSolver)

        assert_states_allclose(self.handler, self.CFDSolver)

    def test_residuals(self):
        assert_residuals_allclose(self.handler, self.CFDSolver, self.ap)
        
    def test_functions(self):
        assert_functions_allclose(self.handler, self.CFDSolver, self.ap)

    def test_forces_and_tractions(self):
        assert_forces_allclose(self.handler, self.CFDSolver)
        assert_tractions_allclose(self.handler, self.CFDSolver)

        # Reset the option
        self.CFDSolver.setOption('forcesAsTractions', True)

        # Make sure we can write the force file.
        self.CFDSolver.writeForceFile('forces.txt')



    # ------------------- Derivative routine checks ----------------------------
    def test_jac_vec_prod_fwd(self):
        assert_fwd_mode_allclose(self.handler, self.CFDSolver)

    def test_jac_vec_prod_bwd(self):
        # assert_fwd_mode_allclose(self.handler, self.CFDSolver)
        assert_bwd_mode_allclose(self.handler, self.CFDSolver)

    def test_dot_products(self):
        assert_dot_products_allclose(self.handler, self.CFDSolver)


    def train(self):
        self.handler.train = True
        self.handler.setRef({})

        
        self.test_solve()
        self.test_restart_read()
        self.test_residuals()
        self.test_functions()
        self.test_forces_and_tractions()
        self.test_jac_vec_prod_fwd()
        self.test_jac_vec_prod_bwd()
        self.test_dot_products()

        trained_ref = self.handler.getRef()
        writeRef(self.ref_file, trained_ref)


if __name__ == '__main__':
    unittest.main()