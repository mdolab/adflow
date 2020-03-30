from __future__ import print_function
import unittest
import numpy
from baseclasses import BaseRegTest

class RegTest10(unittest.TestCase):
    '''
    Test 10: MDO tutorial -- Solve-CL Test
    '''
    N_PROCS = 4

    def setUp(self):
        self.ref_file = 'reg_tests/ref/test10.ref'

    def train(self):
        with BaseRegTest(self.ref_file, train=True) as handler:
            self.regression_test(handler)

    def test(self):
        with BaseRegTest(self.ref_file, train=False) as handler:
            self.regression_test(handler)

    def regression_test(self, handler):
        '''
        This is where the actual testing happens.
        '''
        import copy
        from baseclasses import AeroProblem
        from reg_tests.commonUtils import adflowDefOpts
        from adflow import ADFLOW
        gridFile = 'input_files/mdo_tutorial_euler_scalar_jst.cgns'

        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
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

        # Setup aeroproblem, cfdsolver, mesh and geometry.
        ap = AeroProblem(name='mdo_tutorial', alpha=1.20, mach=0.80,
                         altitude=10000.0, areaRef=45.5, chordRef=3.25)

        # Create the solver
        CFDSolver = ADFLOW(options=options, debug=True)

        # Run CL solve
        CFDSolver.solveCL(ap, 0.475, alpha0=1.20, delta=0.025, tol=1e-4, autoReset=False)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, evalFuncs=['cl'])
        handler.root_add_val(funcs['mdo_tutorial_cl'] - 0.475, 1e-4, 1e-4)

