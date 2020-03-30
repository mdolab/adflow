from __future__ import print_function
import unittest
import numpy
from baseclasses import BaseRegTest

class RegTest16(unittest.TestCase):
    '''
    Test 16: Euler Convergenent Nozzle -- Flow Properties Integration
    '''
    N_PROCS = 4

    def setUp(self):
        self.ref_file = 'reg_tests/ref/test16.ref'

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
        gridFile = 'input_files/euler_conv_nozzle.cgns'

        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
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
        })

        # Setup aeroproblem
        ap = AeroProblem(name='conv_nozzle', alpha=00.0,  mach=0.25, T=500, P=79326.7,
                areaRef=1., chordRef=2., R=287.87,
                evalFuncs=['mdot', 'mdot_up', 'mdot_down',
                            'mavgptot_up', 'mavgptot_down',
                            'aavgptot_up', 'aavgptot_down',
                            'mavgttot_up', 'mavgttot_down',
                            'mavgps_up', 'mavgps_down',
                            'aavgps_up', 'aavgps_down',
                            'mavgmn_up', 'mavgmn_down',
                            'thrust', 'thrust_pressure',
                            'thrust_viscous', 'thrust_momentum'])

        # Create the solver
        CFDSolver = ADFLOW(options=options)
        CFDSolver.addFamilyGroup('upstream',['INFLOW'])
        CFDSolver.addFamilyGroup('downstream',['OUTFLOW'])
        CFDSolver.addFamilyGroup('all_flow',['INFLOW', 'OUTFLOW'])
        CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
        CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")

        CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
        CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")

        CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
        CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")

        CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
        CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")

        CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
        CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")

        CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
        CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")

        CFDSolver.addFunction('mavgmn', 'downstream', name="mavgmn_down")
        CFDSolver.addFunction('mavgmn', 'upstream', name="mavgmn_up")

        CFDSolver.addFunction('drag', 'all_flow', name="thrust") # this naming makes it seem like wishful thinking

        CFDSolver.addFunction('dragpressure', 'all_flow', name="thrust_pressure")
        CFDSolver.addFunction('dragviscous', 'all_flow', name="thrust_viscous")
        CFDSolver.addFunction('dragmomentum', 'all_flow', name="thrust_momentum")

        # Run test
        CFDSolver(ap)

        # Check the residual
        res = CFDSolver.getResidual(ap)
        totalR0, totalRStart, totalRFinal = CFDSolver.getResNorms()
        res /= totalR0

        handler.par_add_norm(res, 1e-10, 1e-10)

        # Get and check the states
        handler.par_add_norm(CFDSolver.getStates(), 1e-10, 1e-10)

        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        handler.root_add_dict(funcs, 1e-10, 1e-10)