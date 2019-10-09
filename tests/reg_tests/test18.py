from __future__ import print_function
import unittest
import numpy
from baseclasses import BaseRegTest

class RegTest18(unittest.TestCase):
    '''
    Test 18:
    '''
    N_PROCS = 1

    def setUp(self):
        self.ref_file = 'reg_tests/ref/test18.ref'

    def train(self):
        with BaseRegTest(self.ref_file, train=True) as handler:
            self.regression_test(handler)

    def test(self):
        with BaseRegTest(self.ref_file, train=False) as handler:
            self.regression_test(handler)

    def regression_test(self, handler, solve=False):
        '''
        This is where the actual testing happens.
        '''
        import copy
        from baseclasses import AeroProblem
        from reg_tests.commonUtils import adflowDefOpts
        from adflow import ADFLOW
        gridFile = 'input_files/conic_conv_nozzle.cgns'

        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'equationType':'euler',
            'smoother':'dadi',
            'nsubiter':3,
            'CFL':4.0,
            'CFLCoarse':1.25,
            'MGCycle':'sg',
            'MGStartLevel':-1,
            'nCyclesCoarse':250,
            'nCycles':1000,
            'nkcfl0':1e10,
            'monitorvariables':['cpu', 'resrho','cl','cd'],
            'volumevariables':['blank'],
            'surfacevariables':['mach', 'cp', 'vx', 'vy','vz', 'blank'],
            'useNKSolver':True,
            'nkswitchtol':.01,
            'nkadpc':True,
            'nkjacobianlag':5,
            'nkouterpreconits':3,
            'nkinnerpreconits':2,
            'L2Convergence':1e-10,
            'L2ConvergenceCoarse':1e-4,
            'adjointl2convergence':1e-6,
            'forcesAsTractions':True,
            'debugzipper':True,
            'nearwalldist':.001,
            'nkls':'none',
            'solutionprecision':'double',
            'adjointsubspacesize':200,
            'outerpreconits':3,
            'zipperSurfaceFamily':'output_fam',
            'flowtype':'internal',
        })

        if not solve:
            options['restartfile'] = gridFile

        # Setup aeroproblem
        ap = AeroProblem(name='nozzle', alpha=90.0, mach=0.5, altitude=0,
                areaRef=1.0, chordRef=1.0, R=287.87,
                evalFuncs=['mdot_up', 'mdot_down', #'mdot_plane',
                        'mavgptot_up', 'mavgptot_down',# 'mavgptot_plane',
                        'aavgptot_up', 'aavgptot_down',# 'aavgptot_plane',
                        'mavgttot_up', 'mavgttot_down',# 'mavgttot_plane',
                        'mavgps_up', 'mavgps_down', #'mavgps_plane'
                        'aavgps_up', 'aavgps_down', #'aavgps_plane'
                        ])

        ap.setBCVar('Pressure',  79326.7, 'downstream')
        ap.addDV('Pressure', family='downstream')

        ap.setBCVar('PressureStagnation',  100000.0, 'upstream')
        ap.addDV('PressureStagnation', family='upstream')

        ap.setBCVar('TemperatureStagnation',  500.0, 'upstream')
        ap.addDV('TemperatureStagnation', family='upstream')

        # Create the solver
        CFDSolver = ADFLOW(options=options)

        CFDSolver.addFamilyGroup('upstream',['inlet'])
        CFDSolver.addFamilyGroup('downstream',['outlet'])
        CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
        CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

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

        CFDSolver.setOption('ncycles',1000)

        # Run test
        CFDSolver(ap)

        # Check the residual
        res = CFDSolver.getResidual(ap)
        handler.par_add_norm(res, 1e-10, 1e-10)

        # Get and check the states
        handler.par_add_norm(CFDSolver.getStates(), 1e-10, 1e-10)

        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        handler.root_add_dict(funcs, 1e-10, 1e-10)