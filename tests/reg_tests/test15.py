from __future__ import print_function
import unittest
import numpy
from baseclasses import BaseRegTest

class RegTest15(unittest.TestCase):
    '''
    Test 15: NACA 0012 2D Time-Accurate, Forced motion, Rigid Rotation of Mesh - DADI Smoother
    '''
    N_PROCS = 4

    def setUp(self):
        self.ref_file = 'reg_tests/ref/test15.ref'

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
        import copy, os
        from baseclasses import AeroProblem
        from reg_tests.commonUtils import adflowDefOpts, defaultFuncList
        from adflow import ADFLOW
        gridFile = 'input_files/naca0012_rans-L2.cgns'

        options = copy.copy(adflowDefOpts)
        k = 0.0808
        M = 0.6
        gamma = 1.4
        R = 287.085
        T = 280.0
        c = 1.0
        alpha_m = 2.77 # 2.89 #2.77 #Modified numbers
        alpha_0 = 2.34 # 2.41 #2.34

        omega = 2*M*numpy.sqrt(gamma*R*T)*k/c
        deltaAlpha = -alpha_0*numpy.pi/180.0

        # Set forcing frequency and other information
        f = 10.0 # [Hz] Forcing frequency of the flow
        period = 1.0/f # [sec]
        nStepPerPeriod = 8
        nPeriods = 1
        nfineSteps = nStepPerPeriod*nPeriods
        dt = period / nStepPerPeriod # [s] The actual timestep

        options.update({
            'gridfile': gridFile,
            'writevolumesolution':False,
            'vis4':.025,
            'vis2':0.5,
            'restrictionrelaxation':.5,
            'smoother':'dadi',
            'equationtype':'RANS',
            'equationmode':'unsteady',
            'timeIntegrationscheme':'bdf',
            'ntimestepsfine':nfineSteps,
            'deltat':dt,
            'nsubiterturb':10,
            'nsubiter':5,
            'useale':False,
            'usegridmotion':True,
            'cfl':2.5,
            'cflcoarse':1.2,
            'ncycles':2000,
            'mgcycle':'3w',
            'mgstartlevel':1,
            'monitorvariables':['cpu','resrho','cl','cd','cmz'],
            'usenksolver':False,
            'l2convergence':1e-6,
            'l2convergencecoarse':1e-4,
            'qmode':True,
            'alphafollowing': False,
            'blockSplitting':True,
            'useblockettes':False,
        })

        # Setup aeroproblem
        ap = AeroProblem(name='0012pitching', alpha=alpha_m,  mach=M, machRef=M,
                reynolds=4800000.0,reynoldsLength=c, T=T, R=R, areaRef=1.0,
                chordRef=c, evalFuncs=['cl','cd','cmz'], xRef=0.25, xRot=0.25,
                degreePol=0, coefPol=[0.0], degreeFourier=1, omegaFourier=omega,
                cosCoefFourier=[0.0,0.0], sinCoefFourier=[deltaAlpha])

        # Create the solver
        CFDSolver = ADFLOW(options=options)
        CFDSolver.addSlices('z',[0.5])

        # Run test
        CFDSolver(ap)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        CFDSolver.checkSolutionFailure(ap, funcs)
        handler.root_add_dict(funcs, 1e-6, 1e-6)