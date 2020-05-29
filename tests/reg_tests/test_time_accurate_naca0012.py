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
from baseclasses import BaseRegTest, AeroProblem

# get the directories that we will use to import some packages in this repo 
baseDir = os.path.dirname(os.path.abspath(__file__))

BaseRegTest.setLocalPaths(baseDir, sys.path)

# import the ADFLOW module this test lives in 
from python.pyADflow import ADFLOW

# import the testing utilities that live a few directories up 
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts

from reg_aeroproblems import ap_naca0012_time_acc 
from reg_test_classes import SolveRegTest, FunctionalsRegTest, AdjointRegTest
from reg_default_options import defaultFuncList


refDir, inputDir, outputDir = BaseRegTest.getLocalDirPaths(baseDir)


class TestSolve(SolveRegTest, unittest.TestCase):
    '''
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 

    Test 15: NACA 0012 2D Time-Accurate, Forced motion, Rigid Rotation of Mesh - DADI Smoother
    '''
    N_PROCS = 4

    def setUp(self):

        self.ref_file = os.path.join(refDir, 'ref15.json')      

        # create the object used to compare the values to the references 
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)
        
        gridFile = os.path.join(inputDir, 'naca0012_rans-L2.cgns.cgns')
        print(gridFile)
        gridFile = '../../python/inputFiles/naca0012_rans-L2.cgns'
        # print(gridFile)
        # quit()
        # import ipdb; ipdb.set_trace()
        f = 10.0 # [Hz] Forcing frequency of the flow
        period = 1.0/f # [sec]
        nStepPerPeriod = 8
        nPeriods = 1
        nfineSteps = nStepPerPeriod*nPeriods
        dt = period / nStepPerPeriod # [s] The actual timestep

        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'outputdirectory':outputDir,

            # 'mgcycle':'sg',
            # 'equationtype':'RANS',
            # 'smoother':'dadi',
            # 'cfl':1.5,
            # 'cflcoarse':1.25,
            # 'resaveraging':'noresaveraging',
            # 'nsubiter':3,
            # 'nsubiterturb':3,
            # 'ncyclescoarse':100,
            # 'ncycles':1000,
            # 'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cmz','yplus','totalr'],
            # 'usenksolver':True,
            # 'l2convergence':1e-14,
            # 'l2convergencecoarse':1e-4,
            # 'nkswitchtol':1e-3,
            # 'adjointl2convergence': 1e-14,
            # 'frozenturbulence':False,

            # 'solutionprecision':'double',
            # 'gridprecision':'double',

            # 'mgcycle':'2w',
            # 'ncyclescoarse':250,
            # 'ncycles':500,
            # 'usenksolver':True,
            # 'nkswitchtol':1e-2,

            # 'l2convergence':1e-14,
            # 'l2convergencecoarse':1e-2,

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
            }
        )

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca0012_time_acc)

        # add the default dvs to the problem 


        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        self.CFDSolver.addSlices('z',[0.5])
      

    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, rtol=1e-8, atol=1e-8)



if __name__ == '__main__':
    unittest.main()