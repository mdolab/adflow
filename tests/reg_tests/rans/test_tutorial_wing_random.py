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

from reg_aeroproblems import ap_tutorial_wing 
from reg_test_classes import SolveRegTest, FunctionalsRegTest, AdjointRegTest
from reg_default_options import defaultFuncList


refDir, inputDir, outputDir = BaseRegTest.getLocalDirPaths(baseDir)



@unittest.skip('no old reg file')
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

        self.ref_file = os.path.join(refDir, 'ref12.json')      

        # create the object used to compare the values to the references 
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)
        
        gridFile = os.path.join(inputDir, 'mdo_tutorial_random_rans_scalar_jst.cgns')


        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'outputdirectory':outputDir,

            'mgcycle':'sg',
            'equationtype':'RANS',
            'smoother':'dadi',
            'cfl':1.5,
            'cflcoarse':1.25,
            'resaveraging':'noresaveraging',
            'nsubiter':3,
            'nsubiterturb':3,
            'ncyclescoarse':100,
            'ncycles':1000,
            'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cmz','yplus','totalr'],
            'usenksolver':True,
            'l2convergence':1e-14,
            'l2convergencecoarse':1e-4,
            'nkswitchtol':1e-3,
            'adjointl2convergence': 1e-14,
            'frozenturbulence':False,

            # 'solutionprecision':'double',
            # 'gridprecision':'double',

            # 'mgcycle':'2w',
            # 'ncyclescoarse':250,
            # 'ncycles':500,
            # 'usenksolver':True,
            # 'nkswitchtol':1e-2,

            # 'l2convergence':1e-14,
            # 'l2convergencecoarse':1e-2,
        })

        # Setup aeroproblem
        self.ap = copy.copy(ap_tutorial_wing)

        # add the default dvs to the problem 
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

@unittest.skip('Broken')
class TestFunctionals(FunctionalsRegTest, unittest.TestCase):
    '''
    Tests that given a flow state the residuals, function, forces/tractions, 
    and jacobian vector products are accurate.

    Test 1: MDO tutorial -- Euler -- Scalar JST
    '''
    N_PROCS = 4


    def setUp(self, train=False):

        self.ref_file = os.path.join(refDir, 'ref7.json')
        
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)


        gridFile = os.path.join(inputDir, 'mdo_tutorial_random_rans_scalar_jst.cgns')

        options = copy.copy(adflowDefOpts)
        options.update({
            'outputdirectory':outputDir,
            'gridfile': gridFile,
            'restartfile': gridFile,
            # 'solutionprecision':'double',
            # 'gridprecision':'double',
            
            # 'mgcycle':'2w',
            
            # 'useblockettes': False,

            'mgcycle':'sg',
            'equationtype':'RANS',
            'smoother':'dadi',
            'cfl':1.5,
            'cflcoarse':1.25,
            'resaveraging':'noresaveraging',
            'nsubiter':3,
            'nsubiterturb':3,
            'ncyclescoarse':100,
            'ncycles':1000,
            'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cmz','yplus','totalr'],
            'usenksolver':True,
            'l2convergence':1e-14,
            'l2convergencecoarse':1e-4,
            'nkswitchtol':1e-3,
            'adjointl2convergence': 1e-14,
            'frozenturbulence':False,
        })


        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=True)


        # Setup aeroproblem
        self.ap = copy.copy(ap_tutorial_wing)

        # add the default dvs to the problem 
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code 
        self.CFDSolver.getResidual(self.ap)


@unittest.skip('no old reg file')
class TestAdjoint(AdjointRegTest, unittest.TestCase):
    '''
    Tests that given a flow state the residuals, function, forces/tractions, 
    and jacobian vector products are accurate.

    Test 1: MDO tutorial -- Euler -- Scalar JST
    '''
    N_PROCS = 4


    def setUp(self, train=False):

        self.ref_file = os.path.join(refDir, 'ref14.json')
        
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)


        gridFile = os.path.join(inputDir, 'mdo_tutorial_random_rans_scalar_jst.cgns')
        ffdFile = os.path.join(inputDir, 'mdo_tutorial_ffd.fmt')

        options = copy.copy(adflowDefOpts)

        options.update({
                    'outputdirectory':outputDir,
                    'gridfile': gridFile,
                    'restartfile': gridFile,
        #             'solutionprecision':'double',
        #             'gridprecision':'double',
                    
        #             'mgcycle':'2w',

        #             'adjointl2convergence': 1e-14,
                    
        #             'useblockettes': False,
                    'mgcycle':'2w',
                    'equationtype':'RANS',
                    'smoother':'dadi',
                    'cfl':1.5,
                    'cflcoarse':1.25,
                    'resaveraging':'noresaveraging',
                    'nsubiter':3,
                    'nsubiterturb':3,
                    'ncyclescoarse':100,
                    'ncycles':1000,
                    'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cmz','yplus','totalr'],
                    'usenksolver':True,
                    'l2convergence':1e-14,
                    'l2convergencecoarse':1e-4,
                    'nkswitchtol':1e-3,
                    'adjointl2convergence': 1e-14,
                    'frozenturbulence':False,
                    'blockSplitting':False,
        })



        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({
            'gridFile': gridFile,
        })
        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=True)


        # Setup aeroproblem
        self.ap = copy.copy(ap_tutorial_wing)
        self.ap.evalFuncs = ['fx', 'mz']

        # add the default dvs to the problem 
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)


        # Setup geometry/mesh
        DVGeo = DVGeometry(ffdFile)
        nTwist = 6
        DVGeo.addRefAxis('wing', Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist),
                                    y=numpy.zeros(nTwist),
                                    z=numpy.linspace(0,14, nTwist), k=2))
        def twist(val, geo):
            for i in range(nTwist):
                geo.rot_z['wing'].coef[i] = val[i]

        DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
        DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
        mesh = USMesh(options=mesh_options)
        self.CFDSolver.setMesh(mesh)
        self.CFDSolver.setDVGeo(DVGeo)

        # propagates the values from the restart file throughout the code 
        self.CFDSolver.getResidual(self.ap)


if __name__ == '__main__':
    unittest.main()