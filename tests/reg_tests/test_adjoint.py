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

from reg_aeroproblems import ap_tutorial_wing , ap_CRM, ap_tutorial_wing_laminar
# from reg_test_classes import SolveRegTest, FunctionalsRegTest, AdjointRegTest


refDir, inputDir, outputDir = BaseRegTest.getLocalDirPaths(baseDir)



class test_objects():
    class RegTest(unittest.TestCase):
        def setUp(self):

            ref_file = os.path.join(refDir, self.ref_file)
            
            ref = utils.readJSONRef(ref_file)
            self.handler = BaseRegTest(ref)


        def train(self):

            self.handler.train = True
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
    # Tutorial scalar JST 
    { "options": {
        'gridfile':os.path.join(inputDir, 'mdo_tutorial_euler_scalar_jst.cgns'),
        'restartfile':os.path.join(inputDir, 'mdo_tutorial_euler_scalar_jst.cgns'),
        'l2convergence': 1e-14,
        'adjointl2convergence': 1e-14,

        'mgcycle':'2w',
        'ncyclescoarse':250,
        'usenksolver':True,
        'useblockettes': False,
        'frozenturbulence':False,
        'blockSplitting':False,
        },
        "ref_file": 'adjoint_euler_scalar_jst_tut_wing.json', 
        'ap':ap_tutorial_wing,
        'evalFuncs':['cl', 'cd'],
        'name': 'single_core',
        'N_PROCS': 1

    }, 
    # Tutorial scalar JST
    { "options": {
        'gridfile':os.path.join(inputDir, 'mdo_tutorial_euler_scalar_jst.cgns'),
        'restartfile':os.path.join(inputDir, 'mdo_tutorial_euler_scalar_jst.cgns'),
        'l2convergence': 1e-14,
        'adjointl2convergence': 1e-14,

        'mgcycle':'2w',
        'ncyclescoarse':250,
        'usenksolver':True,
        'useblockettes': False,
        'frozenturbulence':False,
        'blockSplitting':False,
        },
        "ref_file": 'adjoint_euler_scalar_jst_tut_wing.json', 
        'ap':ap_tutorial_wing,
        'evalFuncs':['cl', 'cd'],
        'N_PROCS': 4
    }, 

    # Tutorial wing RANS
    { "options": {
        'gridfile':os.path.join(inputDir, 'mdo_tutorial_rans_scalar_jst.cgns'),
        'restartfile':os.path.join(inputDir, 'mdo_tutorial_rans_scalar_jst.cgns'),
        # 'l2convergence': 1e-14,
        # 'l2convergencecoarse':1e-4,
        # 'ncycles':10000,

        # 'mgcycle':'sg',
        # 'equationtype':'RANS',
        # 'smoother':'dadi',

        # 'ncyclescoarse':100,

        # 'usenksolver':True,
        
        # 'nkswitchtol':1e-3,
        # 'equationtype':'RANS',
        # 'frozenturbulence':False,
        
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

        },
        "ref_file": 'adjoint_rans_tut_wing.json', 
        'ap':ap_tutorial_wing,
        'evalFuncs': ['fx', 'mz']
    }, 
    
])  
class TestAdjoint(test_objects.RegTest):
    '''
    Tests that given a flow state the residuals, function, forces/tractions, 
    and jacobian vector products are accurate.

    Test 1: MDO tutorial -- Euler -- Scalar JST
    '''
    N_PROCS = 4

    options = None
    ap = None
    ref_file = None

    def train(self):
        if self.options == None and self.ap == None and self.ref_file == None:
                # return imeediately when the train method is being called on the based class and NOT the 
                # classes created using paramertized 
                # this will happen when testing, but will hopefully be fixed down the line
                return
        super().train()

    def setUp(self, train=False):
        if self.options == None and self.ap == None and self.ref_file == None:
                # return immediately when the setup method is being called on the based class and NOT the 
                # classes created using paramertized 
                # this will happen when testing, but will hopefully be fixed down the line
                return

        super().setUp()


        options = copy.copy(adflowDefOpts)
        options.update(self.options)


        ffdFile = os.path.join(inputDir, 'mdo_tutorial_ffd.fmt')

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({
            'gridFile': options['gridfile']
        })
        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=True)


        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

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

    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap)
        
    def test_adjoint(self):
        utils.assert_adjoint_sens_allclose(self.handler, self.CFDSolver, self.ap)

