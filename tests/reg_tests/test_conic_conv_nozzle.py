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
from baseclasses import BaseRegTest

# get the directories that we will use to import some packages in this repo 
baseDir = os.path.dirname(os.path.abspath(__file__))

BaseRegTest.setLocalPaths(baseDir, sys.path)

# import the ADFLOW module this test lives in 
from python.pyADflow import ADFLOW

# import the testing utilities that live a few directories up 
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts

from reg_aeroproblems import ap_conic_conv_nozzle 
from reg_test_classes import SolveRegTest, FunctionalsRegTest, AdjointRegTest


refDir, inputDir, outputDir = BaseRegTest.getLocalDirPaths(baseDir)




class TestSolveIntegrationPlane(SolveRegTest, unittest.TestCase):
    '''
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 

    Test 9: MDO tutorial -- Euler -- Solution Test
    '''
    N_PROCS = 4

    def setUp(self):

        self.ref_file = os.path.join(refDir, 'ref17.json')      

        # create the object used to compare the values to the references 
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)
        
        gridFile = os.path.join(inputDir, 'conic_conv_nozzle_mb.cgns')
        planeFile = os.path.join(inputDir, 'integration_plane_viscous.fmt')


        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'outputdirectory':outputDir,

            # Physics Parameters
            'equationType':'euler',
            'smoother':'dadi',
            'nsubiter':3,
            'CFL':4.0,
            'CFLCoarse':1.25,
            'MGCycle':'3w',
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
            # Convergence Parameters
            'L2Convergence':1e-10,
            'L2ConvergenceCoarse':1e-2,
            'adjointl2convergence':1e-6,
            'forcesAsTractions':True,
            'debugzipper':True,
            'nearwalldist':.001,
            # 'nkls':'none',
            'solutionprecision':'double',
            'adjointsubspacesize':200,
            'outerpreconits':3,
            'zipperSurfaceFamily':'output_fam',
            'flowtype':'internal',
        })

        # Setup aeroproblem
        self.ap = copy.copy(ap_conic_conv_nozzle)
        self.ap.evalFuncs.extend(['mdot_plane',
                                  'mavgptot_plane',
                                  'aavgptot_plane',
                                  'mavgttot_plane',
                                  'mavgps_plane',
                                  'aavgps_plane'])

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        self.CFDSolver.addIntegrationSurface(planeFile, 'viscous_plane')
        self.CFDSolver.finalizeUserIntegrationSurfaces()

        self.CFDSolver.addFamilyGroup('upstream',['inlet'])
        self.CFDSolver.addFamilyGroup('downstream',['outlet'])
        self.CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
        self.CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

        self.CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
        self.CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
        self.CFDSolver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

        self.CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
        self.CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
        self.CFDSolver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

        self.CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
        self.CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")
        self.CFDSolver.addFunction('aavgptot', 'viscous_plane', name="aavgptot_plane")

        self.CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
        self.CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
        self.CFDSolver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

        self.CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
        self.CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")
        self.CFDSolver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

        self.CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
        self.CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")
        self.CFDSolver.addFunction('aavgps', 'viscous_plane', name="aavgps_plane")




@unittest.skip('res is less than needed')
class TestSolveOverset(SolveRegTest, unittest.TestCase):
    '''
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 

    Test 9: MDO tutorial -- Euler -- Solution Test
    '''
    N_PROCS = 4

    def setUp(self):

        self.ref_file = os.path.join(refDir, 'ref18.json')      

        # create the object used to compare the values to the references 
        ref = utils.readJSONRef(self.ref_file)
        self.handler = BaseRegTest(ref, train=False)
        
        gridFile = os.path.join(inputDir, 'conic_conv_nozzle.cgns')


        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'outputdirectory':outputDir,

            # Physics Parameters
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
            # Convergence Parameters
            'L2Convergence':1e-10,
            'L2ConvergenceCoarse':1e-4,
            'adjointl2convergence':1e-6,
            'forcesAsTractions':True,
            'debugzipper':True,
            'nearwalldist':.001,
            # 'nkls':'none',
            'solutionprecision':'double',
            'adjointsubspacesize':200,
            'outerpreconits':3,
            'zipperSurfaceFamily':'output_fam',
            'flowtype':'internal',
        })

        # Setup aeroproblem
        self.ap = copy.copy(ap_conic_conv_nozzle)

        self.ap.setBCVar('Pressure',  79326.7, 'downstream')
        self.ap.addDV('Pressure', family='downstream')

        self.ap.setBCVar('PressureStagnation',  100000.0, 'upstream')
        self.ap.addDV('PressureStagnation', family='upstream')

        self.ap.setBCVar('TemperatureStagnation',  500.0, 'upstream')
        self.ap.addDV('TemperatureStagnation', family='upstream')

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)


        self.CFDSolver.addFamilyGroup('upstream',['inlet'])
        self.CFDSolver.addFamilyGroup('downstream',['outlet'])
        self.CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
        self.CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

        self.CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
        self.CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")

        self.CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
        self.CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")

        self.CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
        self.CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")

        self.CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
        self.CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")

        self.CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
        self.CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")

        self.CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
        self.CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")





if __name__ == '__main__':
    unittest.main()