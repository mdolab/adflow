from __future__ import print_function
import unittest
import numpy as np
import copy
from baseclasses import AeroProblem
import sys
from pprint import pprint as pp
import os
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from adflow.pyADflow_C import ADFLOW_C
# from adflow.pyADflow import ADFLOW
from mpi4py import MPI




class TestNK(unittest.TestCase):
    "class to hold the common setup method for the tests"

    def setUp(self):
        gridFile = os.path.join(baseDir, '../../inputFiles/cube.cgns')
        self.aero_options = {
            # I/O Parameters
            'gridfile': gridFile,
            'outputdirectory':'../output_files',
            # "restartfile": gridFile,
            "equationType": "RANS",
            'monitorvariables':['cpu', 'resrho', 'resturb', 'cl', 'cd'],
            "writevolumesolution": False,
            'writesurfacesolution': False,


            "MGCycle": "sg",
            "MGStartLevel": -1,
            "useANKSolver": False,
            "useNKSolver": True,

            'L2Convergence':1e-14,
            'nCycles':500,
            'adjointl2convergence': 1e-14,
            'useblockettes':False
        }


        # Aerodynamic problem description (taken from test 17)
        self.ap = AeroProblem(
            name="cube",
            V=32,  # m/s
            T=273 + 60,  # kelvin
            P=93e3,  # pa
            areaRef=1.0,  # m^2
            chordRef=1.0,  # m^2
            evalFuncs=["cd"],
            alpha=0.0,
            beta=0.00,
            xRef=0.0,
            yRef=0.0,
            zRef=0.0,
        )

        self.ap.alpha = 1e-200j

        self.CFDSolver = ADFLOW_C(options=self.aero_options)
        # self.CFDSolver.getResidual(self.ap)
        # self.CFDSolver(self.ap)

    def test_solve(self):
        self.CFDSolver(self.ap)



#endregion


if __name__ == '__main__':
    unittest.main()
