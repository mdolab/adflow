# built-ins
import unittest
import numpy as np
import os
import sys
import copy

# MACH classes
from adflow import ADFLOW

# import the testing utilities
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts

from reg_aeroproblems import ap_conic_conv_nozzle
from reg_test_classes import test_objects


baseDir = os.path.dirname(os.path.abspath(__file__))


class TestSolveIntegrationPlane(test_objects.RegTest):
    """
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test

    Test 9: MDO tutorial -- Euler -- Solution Test
    """

    N_PROCS = 2

    options = {
        "gridfile": os.path.join(baseDir, "../../inputFiles/conic_conv_nozzle_mb.cgns"),
        "outputdirectory": os.path.join(baseDir, "../output_files"),
        # Physics Parameters
        "equationType": "euler",
        "smoother": "dadi",
        "nsubiter": 3,
        "CFL": 4.0,
        "CFLCoarse": 1.25,
        "MGCycle": "3w",
        "MGStartLevel": -1,
        "nCyclesCoarse": 250,
        "nCycles": 1000,
        "nkcfl0": 1e10,
        "monitorvariables": ["cpu", "resrho", "cl", "cd"],
        "volumevariables": ["blank"],
        "surfacevariables": ["mach", "cp", "vx", "vy", "vz", "blank"],
        "useNKSolver": True,
        "nkswitchtol": 0.01,
        "nkadpc": True,
        "nkjacobianlag": 5,
        "nkouterpreconits": 3,
        "nkinnerpreconits": 2,
        # Convergence Parameters
        "L2Convergence": 1e-10,
        "L2ConvergenceCoarse": 1e-2,
        "adjointl2convergence": 1e-6,
        "forcesAsTractions": True,
        "debugzipper": True,
        "nearwalldist": 0.001,
        # 'nkls':'none',
        "solutionprecision": "double",
        "adjointsubspacesize": 200,
        "outerpreconits": 3,
        "zipperSurfaceFamily": "output_fam",
        "flowtype": "internal",
    }
    ref_file = "solve_conic_mb.json"

    def setUp(self):

        super().setUp()

        self.ap = copy.copy(ap_conic_conv_nozzle)
        options = copy.copy(adflowDefOpts)
        options.update(self.options)

        # Setup aeroproblem

        planeFile = os.path.join(baseDir, "../../inputFiles/integration_plane_viscous.fmt")

        options = copy.copy(adflowDefOpts)
        options.update(self.options)

        # Setup aeroproblem
        self.ap.evalFuncs.extend(
            ["mdot_plane", "mavgptot_plane", "aavgptot_plane", "mavgttot_plane", "mavgps_plane", "aavgps_plane"]
        )

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        self.CFDSolver.addIntegrationSurface(planeFile, "viscous_plane")
        self.CFDSolver.finalizeUserIntegrationSurfaces()

        self.CFDSolver.addFamilyGroup("upstream", ["inlet"])
        self.CFDSolver.addFamilyGroup("downstream", ["outlet"])
        self.CFDSolver.addFamilyGroup("all_flow", ["inlet", "outlet"])
        self.CFDSolver.addFamilyGroup("output_fam", ["all_flow", "allWalls"])

        self.CFDSolver.addFunction("mdot", "upstream", name="mdot_up")
        self.CFDSolver.addFunction("mdot", "downstream", name="mdot_down")
        self.CFDSolver.addFunction("mdot", "viscous_plane", name="mdot_plane")

        self.CFDSolver.addFunction("mavgptot", "downstream", name="mavgptot_down")
        self.CFDSolver.addFunction("mavgptot", "upstream", name="mavgptot_up")
        self.CFDSolver.addFunction("mavgptot", "viscous_plane", name="mavgptot_plane")

        self.CFDSolver.addFunction("aavgptot", "downstream", name="aavgptot_down")
        self.CFDSolver.addFunction("aavgptot", "upstream", name="aavgptot_up")
        self.CFDSolver.addFunction("aavgptot", "viscous_plane", name="aavgptot_plane")

        self.CFDSolver.addFunction("mavgttot", "downstream", name="mavgttot_down")
        self.CFDSolver.addFunction("mavgttot", "upstream", name="mavgttot_up")
        self.CFDSolver.addFunction("mavgttot", "viscous_plane", name="mavgttot_plane")

        self.CFDSolver.addFunction("mavgps", "downstream", name="mavgps_down")
        self.CFDSolver.addFunction("mavgps", "upstream", name="mavgps_up")
        self.CFDSolver.addFunction("mavgps", "viscous_plane", name="mavgps_plane")

        self.CFDSolver.addFunction("aavgps", "downstream", name="aavgps_down")
        self.CFDSolver.addFunction("aavgps", "upstream", name="aavgps_up")
        self.CFDSolver.addFunction("aavgps", "viscous_plane", name="aavgps_plane")

    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap)
        utils.assert_states_allclose(self.handler, self.CFDSolver)
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap)


class TestSolveOverset(test_objects.RegTest):
    """
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 17 and 18
    """

    N_PROCS = 2

    options = {
        "gridfile": os.path.join(baseDir, "../../inputFiles/conic_conv_nozzle.cgns"),
        "outputdirectory": os.path.join(baseDir, "../output_files"),
        # Physics Parameters
        "equationType": "euler",
        "smoother": "dadi",
        "nsubiter": 3,
        "CFL": 4.0,
        "CFLCoarse": 1.25,
        "MGCycle": "sg",
        "MGStartLevel": -1,
        "nCyclesCoarse": 250,
        "nCycles": 1000,
        "nkcfl0": 1e10,
        "monitorvariables": ["cpu", "resrho", "cl", "cd"],
        "volumevariables": ["blank"],
        "surfacevariables": ["mach", "cp", "vx", "vy", "vz", "blank"],
        "useNKSolver": True,
        "nkswitchtol": 0.01,
        "nkadpc": True,
        "nkjacobianlag": 5,
        "nkouterpreconits": 3,
        "nkinnerpreconits": 2,
        # Convergence Parameters
        "L2Convergence": 1e-10,
        "L2ConvergenceCoarse": 1e-4,
        "adjointl2convergence": 1e-6,
        "forcesAsTractions": True,
        "debugzipper": True,
        "nearwalldist": 0.001,
        # 'nkls':'none',
        "solutionprecision": "double",
        "adjointsubspacesize": 200,
        "outerpreconits": 3,
        "zipperSurfaceFamily": "output_fam",
        "flowtype": "internal",
        "blocksplitting": True,
    }
    ap = copy.copy(ap_conic_conv_nozzle)
    ref_file = "solve_conic_overset.json"

    def setUp(self):
        super().setUp()

        options = copy.copy(adflowDefOpts)
        options.update(self.options)

        # Setup aeroproblem

        self.ap.setBCVar("Pressure", 79326.7, "downstream")
        self.ap.addDV("Pressure", family="downstream")

        self.ap.setBCVar("PressureStagnation", 100000.0, "upstream")
        self.ap.addDV("PressureStagnation", family="upstream")

        self.ap.setBCVar("TemperatureStagnation", 500.0, "upstream")
        self.ap.addDV("TemperatureStagnation", family="upstream")

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

        self.CFDSolver.addFamilyGroup("upstream", ["inlet"])
        self.CFDSolver.addFamilyGroup("downstream", ["outlet"])
        self.CFDSolver.addFamilyGroup("all_flow", ["inlet", "outlet"])
        self.CFDSolver.addFamilyGroup("output_fam", ["all_flow", "allWalls"])

        self.CFDSolver.addFunction("mdot", "upstream", name="mdot_up")
        self.CFDSolver.addFunction("mdot", "downstream", name="mdot_down")

        self.CFDSolver.addFunction("mavgptot", "downstream", name="mavgptot_down")
        self.CFDSolver.addFunction("mavgptot", "upstream", name="mavgptot_up")

        self.CFDSolver.addFunction("aavgptot", "downstream", name="aavgptot_down")
        self.CFDSolver.addFunction("aavgptot", "upstream", name="aavgptot_up")

        self.CFDSolver.addFunction("mavgttot", "downstream", name="mavgttot_down")
        self.CFDSolver.addFunction("mavgttot", "upstream", name="mavgttot_up")

        self.CFDSolver.addFunction("mavgps", "downstream", name="mavgps_down")
        self.CFDSolver.addFunction("mavgps", "upstream", name="mavgps_up")

        self.CFDSolver.addFunction("aavgps", "downstream", name="aavgps_down")
        self.CFDSolver.addFunction("aavgps", "upstream", name="aavgps_up")

    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap)
        utils.assert_states_allclose(self.handler, self.CFDSolver)
        # Check the residual
        res = self.CFDSolver.getResidual(self.ap)
        totalR0 = self.CFDSolver.getFreeStreamResidual(self.ap)
        res /= totalR0

        reducedSum = self.CFDSolver.comm.reduce(np.sum(res ** 2))
        if self.CFDSolver.comm.rank == 0:
            self.assertLessEqual(np.sqrt(reducedSum), self.options["L2Convergence"])


if __name__ == "__main__":
    unittest.main()