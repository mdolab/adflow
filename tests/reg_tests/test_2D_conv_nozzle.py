# built-ins
import unittest
import os
import copy

# MACH classes
from adflow import ADFLOW

import reg_test_utils as utils
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_2D_conv_nozzle
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestSolve(reg_test_classes.RegTest):
    """
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regresson test 16
    """

    N_PROCS = 2
    ref_file = "solve_2D_conv_nozzle.json"
    options = {
        "gridfile": os.path.join(baseDir, "../../inputFiles/euler_conv_nozzle.cgns"),
        "solutionprecision": "double",
        "gridprecision": "double",
        "liftIndex": 2,
        "CFL": 3.0,
        "CFLCoarse": 1.5,
        "MGCycle": "2w",
        "MGStartLevel": 2,
        "nCyclesCoarse": 500,
        "nCycles": 2500,
        "monitorvariables": ["resrho", "cl", "cd", "yplus"],
        "nsubiterturb": 3,
        "useNKSolver": True,
        "NKSubSpaceSize": 60,
        "L2Convergence": 1e-14,
        "L2ConvergenceCoarse": 1e-2,
        "NKSwitchTol": 1e-2,
        "nkadpc": False,
        "vis4": 0.006,
        "vis2": 0.0,
        "blocksplitting": True,
        "solutionPrecision": "double",
        "flowtype": "internal",
    }
    ap = ap_2D_conv_nozzle

    def setUp(self):
        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

        self.CFDSolver.addFamilyGroup("upstream", ["INFLOW"])
        self.CFDSolver.addFamilyGroup("downstream", ["OUTFLOW"])
        self.CFDSolver.addFamilyGroup("all_flow", ["INFLOW", "OUTFLOW"])
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

        self.CFDSolver.addFunction("mavgmn", "downstream", name="mavgmn_down")
        self.CFDSolver.addFunction("mavgmn", "upstream", name="mavgmn_up")

        self.CFDSolver.addFunction("drag", "all_flow", name="thrust")  # this naming makes it seem like wishful thinking

        self.CFDSolver.addFunction("dragpressure", "all_flow", name="thrust_pressure")
        self.CFDSolver.addFunction("dragviscous", "all_flow", name="thrust_viscous")
        self.CFDSolver.addFunction("dragmomentum", "all_flow", name="thrust_momentum")

    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-9)
        utils.assert_states_allclose(self.handler, self.CFDSolver, tol=1e-10)
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)


if __name__ == "__main__":
    unittest.main()