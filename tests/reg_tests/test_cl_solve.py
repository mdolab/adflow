# built-ins
import unittest
import os
import copy

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_tutorial_wing
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))


class TestSolve(reg_test_classes.RegTest):
    """
    Tests that ADflow can change the dvs to achieve the target cl to the given tolerance.

    Compares the norm of the residuals, norm of states, and the functionals.

    based on the old regression test 10

    """

    N_PROCS = 2
    ref_file = "solve_cl.json"

    def setUp(self):
        super().setUp()
        gridFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns")
        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(
            {
                "gridfile": gridFile,
                "solutionprecision": "double",
                "gridprecision": "double",
                "mgcycle": "2w",
                "ncyclescoarse": 250,
                "ncycles": 500,
                "usenksolver": True,
                "nkswitchtol": 1e-2,
                "l2convergence": 1e-14,
                "l2convergencecoarse": 1e-2,
                "computecavitation": True,
            }
        )

        # Setup aeroproblem
        self.ap = copy.copy(ap_tutorial_wing)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

    def test_solve(self):
        self.CFDSolver.solveCL(self.ap, 0.475, alpha0=1.20, delta=0.025, tol=1e-4, autoReset=False)
        self.assert_solution_failure()
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["cl", "cavitation"])
        self.handler.root_add_val("CL-CL*", funcs["mdo_tutorial_cl"] - 0.475, rtol=1e-4, atol=1e-4)
        self.handler.root_add_val("cavitation", funcs["mdo_tutorial_cavitation"], rtol=1e-4, atol=1e-4)


if __name__ == "__main__":
    unittest.main()