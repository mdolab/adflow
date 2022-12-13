# built-ins
import unittest
import numpy as np
import os
import copy

# MACH classes
from adflow import ADFLOW
from baseclasses.utils import Error

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
                "mgcycle": "sg",
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

    def test_basic_solve(self):

        self.CFDSolver.solveCL(self.ap, 0.475, alpha0=1.20, delta=0.025, tol=1e-4, autoReset=False)
        self.assert_solution_failure()
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["cl", "cavitation"])
        self.handler.root_add_val("CL-CL*", funcs["mdo_tutorial_cl"] - 0.475, rtol=1e-4, atol=1e-4)
        self.handler.root_add_val("cavitation", funcs["mdo_tutorial_cavitation"], rtol=1e-4, atol=1e-4)

    def test_optional_features(self):

        CLStar = 0.475
        tol = 1e-6
        clSolveRes = self.CFDSolver.solveCL(
            self.ap,
            CLStar,
            alpha0=1.20,
            CLalphaGuess=0.12,  # the actual value is sth like 0.139, this is offset on purpose
            tol=tol,
            relaxCLa=0.9,
            L2ConvRel=[1e-4, 1e-3],
        )
        self.assert_solution_failure()
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["cl", "cavitation"])

        np.testing.assert_allclose(CLStar, funcs["mdo_tutorial_cl"], rtol=tol)
        np.testing.assert_equal(True, clSolveRes["converged"])

    def test_clsolve_l2_failure(self):

        CLStar = 0.475
        tol = 1e-6
        clSolveRes = self.CFDSolver.solveCL(
            self.ap,
            CLStar,
            maxIter=8,
            alpha0=1.20,
            CLalphaGuess=0.12,
            tol=tol,
            relaxCLa=0.9,
            L2ConvRel=[1e-4, 1e-2],  # the second value is high to force a failure in the CFD convergence
        )
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["cl", "cavitation"])

        # this test is designed to fail with the CFD L2 convergence, but the CL result
        # will meet the error tolerance. we still want the cl solve to return a fail flag
        # since we failed to fully converge the CFD!
        np.testing.assert_allclose(CLStar, funcs["mdo_tutorial_cl"], rtol=tol)
        np.testing.assert_equal(False, clSolveRes["converged"])

    def test_clsolve_cl_failure(self):

        CLStar = 0.475
        tol = 1e-6
        clSolveRes = self.CFDSolver.solveCL(
            self.ap,
            CLStar,
            maxIter=3,
            alpha0=1.20,
            CLalphaGuess=0.12,
            tol=tol,
            relaxCLa=0.9,
            L2ConvRel=1e-4,
        )
        self.assert_solution_failure()
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["cl", "cavitation"])

        # this test should run out of iterations, so clsolve should return fail
        np.testing.assert_equal(False, clSolveRes["converged"])

    def test_clsolve_cfd_stall(self):
        self.CFDSolver.setOption("nCycles", 1)

        CLStar = 0.475
        tol = 1e-6

        # because errOnStall is set to True, the clsolve should raise an error
        with np.testing.assert_raises(Error):
            self.CFDSolver.solveCL(
                self.ap,
                CLStar,
                maxIter=10,
                alpha0=1.20,
                CLalphaGuess=0.12,
                tol=tol,
                relaxCLa=0.9,
                L2ConvRel=1e-4,
                errOnStall=True,
            )


if __name__ == "__main__":
    unittest.main()
