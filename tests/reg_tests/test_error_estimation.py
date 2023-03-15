# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
from numpy.testing import assert_allclose

# MACH classes
from adflow import ADFLOW

# import the testing utilities
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_tutorial_wing
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))
evalFuncs = ["fx", "mz", "cl", "cd", "cmz", "lift", "drag"]


@parameterized_class(
    [
        # fully converged
        {
            "name": "euler_wing_full",
            "options": {
                "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
                "L2Convergence": 1e-14,
                "adjointL2Convergence": 1e-14,
                "mgcycle": "sg",
                "useNKSolver": True,
                "NKSwitchtol": 1e-5,
            },
            "ref_file": "error_euler_wing_full.json",
            "aero_prob": ap_tutorial_wing,
            "evalFuncs": evalFuncs,
        },
        # partially converged
        {
            "name": "euler_wing_partial",
            "options": {
                "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
                "L2Convergence": 1e-6,
                "adjointL2Convergence": 1e-14,
                "mgcycle": "sg",
                "ncycles": 5000,
                "useNKSolver": True,
                "NKSwitchtol": 1e-5,
            },
            "ref_file": "error_euler_wing_partial.json",
            "aero_prob": ap_tutorial_wing,
            "evalFuncs": evalFuncs,
        },
    ]
)
class TestError(reg_test_classes.RegTest):
    """
    Tests that ADflow can converge the given test problems to the given tolerance.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.ap = copy.deepcopy(self.aero_prob)
        self.ap.evalFuncs = self.evalFuncs

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

    def test_error(self):
        # do the solve
        self.CFDSolver(self.ap)
        self.assert_solution_failure()

        # solve the adjoint
        funcsSens = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSens)
        self.assert_adjoint_failure()

        # compute the error
        funcsError = {}
        self.CFDSolver.solveErrorEstimate(self.ap, funcsError)

        # if fully converged, check that the error is very small
        if "full" in self.name:
            for value in funcsError.values():
                assert_allclose(value, 0, atol=1e-7)

        self.handler.root_add_dict("errors", funcsError, tol=1e-5)


if __name__ == "__main__":
    unittest.main()
