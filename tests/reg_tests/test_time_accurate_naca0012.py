# built-ins
import unittest
import numpy
import os
import copy


# MACH testing class
from adflow import ADFLOW

# import the testing utilities that live a few directories up
import reg_test_utils as utils

from reg_default_options import adflowDefOpts

from reg_aeroproblems import ap_naca0012_time_acc
from reg_test_classes import test_objects


baseDir = os.path.dirname(os.path.abspath(__file__))


class TestSolve(test_objects.RegTest):
    """
    Tests that ADflow can converge the wing from the mdo tutorial using the euler
    equation to the required accuracy as meassure by the norm of the residuals,
    and states, and the accuracy of the functions

    based on the old regression test 15
    """

    N_PROCS = 2
    ref_file = "solve_rans_time_acc_naca0012.json"

    def setUp(self):
        super().setUp()

        gridFile = os.path.join(baseDir, "../../inputFiles/naca0012_rans-L2.cgns")

        f = 10.0  # [Hz] Forcing frequency of the flow
        period = 1.0 / f  # [sec]
        nStepPerPeriod = 8
        nPeriods = 1
        nfineSteps = nStepPerPeriod * nPeriods
        dt = period / nStepPerPeriod  # [s] The actual timestep

        options = copy.copy(adflowDefOpts)
        options.update(
            {
                "gridfile": gridFile,
                "outputdirectory": os.path.join(baseDir, "../output_files"),
                "writevolumesolution": False,
                "vis4": 0.025,
                "vis2": 0.5,
                "restrictionrelaxation": 0.5,
                "smoother": "dadi",
                "equationtype": "RANS",
                "equationmode": "unsteady",
                "timeIntegrationscheme": "bdf",
                "ntimestepsfine": nfineSteps,
                "deltat": dt,
                "nsubiterturb": 10,
                "nsubiter": 5,
                "useale": False,
                "usegridmotion": True,
                "cfl": 2.5,
                "cflcoarse": 1.2,
                "ncycles": 2000,
                "mgcycle": "3w",
                "mgstartlevel": 1,
                "monitorvariables": ["cpu", "resrho", "cl", "cd", "cmz"],
                "usenksolver": False,
                "l2convergence": 1e-6,
                "l2convergencecoarse": 1e-4,
                "qmode": True,
                "alphafollowing": False,
                "blockSplitting": True,
                "useblockettes": False,
            }
        )

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca0012_time_acc)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        self.CFDSolver.addSlices("z", [0.5])

    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, rtol=1e-8, atol=1e-8)


if __name__ == "__main__":
    unittest.main()