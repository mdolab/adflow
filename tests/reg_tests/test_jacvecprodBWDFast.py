# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts, defaultAeroDVs

from reg_aeroproblems import ap_tutorial_wing
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
    # scalar JST
    {
        "name": "euler_scalar_jst_tut_wing_1core",
        "options": {
            "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "l2convergence": 1e-14,
            "mgcycle": "2w",
            "ncyclescoarse": 250,
            "usenksolver": True,
            "useblockettes": False,
        },
        "ref_file": "funcs_euler_scalar_jst_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "N_PROCS": 1,
    },
]


@parameterized_class(test_params)
class TestJacVecBWDFast(reg_test_classes.RegTest):
    """
    Tests that given a flow state the residuals, function, forces/tractions,
    and jacobian vector products are accurate.

    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when testing, but will hopefully be fixed down the line
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        # Create the solver
        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    # ------------------- Derivative routine checks ----------------------------
    def test_BWD(self):
        #

        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBar = self.CFDSolver.computeJacobianVectorProductBwd(
            resBar=dwBar,
            wDeriv=True,
        )

        wBarfast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBar, wBarfast, atol=1e-16, err_msg="w wrt res")

    def test_repeated_calls(self):

        dwBar = self.CFDSolver.getStatePerturbation(314)

        wBarfast1 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        wBarfast2 = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=dwBar)

        np.testing.assert_allclose(wBarfast1, wBarfast2, atol=1e-16, err_msg="w wrt res double call")


if __name__ == "__main__":
    unittest.main()
