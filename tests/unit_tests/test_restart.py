import unittest
import os
from adflow import ADFLOW
import numpy as np
import sys

baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir, "../reg_tests"))

from reg_aeroproblems import ap_tutorial_wing  # noqa E402


class BasicTests(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        gridFile = "input_files/mdo_tutorial_euler_scalar_jst.cgns"
        self.options = {
            "gridfile": os.path.join(baseDir, "../../", gridFile),
            "restartFile": os.path.join(baseDir, "../../", gridFile),
            "MGCycle": "sg",
            "equationType": "Euler",
        }

    def test_import(self):
        CFDSolver = ADFLOW(options=self.options, debug=False)
        res = CFDSolver.getResidual(ap_tutorial_wing)
        res_norm = np.linalg.norm(res)
        np.testing.assert_allclose(res_norm, 0.0, atol=1e-11, err_msg="residual")

    def test_import_block_splitting(self):
        self.options["partitionLikeNProc"] = 50
        CFDSolver = ADFLOW(options=self.options, debug=False)
        res = CFDSolver.getResidual(ap_tutorial_wing)
        res_norm = np.linalg.norm(res)
        np.testing.assert_allclose(res_norm, 0.0, atol=1e-11, err_msg="residual")


if __name__ == "__main__":
    unittest.main()
