# built-ins
import unittest
import os
import copy
import numpy as np

# MACH classes
from adflow import ADFLOW

from reg_default_options import adflowDefOpts, defaultAeroDVs

from reg_aeroproblems import ap_tutorial_wing

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestCenterOfForce(unittest.TestCase):
    """
    Test center of force outputs' consistency with moment outputs,
    which are all computed separately in the fortran layer.

    """

    N_PROCS = 2

    def setUp(self):

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(
            {
                "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
                "mgcycle": "sg",
                "equationtype": "RANS",
                "smoother": "DADI",
                "cfl": 1.5,
                "cflcoarse": 1.25,
                "resaveraging": "never",
                "nsubiter": 3,
                "nsubiterturb": 3,
                "ncyclescoarse": 100,
                "ncycles": 1000,
                "monitorvariables": ["cpu", "resrho", "resturb", "cl", "cd", "cmz", "yplus", "totalr"],
                "usenksolver": True,
                "l2convergence": 1e-14,
                "l2convergencecoarse": 1e-4,
                "nkswitchtol": 1e-3,
                "adjointl2convergence": 1e-14,
                "frozenturbulence": False,
            }
        )

        # Create the solver
        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)
        self.ap = copy.deepcopy(ap_tutorial_wing)
        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def test_cof(self):
        # We should be able to match the moment calculations
        # using the forces and center of force values.
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        # forces from CFD
        fx = funcs[f"{self.ap.name}_fx"]
        fy = funcs[f"{self.ap.name}_fy"]
        fz = funcs[f"{self.ap.name}_fz"]

        # moments from CFD
        mx = funcs[f"{self.ap.name}_mx"]
        my = funcs[f"{self.ap.name}_my"]
        mz = funcs[f"{self.ap.name}_mz"]

        # center of x-force
        cofxx = funcs[f"{self.ap.name}_cofxx"]
        cofxy = funcs[f"{self.ap.name}_cofxy"]
        cofxz = funcs[f"{self.ap.name}_cofxz"]

        # center of y-force
        cofyx = funcs[f"{self.ap.name}_cofyx"]
        cofyy = funcs[f"{self.ap.name}_cofyy"]
        cofyz = funcs[f"{self.ap.name}_cofyz"]

        # center of z-force
        cofzx = funcs[f"{self.ap.name}_cofzx"]
        cofzy = funcs[f"{self.ap.name}_cofzy"]
        cofzz = funcs[f"{self.ap.name}_cofzz"]

        # reference point for the AP
        xref = self.ap.xRef
        yref = self.ap.yRef
        zref = self.ap.zRef

        mymx = fz * (cofzy - yref) - fy * (cofyz - zref)
        mymy = fx * (cofxz - zref) - fz * (cofzx - xref)
        mymz = fy * (cofyx - xref) - fx * (cofxy - yref)

        np.testing.assert_allclose(mx, mymx, rtol=1e-9)
        np.testing.assert_allclose(my, mymy, rtol=1e-9)
        np.testing.assert_allclose(mz, mymz, rtol=1e-9)


if __name__ == "__main__":
    unittest.main()
