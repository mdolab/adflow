# built-ins
import os
import copy
from parameterized import parameterized_class

# MACH classes
from adflow import ADFLOW

# import the testing utilities
import reg_test_utils as utils
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_naca0012_time_acc
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))


@parameterized_class(
    [
        # no restart
        {
            "name": "correction_disabled",
            "options": {
                "infChangeCorrection": False,
            },
            "ref_file": "flow_restart_correction_disabled.json",
            "aero_prob": ap_naca0012_time_acc,
        },
        # offset
        {
            "name": "correction_offset",
            "options": {
                "infChangeCorrection": True,
                "infChangeCorrectionType": "offset",
            },
            "ref_file": "flow_restart_correction_offset.json",
            "aero_prob": ap_naca0012_time_acc,
        },
        # rotate
        {
            "name": "correction_rotate",
            "options": {
                "infChangeCorrection": True,
                "infChangeCorrectionType": "rotate",
            },
            "ref_file": "flow_restart_correction_rotate.json",
            "aero_prob": ap_naca0012_time_acc,
        },
    ]
)
class TestRestart(reg_test_classes.RegTest):
    """
    Test different flow restart algorithms
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        super().setUp()
        gridFile = os.path.join(baseDir, "../../input_files/naca0012_rans-L2.cgns")
        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(
            {
                "gridfile": gridFile,
                "mgcycle": "sg",
                "equationtype": "RANS",
                "smoother": "DADI",
                "cfl": 1.5,
                "cflcoarse": 1.25,
                "solutionprecision": "double",
                "gridprecision": "double",
                "mgcycle": "sg",
                "anksecondordswitchtol": 1e-2,
                "ankcoupledswitchtol": 1e-4,
                "l2convergence": 1e-13,
                "writeVolumeSolution": False,
                "writeSurfaceSolution": False,
            }
        )
        options.update(self.options)

        self.ap = copy.deepcopy(self.aero_prob)
        self.ap.addDV("alpha", value=0.0)
        self.ap.addDV("mach", value=0.6)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

    def test_restart_correction(self):
        # start at aoa=0
        self.ap.setDesignVars({f"alpha_{self.ap.name}": 0.0})
        # do the solve
        self.CFDSolver(self.ap)

        # change AoA and Mach number
        self.ap.setDesignVars({f"alpha_{self.ap.name}": 4.0})
        self.ap.setDesignVars({f"mach_{self.ap.name}": 0.8})

        # check the changes in residual. the getRes method calls setAeroProblem,
        # which runs the correction update if requested.
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)
