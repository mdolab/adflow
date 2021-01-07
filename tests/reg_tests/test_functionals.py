# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class

# MACH classes
from adflow import ADFLOW

# import the testing utilities
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs

from reg_aeroproblems import ap_tutorial_wing, ap_CRM, ap_tutorial_wing_laminar
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))


@parameterized_class(
    [
        # scalar JST
        {
            "name": "euler_scalar_jst_tut_wing_1core",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_scalar_jst.cgns"),
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
        # scalar JST
        {
            "name": "euler_scalar_jst_tut_wing",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_scalar_jst.cgns"),
                "l2convergence": 1e-14,
                "mgcycle": "2w",
                "ncyclescoarse": 250,
                "usenksolver": True,
                "useblockettes": False,
            },
            "ref_file": "funcs_euler_scalar_jst_tut_wing.json",
            "aero_prob": ap_tutorial_wing,
        },
        # Matrix JST
        {
            "name": "euler_matrix_jst_tut_wing",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_matrix.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_matrix.cgns"),
                "mgcycle": "2w",
                "ncyclescoarse": 250,
                "usenksolver": True,
                "nkswitchtol": 1e-2,
                "vis4": 0.1,
                "discretization": "central plus matrix dissipation",
                "coarsediscretization": "central plus matrix dissipation",
                "l2convergence": 1e-14,
                "useblockettes": False,
            },
            "ref_file": "funcs_euler_matrix_jst_tut_wing.json",
            "aero_prob": ap_tutorial_wing,
        },
        # Upwind
        {
            "name": "euler_upwind_tut_wing",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_upwind.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_euler_upwind.cgns"),
                "mgcycle": "2w",
                "ncyclescoarse": 250,
                "usenksolver": True,
                "nkswitchtol": 1e-2,
                "vis4": 0.1,
                "discretization": "upwind",
                "useblockettes": False,
                "l2convergence": 1e-14,
            },
            "ref_file": "funcs_euler_upwind_tut_wing.json",
            "aero_prob": ap_tutorial_wing,
        },
        # Tutorial wing random block order
        {
            "name": "euler_scalar_jst_rand_tut_wing",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_random_euler_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_random_euler_scalar_jst.cgns"),
                "l2convergence": 1e-14,
                "mgcycle": "2w",
                "ncyclescoarse": 250,
                "usenksolver": True,
                "useblockettes": False,
            },
            "ref_file": "funcs_euler_scalar_jst_rand_tut_wing.json",
            "aero_prob": ap_tutorial_wing,
        },
        # Tutorial wing laminar
        {
            "name": "laminar_tut_wing",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_viscous_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_viscous_scalar_jst.cgns"),
                "l2convergence": 1e-14,
                "l2convergencecoarse": 1e-4,
                "ncycles": 1000,
                "cfl": 1.5,
                "cflcoarse": 1.25,
                "mgcycle": "2w",
                "ncyclescoarse": 250,
                "usenksolver": True,
                "nkswitchtol": 1e-2,
                "equationtype": "laminar ns",
                "useblockettes": False,
            },
            "ref_file": "funcs_laminar_tut_wing.json",
            "aero_prob": ap_tutorial_wing_laminar,
        },
        # Tutorial wing RANS
        {
            "name": "rans_tut_wing",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_rans_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_rans_scalar_jst.cgns"),
                "mgcycle": "sg",
                "equationtype": "rans",
                "smoother": "dadi",
                "cfl": 1.5,
                "cflcoarse": 1.25,
                "resaveraging": "noresaveraging",
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
            },
            "ref_file": "funcs_rans_tut_wing.json",
            "aero_prob": ap_tutorial_wing,
        },
        # Tutorial wing random RANS
        {
            "name": "rans_rand_tut_wing",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_random_rans_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/mdo_tutorial_random_rans_scalar_jst.cgns"),
                "mgcycle": "sg",
                "equationtype": "rans",
                "smoother": "dadi",
                "cfl": 1.5,
                "cflcoarse": 1.25,
                "resaveraging": "noresaveraging",
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
            },
            "ref_file": "funcs_rans_rand_tut_wing.json",
            "aero_prob": ap_tutorial_wing,
        },
        # CRM WBT
        {
            "name": "euler_scalar_jst_CRM_WBT",
            "options": {
                "gridfile": os.path.join(baseDir, "../../inputFiles/CRM_wbt_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../inputFiles/CRM_wbt_scalar_jst.cgns"),
                "mgcycle": "sg",
                "cfl": 1.5,
                "cflcoarse": 1.25,
                "resaveraging": "noresaveraging",
                "ncycles": 1000,
                "monitorvariables": ["resrho", "cl", "cd", "cmy", "yplus", "totalr"],
                "usenksolver": True,
                "l2convergence": 1e-14,
                "l2convergencecoarse": 1e-4,
                "nkswitchtol": 1e-1,
                "adjointl2convergence": 1e-14,
                "liftindex": 3,
                "useblockettes": False,
            },
            "ref_file": "funcs_euler_scalar_jst_CRM_WBT.json",
            "aero_prob": ap_CRM,
        },
    ]
)
class TestFunctionals(reg_test_classes.RegTest):
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

    def test_restart_read(self):
        utils.assert_problem_size_equal(self.handler, self.CFDSolver, tol=1e-10)
        utils.assert_states_allclose(self.handler, self.CFDSolver, tol=1e-10)

    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_functions(self):
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-9)

    def test_forces_and_tractions(self):
        utils.assert_forces_allclose(self.handler, self.CFDSolver, tol=1e-10)
        utils.assert_tractions_allclose(self.handler, self.CFDSolver, tol=1e-10)

        # Reset the option
        self.CFDSolver.setOption("forcesAsTractions", True)

        # Make sure we can write the force file.
        forces_file = os.path.join(self.CFDSolver.getOption("outputdirectory"), "forces.txt")
        self.CFDSolver.writeForceFile(forces_file)

    # ------------------- Derivative routine checks ----------------------------
    def test_jac_vec_prod_fwd(self):
        utils.assert_fwd_mode_allclose(self.handler, self.CFDSolver, self.ap, tol=5e-9)

    def test_jac_vec_prod_bwd(self):
        utils.assert_bwd_mode_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_dot_products(self):
        utils.assert_dot_products_allclose(self.handler, self.CFDSolver, tol=1e-10)


if __name__ == "__main__":
    unittest.main()