# built-ins
from importlib import util
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np

# MACH classes
from adflow import ADFLOW, ADFLOW_C

from reg_default_options import adflowDefOpts, defaultAeroDVs
import reg_test_utils as utils
from baseclasses.testing import getTol

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
class TestJacVecFwd(reg_test_classes.RegTest):
    """
    Tests that given a flow state the FWD jacobian vector products are agree the privous values recorded in the ref file.
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
    def test_wDot(self):
        utils.assert_fwd_mode_wdot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=5e-9)

    def test_xVDot(self):
        utils.assert_fwd_mode_xVDot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=1e-10)

    def test_xDvDot(self):
        utils.assert_fwd_mode_xDvDot_allclose(self.handler, self.CFDSolver, self.ap, seed=1.0, tol=1e-10)


@parameterized_class(test_params)
class TestJacVecFwdFD(reg_test_classes.RegTest):
    """
    Tests that given a flow state the FWD jacobian vector products are agree with FD.
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
    def test_wDot(self):
        # perturb each input and check that the outputs match the FD to with in reason
        wDot = self.CFDSolver.getStatePerturbation(321)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True
        )
        resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=1e-8
        )

        np.testing.assert_allclose(resDot_FD, resDot, rtol=8e-4, err_msg="residual")

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func], rtol=1e-5, err_msg=func)

        np.testing.assert_allclose(fDot_FD, fDot, rtol=5e-4, err_msg="forces")

    def test_xVDot(self):
        # perturb each input and check that the outputs match the FD to with in reason
        xVDot = self.CFDSolver.getSpatialPerturbation(314)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True
        )

        resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=1e-8
        )

        idx_max = np.argmax((resDot_FD - resDot) / resDot)
        print(resDot[idx_max], resDot_FD[idx_max])

        np.testing.assert_allclose(resDot_FD, resDot, atol=5e-4, err_msg="residual")

        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func], rtol=5e-6, err_msg=func)

        np.testing.assert_allclose(fDot_FD, fDot, rtol=5e-4, err_msg="forces")

    def test_xDvDot(self):
        # perturb each input and check that the outputs match the FD to with in reason
        step_size = {
            "alpha": 1e-4,
            "beta": 1e-5,
            "mach": 1e-5,
            "P": 1e-1,
            "T": 1e-4,
            "xRef": 1e-5,
            "yRef": 1e-5,
            "zRef": 1e-5,
        }

        for aeroDV in self.ap.DVs.values():
            key = aeroDV.key
            xDvDot = {key: 1.0}

            resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True
            )

            resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="FD", h=step_size[key]
            )

            # the tolerances here are loose becuase different ouputs have different optimal steps
            np.testing.assert_allclose(resDot_FD, resDot, atol=5e-5, err_msg=f"residual wrt {key}")

            for func in funcsDot:
                if np.abs(funcsDot[func]) <= 1e-16:
                    np.testing.assert_allclose(
                        funcsDot_FD[func], funcsDot[func], atol=5e-5, err_msg=f"{func} wrt {key}"
                    )
                else:
                    np.testing.assert_allclose(
                        funcsDot_FD[func], funcsDot[func], rtol=1e-3, err_msg=f"{func} wrt {key}"
                    )

            np.testing.assert_allclose(fDot_FD, fDot, atol=5e-7, err_msg=f"forces wrt {key}")


@parameterized_class(test_params)
class TestJacVecFwdCS(reg_test_classes.CmplxRegTest):
    """
    Tests jacobian vector products against CS.

    """

    N_PROCS = 2

    h = 1e-40

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

        # Setup aeroproblem
        self.ap = copy.deepcopy(self.aero_prob)
        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver = ADFLOW_C(options=options, debug=True)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    # ------------------- Derivative routine checks ----------------------------
    def cmplx_test_wDot(self):

        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # perturb each input and check that the outputs match the FD to with in reason
        wDot = self.CFDSolver.getStatePerturbation(314)

        resDot_CS, funcsDot_CS, fDot_CS = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        rtol, atol = getTol(tol=1e-9)

        self.handler.root_print("||dR/dw * wDot||")
        self.handler.par_add_norm("||dR/dw * wDot||", resDot_CS, rtol=rtol, atol=atol)

        self.handler.root_print("dFuncs/dw * wDot")
        self.handler.root_add_dict("dFuncs/dw * wDot", funcsDot_CS, rtol=rtol, atol=atol)

        self.handler.root_print("||dF/dw * wDot||")
        self.handler.par_add_norm("||dF/dw * wDot||", fDot_CS, rtol=rtol, atol=atol)

    def cmplx_test_xVDot(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # perturb each input and check that the outputs match the FD to with in reason
        xVDot = self.CFDSolver.getSpatialPerturbation(314)

        rtol, atol = getTol(tol=1e-10)

        resDot_cs, funcsDot_cs, fDot_cs = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        self.handler.root_print("||dR/dXv * xVDot||")
        self.handler.par_add_norm("||dR/dXv * xVDot||", resDot_cs, rtol=rtol, atol=atol)

        # These can be finiky sometimes so a bigger tolerance.
        self.handler.root_print("dFuncs/dXv * xVDot")
        self.handler.root_add_dict("dFuncs/dXv * xVDot", funcsDot_cs, rtol=rtol * 10, atol=atol * 10)

        self.handler.root_print("||dF/dXv * xVDot||")
        self.handler.par_add_norm("||dF/dXv * xVDot||", fDot_cs, rtol=rtol, atol=atol)

    def cmplx_test_xDvDot(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return
        rtol, atol = getTol(tol=1e-10)

        # perturb each input and check that the outputs match the FD to with in reason
        for aeroDV in self.ap.DVs.values():
            key = aeroDV.key
            self.handler.root_print("  -> %s" % key)
            xDvDot = {key: 1.0}

            resDot_cs, funcsDot_cs, fDot_cs = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
            )

            self.handler.root_print("||dR/d%s||" % key)
            self.handler.par_add_norm("||dR/d%s||" % key, resDot_cs, rtol=rtol, atol=atol)

            self.handler.root_print("dFuncs/d%s" % key)
            self.handler.root_add_dict("dFuncs/d%s" % key, funcsDot_cs, rtol=rtol, atol=atol)

            self.handler.root_print("||dF/d%s||" % key)
            self.handler.par_add_norm("||dF/d%s||" % key, fDot_cs, rtol=rtol, atol=atol)


if __name__ == "__main__":
    unittest.main()
