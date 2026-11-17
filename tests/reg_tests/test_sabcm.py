# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np
from mpi4py import MPI

# MACH classes
from adflow import ADFLOW, ADFLOW_C

from reg_default_options import adflowDefOpts, defaultAeroDVs
import reg_test_utils as utils
from baseclasses.testing import getTol

from reg_aeroproblems import ap_tutorial_wing
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))

# The SA-BCM residuals and derivatives are checked at two states built from the
# converged fully turbulent SA solution of the tutorial wing. Neither is a converged
# SA-BCM solution, but the residuals and the Jacobian-vector products are valid at
# any state.
# * "turbulent": the SA solution itself. In the boundary layer Term2 saturates the
#   intermittency, so this state exercises Term2 and the intermittency blend.
# * "transitional": the SA solution with nuTilde scaled by 0.01 and Tu = 2 %. The
#   boundary layer is then laminar or transitional and Term1 (Re_theta against the
#   critical Re_theta) is active as well.
sabcmOptions = {
    "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
    "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
    "equationType": "RANS",
    "useSABCM": True,
    "turbIntensityInf": 0.005,
}

test_params = [
    {
        "name": "rans_sabcm_smooth_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": True},
        "ref_file": "funcs_rans_sabcm_smooth_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 1.0,
    },
    {
        "name": "rans_sabcm_original_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": False},
        "ref_file": "funcs_rans_sabcm_original_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 1.0,
    },
    {
        "name": "rans_sabcm_smooth_transitional_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": True, "turbIntensityInf": 0.02},
        "ref_file": "funcs_rans_sabcm_smooth_transitional_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 0.01,
    },
    {
        "name": "rans_sabcm_original_transitional_tut_wing",
        "options": {**sabcmOptions, "SABCMSmooth": False, "turbIntensityInf": 0.02},
        "ref_file": "funcs_rans_sabcm_original_transitional_tut_wing.json",
        "aero_prob": copy.deepcopy(ap_tutorial_wing),
        "nuTildeScale": 0.01,
    },
]


def setTestState(CFDSolver, ap, nuTildeScale):
    """Propagate the restart state and scale nuTilde to build the test state."""
    CFDSolver.getResidual(ap)
    if nuTildeScale != 1.0:
        nw = CFDSolver.adflow.flowvarrefstate.nw
        w = CFDSolver.getStates()
        w.reshape(-1, nw)[:, -1] *= nuTildeScale
        CFDSolver.setStates(w)
        CFDSolver.getResidual(ap)


@parameterized_class(test_params)
class TestSABCM(reg_test_classes.RegTest):
    """
    Tests that given a flow state the residuals and the FWD jacobian vector products
    of the SA-BCM model agree with the values recorded in the ref file.
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        setTestState(self.CFDSolver, self.ap, self.nuTildeScale)

    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

        # The model only changes the turbulence equation, so check its residual on its own as well
        rtol, atol = getTol(tol=1e-10)
        res = self.CFDSolver.getResidual(self.ap) / self.CFDSolver.getFreeStreamResidual(self.ap)
        nw = self.CFDSolver.adflow.flowvarrefstate.nw
        self.handler.root_print("Norm of turbulence residual")
        self.handler.par_add_norm("Norm of turbulence residual", res.reshape(-1, nw)[:, -1], rtol=rtol, atol=atol)

    def test_wDot(self):
        utils.assert_fwd_mode_wdot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=5e-9)

    def test_xVDot(self):
        utils.assert_fwd_mode_xVDot_allclose(self.handler, self.CFDSolver, self.ap, seed=314, tol=1e-10)


@parameterized_class(test_params)
class TestSABCMJacVecFwdCS(reg_test_classes.CmplxRegTest):
    """
    Tests the FWD jacobian vector products of the SA-BCM model against complex step.
    """

    N_PROCS = 2

    h = 1e-40

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.CFDSolver = ADFLOW_C(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        setTestState(self.CFDSolver, self.ap, self.nuTildeScale)

    def cmplx_test_wDot(self):
        if not hasattr(self, "name"):
            return

        wDot = self.CFDSolver.getStatePerturbation(314)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        rtol, atol = getTol(tol=5e-9)

        self.handler.root_print("||dR/dw * wDot||")
        self.handler.par_add_norm("||dR/dw * wDot||", resDot, rtol=rtol, atol=atol)

        self.handler.root_print("dFuncs/dw * wDot")
        self.handler.root_add_dict("dFuncs/dw * wDot", funcsDot, rtol=rtol, atol=atol)

        self.handler.root_print("||dF/dw * wDot||")
        self.handler.par_add_norm("||dF/dw * wDot||", fDot, rtol=rtol, atol=atol)

    def cmplx_test_xVDot(self):
        if not hasattr(self, "name"):
            return

        xVDot = self.CFDSolver.getSpatialPerturbation(314)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode="CS", h=self.h
        )

        rtol, atol = getTol(tol=1e-10)

        self.handler.root_print("||dR/dXv * xVDot||")
        self.handler.par_add_norm("||dR/dXv * xVDot||", resDot, rtol=rtol, atol=atol)

        # These can be finiky sometimes so a bigger tolerance.
        self.handler.root_print("dFuncs/dXv * xVDot")
        self.handler.root_add_dict("dFuncs/dXv * xVDot", funcsDot, rtol=rtol * 10, atol=atol * 10)

        self.handler.root_print("||dF/dXv * xVDot||")
        self.handler.par_add_norm("||dF/dXv * xVDot||", fDot, rtol=rtol, atol=atol)


@parameterized_class(test_params)
class TestSABCMJacVecBwd(unittest.TestCase):
    """
    Tests that the reverse-mode jacobian vector products of the SA-BCM model are
    consistent with the forward mode (dot-product test).
    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            return

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)

        # propagates the values from the restart file throughout the code
        setTestState(self.CFDSolver, self.ap, self.nuTildeScale)

    def test_dot_product(self):
        wDot = self.CFDSolver.getStatePerturbation(314)
        resBar = self.CFDSolver.getStatePerturbation(271)

        resDot = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
        wBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=resBar, wDeriv=True)
        wBarFast = self.CFDSolver.computeJacobianVectorProductBwdFast(resBar=resBar)

        fwdProd = self.CFDSolver.comm.allreduce(np.dot(resBar, resDot), op=MPI.SUM)
        bwdProd = self.CFDSolver.comm.allreduce(np.dot(wBar, wDot), op=MPI.SUM)
        bwdFastProd = self.CFDSolver.comm.allreduce(np.dot(wBarFast, wDot), op=MPI.SUM)

        np.testing.assert_allclose(bwdProd, fwdProd, rtol=1e-10, err_msg="BWD dot product")
        np.testing.assert_allclose(bwdFastProd, fwdProd, rtol=1e-10, err_msg="BWDFast dot product")


if __name__ == "__main__":
    unittest.main()
