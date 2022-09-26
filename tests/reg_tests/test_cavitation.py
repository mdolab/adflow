# built-ins
import unittest
import numpy as np
import os
import copy

# MACH classes
from adflow import ADFLOW
from adflow import ADFLOW_C

# import the testing utilities that live a few directories up
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_naca0012_cavitation
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))


class CavitationBasicTests(reg_test_classes.RegTest):
    """
    Tests for the cavitation metrics.
    """

    N_PROCS = 2
    ref_file = "cavitation_tests.json"

    def setUp(self):
        super().setUp()

        gridFile = os.path.join(baseDir, "../../input_files/naca0012_rans-L2.cgns")

        options = copy.copy(adflowDefOpts)
        options.update(
            {
                "gridfile": gridFile,
                "outputdirectory": os.path.join(baseDir, "../output_files"),
                "writevolumesolution": False,
                "writesurfacesolution": False,
                "writetecplotsurfacesolution": False,
                "mgcycle": "sg",
                "ncycles": 1000,
                "useanksolver": True,
                "usenksolver": True,
                "anksecondordswitchtol": 1e-2,
                "nkswitchtol": 1e-6,
                "volumevariables": ["temp", "mach", "resrho", "cp"],
                "equationType": "RANS",
                "l2convergence": 1e-13,
                "adjointl2convergence": 1e-13,
                "computeCavitation": True,
                "cavitationNumber": 1.0,
                "cpminrho": 1e3,
            }
        )

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca0012_cavitation)
        # change the name
        self.ap.name = "naca0012_rans_2D"

        # add aoa DV
        self.ap.addDV("alpha", name="alpha")

        CFDSolver = ADFLOW(options=options)

        self.CFDSolver = CFDSolver

    def test_cavitation_metrics_and_derivatives(self):
        "Tests if the correct amount of momentum is added to the flow by the actuator"

        evalFuncs = ["cavitation", "cpmin"]

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        ################
        # TEST ALL FUNCS
        ################
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=evalFuncs)
        self.handler.root_print("cavitation outputs")
        self.handler.root_add_dict("cavitation outputs", funcs, rtol=1e-10, atol=1e-10)

        #############
        # TEST TOTALS
        #############
        funcsSens = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSens, evalFuncs=evalFuncs)
        self.handler.root_print("cavitation totals")
        self.handler.root_add_dict("cavitation totals", funcsSens, rtol=1e-10, atol=1e-10)

        for func_name in evalFuncs:
            ##################
            # GEOMETRIC TOTALS
            ##################

            # random state and volume coordinate perturbations for tests below
            wDot = self.CFDSolver.getStatePerturbation(314)
            xVDot = self.CFDSolver.getSpatialPerturbation(314)

            # Now that we have the adjoint solution for both functionals
            # also get a total derivative wrt a random spatial perturbation DV
            psi = -self.CFDSolver.getAdjoint(func_name)
            funcsBar = self.CFDSolver._getFuncsBar(func_name)

            # this is the reverse seed up to the volume coordinates
            xVBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=psi, funcsBar=funcsBar, xVDeriv=True)

            # dot product these two vectors to get a total derivative
            dotLocal = np.dot(xVDot, xVBar)
            self.handler.par_add_sum(
                f"total {func_name} derivative wrt random volume perturbation", dotLocal, rtol=1e-10
            )

            ##################
            # DOT PRODUCT TEST
            ##################
            fDot_w = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, funcDeriv=True)
            fDot_xv = self.CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, funcDeriv=True)
            funcsBar = {func_name: 1.0}
            wBar, xVBar = self.CFDSolver.computeJacobianVectorProductBwd(funcsBar=funcsBar, wDeriv=True, xVDeriv=True)

            # do the dot product. we test both state partials and volume coordinate partials
            # state
            dotLocal1 = np.dot(wDot, wBar)
            dotLocal2 = fDot_w[func_name] / self.CFDSolver.comm.size

            self.handler.par_add_sum(f"Dot product test for w -> {func_name}", dotLocal1, rtol=1e-10)
            self.handler.par_add_sum(f"Dot product test for w -> {func_name}", dotLocal2, rtol=1e-10, compare=True)

            # volume coords
            dotLocal1 = np.dot(xVDot, xVBar)
            dotLocal2 = fDot_xv[func_name] / self.CFDSolver.comm.size
            self.handler.par_add_sum(f"Dot product test for xV -> {func_name}", dotLocal1, rtol=1e-10)
            self.handler.par_add_sum(f"Dot product test for xV -> {func_name}", dotLocal2, rtol=1e-10, compare=True)


class CavitationCmplxTests(reg_test_classes.CmplxRegTest):
    """
    Complex step tests for the cavitation functions.
    """

    N_PROCS = 2
    ref_file = "cavitation_tests.json"
    h = 1e-40

    def setUp(self):
        super().setUp()

        gridFile = os.path.join(baseDir, "../../input_files/naca0012_rans-L2.cgns")

        options = copy.copy(adflowDefOpts)
        options.update(
            {
                "gridfile": gridFile,
                "outputdirectory": os.path.join(baseDir, "../output_files"),
                "writevolumesolution": False,
                "writesurfacesolution": False,
                "writetecplotsurfacesolution": False,
                "mgcycle": "sg",
                "ncycles": 1000,
                "useanksolver": True,
                "usenksolver": True,
                "anksecondordswitchtol": 1e-2,
                "nkswitchtol": 1e-6,
                "volumevariables": ["temp", "mach", "resrho", "cp"],
                "equationType": "RANS",
                "l2convergence": 1e-13,
                "adjointl2convergence": 1e-13,
                "computeCavitation": True,
                "cavitationNumber": 1.0,
                "NKUseEW": False,
                "NKLinearSolveTol": 1e-4,
            }
        )

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca0012_cavitation)
        # change the name
        self.ap.name = "naca0012_rans_2D"

        # add aoa DV
        self.ap.addDV("alpha", name="alpha")

        CFDSolver = ADFLOW_C(options=options, debug=True)

        self.CFDSolver = CFDSolver

    def cmplx_test_cavitation_adjoints(self):

        aDV = {"alpha": self.ap.alpha}
        funcs = {}
        funcsSensCS = {}

        # run the baseline evaluation
        self.CFDSolver(self.ap)
        # check if solution failed
        self.assert_solution_failure()
        self.CFDSolver.evalFunctions(self.ap, funcs)

        funcs_plus = {}
        for dv in ["alpha", "vol_perturbation"]:
            if dv == "alpha":
                # save the old alpha
                dvsave = aDV["alpha"]
                # perturb
                aDV["alpha"] += self.h * 1j
                self.ap.setDesignVars(aDV)
            elif dv == "vol_perturbation":
                xVDot = self.CFDSolver.getSpatialPerturbation(314)
                # grid_save = np.zeros_like(xVDot)
                # get the original grid
                grid_save = self.CFDSolver.adflow.warping.getgrid(self.CFDSolver.getSpatialSize())
                # perturb using the random seed
                # this is very intrusive but it works and its testing these functions sensitivities
                # wrt geometric DVs w/o actually having a mesh or dvgeo object
                self.CFDSolver.adflow.warping.setgrid(grid_save + self.h * 1j * xVDot)
                self.CFDSolver.adflow.preprocessingapi.updatecoordinatesalllevels()
                self.CFDSolver.adflow.walldistance.updatewalldistancealllevels()
                self.CFDSolver.adflow.preprocessingapi.updatemetricsalllevels()
                self.CFDSolver.adflow.preprocessingapi.updategridvelocitiesalllevels()

            # call solver again
            # we can also not re-set the flow and call the solver a few times until complex parts converge
            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver(self.ap)  # the complex residuals dont converge well for the coordinate perturbation

            # check if solution failed
            self.assert_solution_failure()

            # save the new funcs in a dict
            funcs_plus[dv] = {}

            # eval functions
            self.CFDSolver.evalFunctions(self.ap, funcs_plus[dv])

            # compute the sens
            funcsSensCS[dv] = {}
            for f in ["cavitation"]:
                funcsSensCS[dv][self.ap[f]] = np.imag(funcs_plus[dv][self.ap[f]]) / self.h

            # reset the DV
            if dv == "alpha":
                aDV["alpha"] = dvsave
                self.ap.setDesignVars(aDV)
            elif dv == "vol_perturbation":
                # set the original grid back
                self.CFDSolver.adflow.warping.setgrid(grid_save)

        ##################
        # TEST DERIVATIVES
        ##################

        # here we compare the CS derivatives to the adjoint values computed by the real test
        # we treat the CS value as the truth, so if this test passes,
        # we assume the adjoint sensitivities are also true

        func_name = "naca0012_rans_2D_cavitation"

        ref_val = self.handler.db["cavitation totals"][func_name]["alpha"]
        np.testing.assert_allclose(funcsSensCS["alpha"][func_name], ref_val, atol=1e-10, rtol=1e-10)

        ref_val = self.handler.db["total cavitation derivative wrt random volume perturbation"]
        np.testing.assert_allclose(funcsSensCS["vol_perturbation"][func_name], ref_val, rtol=1e-4)


if __name__ == "__main__":
    unittest.main()
