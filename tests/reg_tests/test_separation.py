# built-ins
import unittest
import numpy as np
import os
import copy

# MACH classes
from adflow import ADFLOW
from adflow import ADFLOW_C

# import the testing utilities
from parameterized import parameterized_class
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_naca0012_separation
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))
gridFile = os.path.join(baseDir, "../../input_files/naca0012_L3_SEP.cgns")
test_params = [
     {
        "name": "naca0012_rans_2D_surfvec_ks",
        "ref_file" :"separation_tests_surfvec_ks.json",
           "options":     {
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
                "l2convergence": 1e-15,
                "adjointl2convergence": 1e-15,
                "sepSensorModel": "surfvec_ks",
                "sepSweepAngleCorrection": 0.0,
            }},
             {
        "name": "naca0012_rans_2D_surfvec",
        "ref_file" :"separation_tests_surfvec.json",
           "options":     {
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
                "l2convergence": 1e-15,
                "adjointl2convergence": 1e-15,
                "sepSensorModel": "surfvec_ks",
                "sepSweepAngleCorrection": 0.0,
            }},
            {
        "name": "naca0012_rans_2D_heaviside",
        "ref_file" :"separation_tests_heaviside.json",
           "options":     {
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
                "l2convergence": 1e-15,
                "adjointl2convergence": 1e-15,
                "sepSensorModel": "heaviside",
                "sepSweepAngleCorrection": 0.0,
            }}]

@parameterized_class(test_params)
class SeparationBasicTests(reg_test_classes.RegTest):
    """
    Tests for the separation metrics.
    """

    N_PROCS = 2
    # ref_file = None
    options = None

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return
        super().setUp()

        options = copy.copy(adflowDefOpts)
        options.update(self.options)

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca0012_separation)
        # change the name
        self.ap.name = self.name

        # add aoa DV
        self.ap.addDV("alpha", name="alpha")

        CFDSolver = ADFLOW(options=options)

        self.CFDSolver = CFDSolver
        self.CFDSolver.addFamilyGroup("wingup", ["upper_skin"])

        self.CFDSolver.addFunction("sepsensor", "wingup")

    def test_separation_metrics_and_derivatives(self):
        evalFuncs = ["sepsensor_wingup", "cl"]

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        ################
        # TEST ALL FUNCS
        ################
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=evalFuncs)
        self.handler.root_print("separation outputs")
        self.handler.root_add_dict("separation outputs", funcs, rtol=5e-10, atol=5e-10)

        #############
        # TEST TOTALS
        #############
        funcsSens = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSens, evalFuncs=evalFuncs)
        self.handler.root_print("separation totals")
        self.handler.root_add_dict("separation totals", funcsSens, rtol=5e-10, atol=5e-10)

        for funcName in evalFuncs:
            ##################
            # GEOMETRIC TOTALS
            ##################

            # random state and volume coordinate perturbations for tests below
            wDot = self.CFDSolver.getStatePerturbation(314)
            xVDot = self.CFDSolver.getSpatialPerturbation(314)

            # Now that we have the adjoint solution for both functionals
            # also get a total derivative wrt a random spatial perturbation DV
            psi = -self.CFDSolver.getAdjoint(funcName)
            funcsBar = self.CFDSolver._getFuncsBar(funcName)

            # this is the reverse seed up to the volume coordinates
            xVBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=psi, funcsBar=funcsBar, xVDeriv=True)

            # dot product these two vectors to get a total derivative
            dotLocal = np.dot(xVDot, xVBar)
            # this is not the best test out there; the final answer does get affected quite a bit
            # by the processor count or architecture. We just have it here to have some test on the
            # separation functionals' derivatives w.r.t. spatial changes in a lazy way.
            self.handler.par_add_sum(f"total {funcName} derivative wrt random volume perturbation", dotLocal, rtol=1e-3)

            ##################
            # DOT PRODUCT TEST
            ##################
            fDot_w = self.CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, funcDeriv=True)
            fDot_xv = self.CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, funcDeriv=True)
            funcsBar = {funcName: 1.0}
            wBar, xVBar = self.CFDSolver.computeJacobianVectorProductBwd(funcsBar=funcsBar, wDeriv=True, xVDeriv=True)

            # do the dot product. we test both state partials and volume coordinate partials
            # state
            dotLocal1 = np.dot(wDot, wBar)
            dotLocal2 = fDot_w[funcName] / self.CFDSolver.comm.size

            self.handler.par_add_sum(f"Dot product test for w -> {funcName}", dotLocal1, rtol=5e-10)
            self.handler.par_add_sum(f"Dot product test for w -> {funcName}", dotLocal2, rtol=5e-10, compare=True)

            # volume coords
            dotLocal1 = np.dot(xVDot, xVBar)
            dotLocal2 = fDot_xv[funcName] / self.CFDSolver.comm.size
            self.handler.par_add_sum(f"Dot product test for xV -> {funcName}", dotLocal1, rtol=5e-10)
            self.handler.par_add_sum(f"Dot product test for xV -> {funcName}", dotLocal2, rtol=5e-10, compare=True)


@parameterized_class(test_params)
class SeparationCmplxTests(reg_test_classes.CmplxRegTest):
    """
    Complex step tests for the separation functions.
    """

    N_PROCS = 2
    ref_file = None
    h = 1e-40

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return
        super().setUp()


        options = copy.copy(adflowDefOpts)
        options.update(self.options)

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca0012_separation)
        # change the name
        self.ap.name = self.name

        # add aoa DV
        self.ap.addDV("alpha", name="alpha")

        CFDSolver = ADFLOW_C(options=options, debug=True)

        self.CFDSolver = CFDSolver

        self.CFDSolver.addFamilyGroup("wingup", ["upper_skin"])

        self.CFDSolver.addFunction("sepsensor", "wingup")

    def cmplx_test_separation_adjoints(self):
        aDV = {"alpha": self.ap.alpha}
        funcs = {}
        funcsSensCS = {}

        # run the baseline evaluation
        self.CFDSolver(self.ap)
        # check if solution failed
        self.assert_solution_failure()
        self.CFDSolver.evalFunctions(self.ap, funcs)

        evalFuncs = ["sepsensor_wingup", "cl"]
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
                # gridSave = np.zeros_like(xVDot)
                # get the original grid
                gridSave = self.CFDSolver.adflow.warping.getgrid(self.CFDSolver.getSpatialSize())
                # perturb using the random seed
                # this is very intrusive but it works and its testing these functions sensitivities
                # wrt geometric DVs w/o actually having a mesh or dvgeo object
                self.CFDSolver.adflow.warping.setgrid(gridSave + self.h * 1j * xVDot)
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
            for f in evalFuncs:
                funcsSensCS[dv][self.ap[f]] = np.imag(funcs_plus[dv][self.ap[f]]) / self.h

            # reset the DV
            if dv == "alpha":
                aDV["alpha"] = dvsave
                self.ap.setDesignVars(aDV)
            elif dv == "vol_perturbation":
                # set the original grid back
                self.CFDSolver.adflow.warping.setgrid(gridSave)

        ##################
        # TEST DERIVATIVES
        ##################

        # here we compare the CS derivatives to the adjoint values computed by the real test
        # we treat the CS value as the truth, so if this test passes,
        # we assume the adjoint sensitivities are also true

        for funcName in evalFuncs:
            fullName = f"naca0012_rans_2D_{funcName}"

            refVal = self.handler.db["separation totals"][fullName]["alpha"]
            np.testing.assert_allclose(funcsSensCS["alpha"][fullName], refVal, atol=1e-6, rtol=1e-6)

            refVal = self.handler.db[f"total {funcName} derivative wrt random volume perturbation"]
            np.testing.assert_allclose(funcsSensCS["vol_perturbation"][fullName], refVal, rtol=1e-3, atol=1e-3)


if __name__ == "__main__":
    unittest.main()
