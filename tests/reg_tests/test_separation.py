# built-ins
import unittest
import numpy as np
import os
import copy

# MACH classes
from adflow import ADFLOW

# import the testing utilities
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_naca0012_separation
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))


class SeparationBasicTests(reg_test_classes.RegTest):
    """
    Tests for the separation metrics.
    """

    N_PROCS = 2
    ref_file = "separation_tests.json"

    def setUp(self):
        super().setUp()

        gridFile = os.path.join(baseDir, "../../input_files/naca0012_L3_SEP.cgns")

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
                "l2convergence": 1e-15,
                "adjointl2convergence": 1e-15,
                "sepSensorModel": "surfvec",
            }
        )

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca0012_separation)
        # change the name
        self.ap.name = "naca0012_rans_2D"

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


if __name__ == "__main__":
    unittest.main()
