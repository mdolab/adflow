"""Regression tests for HPDDM adjoint solver types.

Tests that GCRODR, BGMRES, and BGCRODR produce sensitivities matching GMRES
within tolerance, and that avoidRedundantAdjoints correctly handles rank-
deficient RHS vectors (e.g. clp+clv=cl).
"""
import os
import sys
import unittest

import numpy as np
from mpi4py import MPI
from parameterized import parameterized_class

from adflow import ADFLOW
from baseclasses import AeroProblem

# Allow running from any directory by inserting the reg_tests dir into path
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, baseDir)
from reg_default_options import adflowDefOpts

gridFile = os.path.join(baseDir, "../../input_files/cube_4x4x4.cgns")


def _makeSolver(adjointSolver, extra=None):
    opts = dict(adflowDefOpts)
    opts.update(
        {
            "gridFile": gridFile,
            "MGCycle": "sg",
            "L2Convergence": 1e-12,
            "adjointL2Convergence": 1e-10,
            "adjointSolver": adjointSolver,
            "adjointSubspaceSize": 100,
            # cube_4x4x4 has only 64 cells; 2-level AMG coarsens to 8 cells which
            # is smaller than the 14-point stencil → MatCreateBAIJ error. Use 1 level.
            "adjointAMGLevels": 1,
            "writeSurfaceSolution": False,
            "writeVolumeSolution": False,
        }
    )
    if extra:
        opts.update(extra)
    return ADFLOW(options=opts, comm=MPI.COMM_WORLD)


def _makeAP():
    return AeroProblem(
        name="cube",
        alpha=1.0,
        mach=0.5,
        altitude=1000.0,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=["cl", "cd", "cmz"],
    )


class TestSolverEquivalence(unittest.TestCase):
    """Verify HPDDM solver types produce sensitivities matching GMRES.

    Each test runs its own GMRES reference because testflo's MPI runner
    executes each test in a separate process, so setUpClass state is not
    shared across test methods.
    """

    N_PROCS = 1

    def _compare(self, adjointSolver, atol=1e-5, rtol=1e-5, extra=None):
        refAP = _makeAP()
        refSolver = _makeSolver("GMRES")
        refSolver(refAP)
        refSens = {}
        refSolver.evalFunctionsSens(refAP, refSens, refAP.evalFuncs)

        ap = _makeAP()
        solver = _makeSolver(adjointSolver, extra)
        solver(ap)
        sens = {}
        solver.evalFunctionsSens(ap, sens, ap.evalFuncs)
        for k in refSens:
            for dvKey in refSens[k]:
                np.testing.assert_allclose(
                    sens[k][dvKey],
                    refSens[k][dvKey],
                    atol=atol,
                    rtol=rtol,
                    err_msg=f"{adjointSolver} mismatch on {k}/{dvKey}",
                )

    def test_gcrodr(self):
        self._compare("GCRODR")

    def test_bgmres(self):
        self._compare("BGMRES")

    def test_bgcrodr(self):
        self._compare("BGCRODR")


@parameterized_class(
    [
        {"adjointSolver": "GMRES", "globalPreconditioner": "additive Schwarz"},
        {"adjointSolver": "GCRODR", "globalPreconditioner": "additive Schwarz"},
        {"adjointSolver": "BGMRES", "globalPreconditioner": "additive Schwarz"},
        {"adjointSolver": "BGCRODR", "globalPreconditioner": "additive Schwarz"},
    ]
)
class TestSweep(unittest.TestCase):
    """Parametric solver x preconditioner sweep with independent and rank-deficient RHS."""

    N_PROCS = 1

    def _make(self, evalFuncs, avoidRedundant=False):
        return _makeSolver(
            self.adjointSolver,
            extra={
                "globalPreconditioner": self.globalPreconditioner,
                "avoidRedundantAdjoints": avoidRedundant,
            },
        )

    def _makeAP(self, evalFuncs):
        return AeroProblem(
            name="cube",
            alpha=1.0,
            mach=0.5,
            altitude=1000.0,
            areaRef=1.0,
            chordRef=1.0,
            evalFuncs=evalFuncs,
        )

    def test_independent_rhs(self):
        evalFuncs = ["cl", "cd", "cmz"]
        ap = self._makeAP(evalFuncs)
        solver = self._make(evalFuncs)
        solver(ap)
        sens = {}
        solver.evalFunctionsSens(ap, sens, evalFuncs)
        self.assertFalse(ap.adjointFailed)
        for f in evalFuncs:
            key = ap.name + "_" + f
            self.assertIn(key, sens)

    def test_rank_deficient_rhs(self):
        # clp + clv = cl, cdp + cdv = cd — 2 exact linear dependencies.
        evalFuncs = ["clp", "clv", "cl", "cdp", "cdv", "cd", "cmz"]
        ap_ref = self._makeAP(evalFuncs)
        ref = _makeSolver("GMRES", extra={"avoidRedundantAdjoints": False})
        ref(ap_ref)
        refSens = {}
        ref.evalFunctionsSens(ap_ref, refSens, evalFuncs)

        ap = self._makeAP(evalFuncs)
        solver = self._make(evalFuncs, avoidRedundant=True)
        solver(ap)
        sens = {}
        solver.evalFunctionsSens(ap, sens, evalFuncs)

        for k in refSens:
            for dvKey, dvVal in refSens[k].items():
                np.testing.assert_allclose(
                    sens[k][dvKey],
                    dvVal,
                    atol=1e-5,
                    rtol=1e-5,
                    err_msg=f"{self.adjointSolver}+{self.globalPreconditioner} mismatch on {k}/{dvKey}",
                )


if __name__ == "__main__":
    unittest.main()
