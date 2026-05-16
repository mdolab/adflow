"""Regression tests for HPDDM adjoint solver types.

Tests that GCRODR, BGMRES, and BGCRODR produce sensitivities matching GMRES
within tolerance, and that avoidRedundantAdjoints correctly handles rank-
deficient RHS vectors (e.g. clp+clv=cl).

Uses naca0012_rans-L2.cgns throughout — large enough for adjointAMGLevels=2
and gives non-trivial sensitivities for all tested functions.
"""
import os
import sys
import unittest
from itertools import product

import numpy as np
from mpi4py import MPI
from parameterized import parameterized_class

from adflow import ADFLOW
from baseclasses import AeroProblem

baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, baseDir)
from reg_default_options import adflowDefOpts

gridFile = os.path.join(baseDir, "../../input_files/naca0012_rans-L2.cgns")

HPDDM_SOLVERS = ["GCRODR", "BGMRES", "BGCRODR"]
PRECONDITIONERS = ["additive Schwarz", "multigrid"]


def _makeSolver(adjointSolver, extra=None):
    opts = dict(adflowDefOpts)
    opts.update(
        {
            "gridFile": gridFile,
            "outputDirectory": os.path.join(baseDir, adflowDefOpts["outputDirectory"]),
            "equationType": "RANS",
            "MGCycle": "sg",
            "L2Convergence": 1e-10,
            "adjointL2Convergence": 1e-8,
            "adjointSolver": adjointSolver,
            "adjointSubspaceSize": 100,
            "adjointAMGLevels": 2,
            "useANKSolver": True,
            "useNKSolver": True,
            "ANKSwitchTol": 10.0,
            "ANKSecondOrdSwitchTol": 1e-3,
            "ANKCoupledSwitchTol": 1e-6,
            "NKSwitchTol": 1e-8,
            "useBlockettes": False,
            "writeSurfaceSolution": False,
            "writeVolumeSolution": False,
        }
    )
    if extra:
        opts.update(extra)
    return ADFLOW(options=opts, comm=MPI.COMM_WORLD)


def _makeAP(evalFuncs=None):
    if evalFuncs is None:
        evalFuncs = ["cl", "cd", "cmz"]
    return AeroProblem(
        name="naca0012",
        alpha=2.77,
        mach=0.6,
        reynolds=4800000.0,
        reynoldsLength=1.0,
        T=280.0,
        R=287.085,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=evalFuncs,
    )


@parameterized_class([{"adjointSolver": s} for s in HPDDM_SOLVERS])
class TestSolverEquivalence(unittest.TestCase):
    """Verify each HPDDM solver type produces sensitivities matching GMRES.

    Each test runs its own GMRES reference because testflo's MPI runner
    executes each test in a separate process, so setUpClass state is not
    shared across test methods.
    """

    N_PROCS = 2

    def test_solver(self):
        refAP = _makeAP()
        refSolver = _makeSolver("GMRES")
        refSolver(refAP)
        refSens = {}
        refSolver.evalFunctionsSens(refAP, refSens, refAP.evalFuncs)

        ap = _makeAP()
        solver = _makeSolver(self.adjointSolver)
        solver(ap)
        sens = {}
        solver.evalFunctionsSens(ap, sens, ap.evalFuncs)
        for k in refSens:
            for dvKey in refSens[k]:
                np.testing.assert_allclose(
                    sens[k][dvKey],
                    refSens[k][dvKey],
                    atol=1e-5,
                    rtol=1e-5,
                    err_msg=f"{self.adjointSolver} mismatch on {k}/{dvKey}",
                )


@parameterized_class(
    [
        {"adjointSolver": s, "globalPreconditioner": pc}
        for s, pc in product(["GMRES"] + HPDDM_SOLVERS, PRECONDITIONERS)
    ]
)
class TestSweep(unittest.TestCase):
    """Parametric solver x preconditioner sweep with independent and rank-deficient RHS."""

    N_PROCS = 2

    def _make(self, evalFuncs, avoidRedundant=False):
        return _makeSolver(
            self.adjointSolver,
            extra={
                "globalPreconditioner": self.globalPreconditioner,
                "avoidRedundantAdjoints": avoidRedundant,
            },
        )

    def test_independent_rhs(self):
        evalFuncs = ["cl", "cd", "cmz"]
        ap = _makeAP(evalFuncs)
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
        ap_ref = _makeAP(evalFuncs)
        ref = _makeSolver("GMRES", extra={"avoidRedundantAdjoints": False})
        ref(ap_ref)
        refSens = {}
        ref.evalFunctionsSens(ap_ref, refSens, evalFuncs)

        ap = _makeAP(evalFuncs)
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
