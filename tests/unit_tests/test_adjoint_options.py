import os
import sys
import unittest

baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(baseDir, "../.."))

from adflow import ADFLOW


class TestAdjointSolverOptions(unittest.TestCase):
    def setUp(self):
        self.opts = {"gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler.cgns")}

    def test_gcrodr_accepted(self):
        self.opts["adjointSolver"] = "GCRODR"
        ADFLOW(options=self.opts)  # Should not raise on enum validation

    def test_bgmres_accepted(self):
        self.opts["adjointSolver"] = "BGMRES"
        ADFLOW(options=self.opts)

    def test_bgcrodr_accepted(self):
        self.opts["adjointSolver"] = "BGCRODR"
        ADFLOW(options=self.opts)

    def test_recycle_size_default(self):
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("adjointRecycleSize"), 5)

    def test_avoid_redundant_default_false(self):
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("avoidRedundantAdjoints"), False)

    def test_deflation_tol_default_zero(self):
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("adjointHPDDMDeflationTol"), 0.0)

    def test_recycle_size_settable(self):
        self.opts["adjointRecycleSize"] = 10
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("adjointRecycleSize"), 10)

    def test_avoid_redundant_settable(self):
        self.opts["avoidRedundantAdjoints"] = True
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("avoidRedundantAdjoints"), True)

    def test_deflation_tol_settable(self):
        self.opts["adjointHPDDMDeflationTol"] = 1e-4
        solver = ADFLOW(options=self.opts)
        self.assertAlmostEqual(solver.getOption("adjointHPDDMDeflationTol"), 1e-4)


if __name__ == "__main__":
    unittest.main()
