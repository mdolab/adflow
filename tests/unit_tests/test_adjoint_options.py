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
        """Given adjointSolver set to 'GCRODR', if ADFLOW is created, then no
        option-validation error is raised.
        """
        self.opts["adjointSolver"] = "GCRODR"
        ADFLOW(options=self.opts)

    def test_bgmres_accepted(self):
        """Given adjointSolver set to 'BGMRES', if ADFLOW is created, then no
        option-validation error is raised.
        """
        self.opts["adjointSolver"] = "BGMRES"
        ADFLOW(options=self.opts)

    def test_bgcrodr_accepted(self):
        """Given adjointSolver set to 'BGCRODR', if ADFLOW is created, then no
        option-validation error is raised.
        """
        self.opts["adjointSolver"] = "BGCRODR"
        ADFLOW(options=self.opts)

    def test_recycle_size_default(self):
        """Given an ADFLOW solver created without specifying adjointRecycleSize,
        if the option is queried, then it returns the default value of 5.
        """
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("adjointRecycleSize"), 5)

    def test_avoid_redundant_default_false(self):
        """Given an ADFLOW solver created without specifying avoidRedundantAdjoints,
        if the option is queried, then it returns False.
        """
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("avoidRedundantAdjoints"), False)

    def test_deflation_tol_default_zero(self):
        """Given an ADFLOW solver created without specifying adjointHPDDMDeflationTol,
        if the option is queried, then it returns 0.0.
        """
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("adjointHPDDMDeflationTol"), 0.0)

    def test_recycle_size_settable(self):
        """Given adjointRecycleSize set to 10, if an ADFLOW solver is created and
        the option is queried, then it returns 10.
        """
        self.opts["adjointRecycleSize"] = 10
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("adjointRecycleSize"), 10)

    def test_avoid_redundant_settable(self):
        """Given avoidRedundantAdjoints set to True, if an ADFLOW solver is created
        and the option is queried, then it returns True.
        """
        self.opts["avoidRedundantAdjoints"] = True
        solver = ADFLOW(options=self.opts)
        self.assertEqual(solver.getOption("avoidRedundantAdjoints"), True)

    def test_deflation_tol_settable(self):
        """Given adjointHPDDMDeflationTol set to 1e-4, if an ADFLOW solver is created
        and the option is queried, then it returns 1e-4.
        """
        self.opts["adjointHPDDMDeflationTol"] = 1e-4
        solver = ADFLOW(options=self.opts)
        self.assertAlmostEqual(solver.getOption("adjointHPDDMDeflationTol"), 1e-4)

    def test_hpddm_unavailable_raises(self):
        """Given a PETSc build without HPDDM support, if ADFLOW is created with
        an HPDDM adjointSolver, then an exception is raised whose message names HPDDM.
        """
        from adflow.adjointLinAlg import isHPDDMAvailable

        if isHPDDMAvailable():
            self.skipTest("HPDDM is available; cannot test the error path")
        self.opts["adjointSolver"] = "GCRODR"
        with self.assertRaises(Exception) as ctx:
            ADFLOW(options=self.opts)
        self.assertIn("HPDDM", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
