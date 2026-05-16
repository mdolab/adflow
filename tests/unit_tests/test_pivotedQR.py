import os
import sys
import unittest

import numpy as np
from petsc4py import PETSc

baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(baseDir, "../.."))

from adflow.adjointLinAlg import pivotedQR, wrapNumpyAsVec


class TestPivotedQR(unittest.TestCase):
    def _vecsFromColumns(self, M):
        return [wrapNumpyAsVec(M[:, jj].copy()) for jj in range(M.shape[1])]

    def test_full_rank_random(self):
        rng = np.random.default_rng(0)
        M = rng.standard_normal((20, 3))
        vecs = self._vecsFromColumns(M)
        indep, dep, C, Q, _ = pivotedQR(vecs, rtol=1e-12)
        self.assertEqual(sorted(indep), [0, 1, 2])
        self.assertEqual(dep, [])
        self.assertIsNone(C)
        for q in Q:
            q.destroy()
        for v in vecs:
            v.destroy()

    def test_one_exact_dependency(self):
        rng = np.random.default_rng(1)
        M = rng.standard_normal((20, 2))
        M_dep = np.concatenate([M, (2 * M[:, 0:1] - M[:, 1:2])], axis=1)
        vecs = self._vecsFromColumns(M_dep)
        indep, dep, C, Q, _ = pivotedQR(vecs, rtol=1e-10)
        self.assertEqual(len(indep), 2)
        self.assertEqual(len(dep), 1)
        recon = C[0, 0] * M_dep[:, indep[0]] + C[1, 0] * M_dep[:, indep[1]]
        np.testing.assert_allclose(recon, M_dep[:, dep[0]], atol=1e-10)
        for q in Q:
            q.destroy()
        for v in vecs:
            v.destroy()

    def test_wind_body_frame(self):
        # CL, CD independent. CFx = cos(a)*CD - sin(a)*CL; CFz = sin(a)*CD + cos(a)*CL.
        rng = np.random.default_rng(2)
        cl = rng.standard_normal(30)
        cd = rng.standard_normal(30)
        a = 0.1
        cfx = np.cos(a) * cd - np.sin(a) * cl
        cfz = np.sin(a) * cd + np.cos(a) * cl
        M = np.column_stack([cl, cd, cfx, cfz])
        vecs = self._vecsFromColumns(M)
        indep, dep, C, Q, _ = pivotedQR(vecs, rtol=1e-12)
        self.assertEqual(len(indep), 2)
        self.assertEqual(len(dep), 2)
        # Verify reconstruction accuracy for each dependent column
        for jj, depIdx in enumerate(dep):
            recon = sum(C[kk, jj] * M[:, indep[kk]] for kk in range(len(indep)))
            np.testing.assert_allclose(recon, M[:, depIdx], atol=1e-10)
        for q in Q:
            q.destroy()
        for v in vecs:
            v.destroy()

    def test_threshold_at_rtol(self):
        # v1 is unit-ish in the x direction, v2 is orthogonal to v1 with
        # norm 1e-8 * ||v1||. After projecting out v1, the residual of v2
        # is 1e-8 * ||v1||. So at rtol=1e-6 (1e-8 < 1e-6) → v2 classified
        # as dependent; at rtol=1e-12 (1e-8 > 1e-12) → v2 classified as
        # independent.
        n = 20
        v1 = np.zeros(n)
        v1[0] = 1.0  # unit vector in x
        v2 = np.zeros(n)
        v2[1] = 1e-8  # orthogonal to v1, norm = 1e-8
        M = np.column_stack([v1, v2])
        for rtol, expectedDep in [(1e-6, 1), (1e-12, 0)]:
            vecs = self._vecsFromColumns(M)
            _, dep, _, Q, _ = pivotedQR(vecs, rtol=rtol)
            self.assertEqual(len(dep), expectedDep)
            for q in Q:
                q.destroy()
            for v in vecs:
                v.destroy()

    def test_single_rhs(self):
        vecs = self._vecsFromColumns(np.ones((10, 1)))
        indep, dep, C, Q, _ = pivotedQR(vecs)
        self.assertEqual(indep, [0])
        self.assertEqual(dep, [])
        self.assertIsNone(C)
        for q in Q:
            q.destroy()
        for v in vecs:
            v.destroy()

    def test_empty_vecs(self):
        indep, dep, C, Q, R = pivotedQR([])
        self.assertEqual(indep, [])
        self.assertEqual(dep, [])
        self.assertIsNone(C)
        self.assertEqual(Q, [])
        self.assertEqual(R.shape, (0, 0))

    def test_all_zero_cols(self):
        M = np.zeros((10, 3))
        M[0, 0] = 1.0  # first col is nonzero; rest are zero
        vecs = self._vecsFromColumns(M)
        indep, dep, C, Q, _ = pivotedQR(vecs, rtol=1e-12)
        self.assertEqual(len(indep), 1)
        self.assertEqual(len(dep), 2)
        for q in Q:
            q.destroy()
        for v in vecs:
            v.destroy()


if __name__ == "__main__":
    unittest.main()
