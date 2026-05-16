"""Helpers for ADflow's adjoint linear-algebra path."""
import numpy as np
from scipy.linalg import solve_triangular
from petsc4py import PETSc

_HPDDM_SOLVERS = {"gcrodr", "bgmres", "bgcrodr"}


def isHPDDMAvailable():
    """Probe whether PETSc was built with HPDDM by trying to set the type."""
    try:
        ksp = PETSc.KSP().create(comm=PETSc.COMM_SELF)
        ksp.setType("hpddm")
        ksp.destroy()
        return True
    except PETSc.Error:
        return False


def adjointSolverNeedsHPDDM(solverName):
    return solverName.lower() in _HPDDM_SOLVERS


def wrapNumpyAsVec(arr, comm=None):
    """Wrap a 1-D numpy array as a PETSc Vec sharing memory (no copy)."""
    if comm is None:
        comm = PETSc.COMM_WORLD
    return PETSc.Vec().createWithArray(arr, comm=comm)


def stackRhsAsMat(vecs, comm=None):
    """Stack a list of PETSc Vecs as columns of a dense Mat (helper for tests).
    Direct port of LinearSolverPlayground/.../loaders.py:stackRhsAsMat.
    """
    if comm is None:
        comm = PETSc.COMM_WORLD
    nRhs = len(vecs)
    nLocal = vecs[0].getLocalSize()
    nGlobal = vecs[0].getSize()
    B = PETSc.Mat().createDense(((nLocal, nGlobal), (nRhs, nRhs)), comm=comm)
    B.setUp()
    for jj, v in enumerate(vecs):
        col = B.getDenseColumnVec(jj, mode="w")
        v.copy(col)
        B.restoreDenseColumnVec(jj, mode="w")
    B.assemblyBegin()
    B.assemblyEnd()
    return B


def pivotedQR(vecs, rtol=1e-12):
    """Distributed column-pivoted modified Gram-Schmidt QR.

    Direct port of LinearSolverPlayground/.../benchmark.py:_pivotedQR.
    Returns (indep_idx, dep_idx, C, Q, R_basis).
    - indep_idx: indices of linearly independent input vecs
    - dep_idx:   indices of dependent input vecs
    - C:         (len(indep), len(dep)) reconstruction coefficients s.t.
                 vecs[dep[j]] ≈ sum_k C[k,j] * vecs[indep[k]]
    - Q:         orthonormal basis Vecs (caller must destroy)
    - R_basis:   (rank, k) upper-triangular factor in pivot order
    """
    k = len(vecs)
    if k == 0:
        return [], [], None, [], np.zeros((0, 0))

    R = np.zeros((k, k))
    perm = list(range(k))
    Q = []
    origNorm2 = np.array([float(v.dot(v)) for v in vecs])
    res2 = origNorm2.copy()
    rank = k
    firstPivotNorm = 0.0

    for step in range(k):
        best = step + int(np.argmax(res2[step:]))
        if best != step:
            perm[step], perm[best] = perm[best], perm[step]
            res2[step], res2[best] = res2[best], res2[step]

        pivotNorm = np.sqrt(max(res2[step], 0.0))
        if step == 0:
            firstPivotNorm = pivotNorm
            if firstPivotNorm == 0.0:
                rank = 0
                break
        if pivotNorm < rtol * firstPivotNorm:
            rank = step
            break

        w = vecs[perm[step]].copy()
        for ll in range(step):
            r_ll = float(Q[ll].dot(w))
            R[ll, step] = r_ll
            w.axpy(-r_ll, Q[ll])
        for ll in range(step):
            corr = float(Q[ll].dot(w))
            R[ll, step] += corr
            w.axpy(-corr, Q[ll])
        R[step, step] = float(w.norm())
        w.scale(1.0 / R[step, step])
        Q.append(w)

        for jj in range(step + 1, k):
            r_sj = float(Q[step].dot(vecs[perm[jj]]))
            newRes2 = res2[jj] - r_sj * r_sj
            if newRes2 < 0.1 * origNorm2[perm[jj]]:
                wTmp = vecs[perm[jj]].copy()
                for ll in range(step + 1):
                    wTmp.axpy(-float(Q[ll].dot(wTmp)), Q[ll])
                newRes2 = max(float(wTmp.dot(wTmp)), 0.0)
                wTmp.destroy()
            res2[jj] = newRes2

    for jj in range(rank, k):
        for ll in range(rank):
            R[ll, jj] = float(Q[ll].dot(vecs[perm[jj]]))

    indepIdx = sorted(perm[:rank])
    depIdx = sorted(perm[rank:])
    rowOrder = [perm[:rank].index(ii) for ii in indepIdx]
    R_basis = np.empty((rank, k))
    R_basis[:, :rank] = R[:rank, :rank][:, rowOrder]

    if depIdx:
        colOrder = [perm[rank:].index(ii) for ii in depIdx]
        R_basis[:, rank:] = R[:rank, rank:k][:, colOrder]
        C_perm = solve_triangular(R[:rank, :rank], R[:rank, rank:k], lower=False)
        C = C_perm[np.ix_(rowOrder, colOrder)]
    else:
        C = None

    return indepIdx, depIdx, C, Q, R_basis
