"""
Test script for time spectral adjoint verification.

Verifies the adjoint of the time spectral residual initialization
using dot product tests (forward vs backward consistency) and
forward mode AD vs finite difference comparison.

Run with: mpirun -np 2 python test_ts_adjoint.py
"""

import os
import sys
import copy
import numpy as np

# Add the test directory to path for reg_default_options
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "reg_tests"))

from adflow import ADFLOW

try:
    from idwarp import USMesh

    HAS_IDWARP = True
except ImportError:
    HAS_IDWARP = False

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

baseDir = os.path.dirname(os.path.abspath(__file__))
inputDir = os.path.join(baseDir, "../input_files")
outputDir = os.path.join(baseDir, "output_files")


def printroot(msg):
    if rank == 0:
        print(msg, flush=True)


def setup_ts_solver():
    """Set up the time spectral solver with mesh deformation."""
    from baseclasses import AeroProblem

    gridFile = os.path.join(inputDir, "naca64A010_euler-L2.cgns")

    if not os.path.isfile(gridFile):
        printroot(f"ERROR: Grid file not found: {gridFile}")
        printroot("Run: cd input_files && bash get-input-files.sh")
        sys.exit(1)

    ntimeintervalsspectral = 3

    options = {
        # Common defaults
        "gridFile": gridFile,
        "outputDirectory": outputDir,
        "writeVolumeSolution": False,
        "writeSurfaceSolution": False,
        # Solver options
        "blocksplitting": True,
        "useblockettes": False,
        "equationtype": "Euler",
        "equationmode": "time spectral",
        "mgcycle": "sg",
        "l2convergence": 1e-13,
        "ncycles": 200000,
        "monitorvariables": ["resrho", "cl"],
        # NK solver
        "usenksolver": True,
        "nkswitchtol": 1e-4,
        "NKSubSpaceSize": 400,
        "applypcsubspacesize": 400,
        # ANK solver
        "useanksolver": True,
        "ankswitchtol": 1e-2,
        "anksubspacesize": 50,
        # Time spectral
        "alphafollowing": False,
        "timeintervals": ntimeintervalsspectral,
        "useexternaldynamicmesh": True,
        "usetsinterpolatedgridvelocity": True,
        # Adjoint
        "adjointl2convergence": 1e-14,
        "adjointmaxiter": 500,
        "nkadpc": False,
    }

    # Create the solver
    CFDSolver = ADFLOW(options=options, debug=True)

    # Setup aeroproblem
    ap = AeroProblem(
        name="64A010pitchingTS",
        alpha=0.0,
        mach=0.796,
        reynolds=12.56e6,
        reynoldsLength=1.0,
        T=305.0,
        R=287.085,
        areaRef=1.0,
        chordRef=1.0,
        evalFuncs=["cl", "cd", "cmz"],
        xRef=0.248,
        xRot=0.248,
        omegaFourier=112.59,
    )
    ap.addDV("alpha")
    ap.addDV("mach")

    if HAS_IDWARP:
        # Deform the mesh for time spectral
        meshOptions = {"gridFile": gridFile}
        mesh = USMesh(options=meshOptions)
        CFDSolver.setMesh(mesh)

        # Motion history - pitching oscillation
        alpha_0 = 1.01
        deltaAlpha = -alpha_0 * np.pi / 180.0
        alpha = np.linspace(0.0, 2.0 * np.pi, ntimeintervalsspectral + 1)[:-1]
        alpha = -np.sin(alpha)
        alpha *= deltaAlpha

        xRot = 0.25

        # Get undeformed points
        MDGroup = CFDSolver.allWallsGroup
        cfdPts0 = CFDSolver.getSurfaceCoordinates(MDGroup, includeZipper=False)
        N_pts = cfdPts0.shape[0]

        # Deform mesh for each time instance
        for sps in range(ntimeintervalsspectral):
            cfdPoints_deformed = np.zeros((N_pts, 3))
            ptch_loc = alpha[sps]
            cc = np.cos(ptch_loc)
            ss = np.sin(ptch_loc)
            for j in range(N_pts):
                cfdPoints_deformed[j, 0] = cc * (cfdPts0[j, 0] - xRot) + ss * cfdPts0[j, 1] + xRot
                cfdPoints_deformed[j, 1] = -ss * (cfdPts0[j, 0] - xRot) + cc * cfdPts0[j, 1]
                cfdPoints_deformed[j, 2] = cfdPts0[j, 2]

            mesh.setSurfaceCoordinates(cfdPoints_deformed)
            mesh.warpMesh()
            m = mesh.getSolverGrid()
            CFDSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

        CFDSolver._updateGeomInfo = True
        CFDSolver.updateGeometryInfo()

    return CFDSolver, ap


def test_dot_product(CFDSolver, ap):
    """
    Dot product test: verifies forward and backward mode consistency.

    For a matrix J = dR/dw:
      <J * wDot, dwBar> == <wDot, J^T * dwBar>

    This is the gold standard test for adjoint correctness.
    """
    printroot("\n" + "=" * 60)
    printroot("DOT PRODUCT TEST: Forward vs Backward mode consistency")
    printroot("=" * 60)

    # Get random perturbation vectors
    wDot = CFDSolver.getStatePerturbation(314)
    dwBar = CFDSolver.getStatePerturbation(101)

    # Forward mode: resDot = J * wDot
    resDot, funcsDot = CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=False
    )

    # Backward mode: wBar = J^T * dwBar
    wBar = CFDSolver.computeJacobianVectorProductBwd(
        resBar=dwBar, wDeriv=True, xVDeriv=False, xDvDerivAero=False
    )

    # Compute dot products locally then reduce
    dotLocal1 = np.sum(resDot * dwBar)
    dotLocal2 = np.sum(wDot * wBar)

    dotGlobal1 = comm.allreduce(dotLocal1, op=MPI.SUM)
    dotGlobal2 = comm.allreduce(dotLocal2, op=MPI.SUM)

    rel_error = abs(dotGlobal1 - dotGlobal2) / (abs(dotGlobal1) + 1e-30)

    printroot(f"  <J*wDot, dwBar>   = {dotGlobal1:.15e}")
    printroot(f"  <wDot, J^T*dwBar> = {dotGlobal2:.15e}")
    printroot(f"  Relative error    = {rel_error:.6e}")

    passed = rel_error < 1e-10
    printroot(f"  Result: {'PASS' if passed else 'FAIL'}")

    return passed


def test_fwd_vs_fd(CFDSolver, ap):
    """
    Forward mode AD vs Finite Difference comparison.

    Verifies dR/dw * wDot computed by AD matches the FD approximation:
      (R(w + h*wDot) - R(w)) / h
    """
    printroot("\n" + "=" * 60)
    printroot("FORWARD MODE AD vs FINITE DIFFERENCE")
    printroot("=" * 60)

    wDot = CFDSolver.getStatePerturbation(314)

    # AD forward mode
    resDot_AD, funcsDot_AD = CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=False
    )

    # Finite difference
    resDot_FD, funcsDot_FD = CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=False, mode="FD", h=1e-8
    )

    # Compare residual derivatives
    norm_AD = np.sqrt(comm.allreduce(np.sum(resDot_AD**2), op=MPI.SUM))
    norm_FD = np.sqrt(comm.allreduce(np.sum(resDot_FD**2), op=MPI.SUM))
    norm_diff = np.sqrt(comm.allreduce(np.sum((resDot_AD - resDot_FD) ** 2), op=MPI.SUM))

    rel_error_res = norm_diff / (norm_AD + 1e-30)
    printroot(f"  Residual derivative:")
    printroot(f"    ||dR_AD||  = {norm_AD:.6e}")
    printroot(f"    ||dR_FD||  = {norm_FD:.6e}")
    printroot(f"    Relative error = {rel_error_res:.6e}")
    res_passed = rel_error_res < 1e-3
    printroot(f"    Result: {'PASS' if res_passed else 'FAIL'}")

    # Compare function derivatives
    func_passed = True
    for func in funcsDot_AD:
        val_AD = funcsDot_AD[func]
        val_FD = funcsDot_FD[func]
        if abs(val_AD) > 1e-15:
            rel_err = abs(val_AD - val_FD) / abs(val_AD)
        else:
            rel_err = abs(val_AD - val_FD)
        ok = rel_err < 1e-4
        func_passed = func_passed and ok
        printroot(f"  {func}: AD={val_AD:.10e}, FD={val_FD:.10e}, rel_err={rel_err:.6e} {'PASS' if ok else 'FAIL'}")

    return res_passed and func_passed


def test_adjoint_sens(CFDSolver, ap):
    """
    Test evalFunctionsSens - the main user-facing adjoint sensitivity method.
    This solves the adjoint system and computes dF/dX for all design variables.
    """
    printroot("\n" + "=" * 60)
    printroot("ADJOINT SENSITIVITY TEST (evalFunctionsSens)")
    printroot("=" * 60)

    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)

    # Check adjoint failure
    funcsCheck = {}
    CFDSolver.checkAdjointFailure(ap, funcsCheck)
    adj_converged = not funcsCheck.get("fail", True)

    printroot(f"  Adjoint converged: {adj_converged}")

    for func_name, sens_dict in funcsSens.items():
        printroot(f"  {func_name}:")
        for dv_name, val in sens_dict.items():
            if isinstance(val, np.ndarray):
                printroot(f"    d/d({dv_name}) = array, norm={np.linalg.norm(val):.6e}")
            else:
                printroot(f"    d/d({dv_name}) = {val:.10e}")

    return adj_converged


def main():
    printroot("=" * 60)
    printroot("TIME SPECTRAL ADJOINT VERIFICATION TEST")
    printroot("=" * 60)

    if not HAS_IDWARP:
        printroot("WARNING: IDWarp not available. Skipping mesh deformation.")
        printroot("The test will run without mesh deformation (rigid grid).")

    # Setup
    printroot("\n--- Setting up solver ---")
    CFDSolver, ap = setup_ts_solver()

    # Solve forward
    printroot("\n--- Running forward solve ---")
    CFDSolver(ap)

    # Check solution
    funcs = {}
    CFDSolver.checkSolutionFailure(ap, funcs)
    if funcs.get("fail", True):
        printroot("ERROR: Forward solve failed!")
        sys.exit(1)

    # Evaluate functions
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    printroot("\nFunction values:")
    for k, v in funcs.items():
        printroot(f"  {k} = {v:.10e}")

    # Initialize residual for adjoint tests
    CFDSolver.getResidual(ap)

    # Run tests
    results = {}

    # Test 1: Dot product test (most important)
    results["dot_product"] = test_dot_product(CFDSolver, ap)

    # Test 2: Forward AD vs FD
    results["fwd_vs_fd"] = test_fwd_vs_fd(CFDSolver, ap)

    # Test 3: Adjoint sensitivity
    results["adjoint_sens"] = test_adjoint_sens(CFDSolver, ap)

    # Summary
    printroot("\n" + "=" * 60)
    printroot("SUMMARY")
    printroot("=" * 60)
    all_passed = True
    for name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        printroot(f"  {name}: {status}")
        all_passed = all_passed and passed

    printroot(f"\nOverall: {'ALL TESTS PASSED' if all_passed else 'SOME TESTS FAILED'}")
    printroot("=" * 60)

    return 0 if all_passed else 1


if __name__ == "__main__":
    ret = main()
    sys.exit(ret)
