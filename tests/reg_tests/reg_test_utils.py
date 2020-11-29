import numpy
from baseclasses.BaseRegTest import getTol

# =============================================================================
#                         Assert Statements
# =============================================================================


def assert_adjoint_sens_allclose(handler, CFDSolver, ap, evalFuncs=None, **kwargs):
    rtol, atol = getTol(**kwargs)
    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens, evalFuncs=None)
    handler.root_print("Eval Functions Sens:")
    handler.root_add_dict("Eval Functions Sens:", funcsSens, rtol=rtol, atol=atol)


def assert_problem_size_equal(handler, CFDSolver, **kwargs):
    rtol, atol = getTol(**kwargs)
    # Now a few simple checks
    handler.root_print("Total number of state DOF")
    handler.par_add_sum("Total number of state DOF", CFDSolver.getStateSize())

    handler.root_print("Total number of adjoint state DOF")
    handler.par_add_sum("Total number of adjoint state DOF", CFDSolver.getAdjointStateSize())

    handler.root_print("Total number of spatial DOF")
    handler.par_add_sum("Total number of spatial DOF", CFDSolver.getSpatialSize())


def assert_functions_allclose(handler, CFDSolver, ap, **kwargs):
    rtol, atol = getTol(**kwargs)
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    handler.root_print("Eval Functions:")
    handler.root_add_dict("Eval Functions:", funcs, rtol=rtol, atol=atol)


def assert_residuals_allclose(handler, CFDSolver, ap, **kwargs):
    rtol, atol = getTol(**kwargs)
    # Check the residual
    res = CFDSolver.getResidual(ap)
    totalR0 = CFDSolver.getFreeStreamResidual(ap)
    res /= totalR0
    handler.root_print("Norm of residual")
    handler.par_add_norm("Norm of residual", res, rtol=rtol, atol=atol)


# def assert_residuals_lessthan(handler, CFDSolver, ap, **kwargs):
#     rtol, atol = getTol(**kwargs)
#     # Check the residual
#     res = CFDSolver.getResidual(ap)
#     totalR0 = CFDSolver.getFreeStreamResidual(ap)
#     res /= totalR0
#     unittest.TestCase.assertLessEqual(totalR0, atol)
#     unittest.TestCase.assertLessEqual(res, rtol)
#     # handler.root_print('Norm of residual')
#     # handler.par_add_norm('Norm of residual', res, rtol=rtol, atol=atol)


def assert_forces_allclose(handler, CFDSolver, **kwargs):
    rtol, atol = getTol(**kwargs)
    CFDSolver.setOption("forcesAsTractions", False)

    forces = CFDSolver.getForces()
    handler.root_print("Sum of Forces x")
    handler.par_add_sum("Sum of Forces x", forces[:, 0], rtol=rtol, atol=atol)
    handler.root_print("Sum of Forces y")
    handler.par_add_sum("Sum of Forces y", forces[:, 1], rtol=rtol, atol=atol)
    handler.root_print("Sum of Forces z")
    handler.par_add_sum("Sum of Forces z", forces[:, 2], rtol=rtol, atol=atol)


def assert_tractions_allclose(handler, CFDSolver, **kwargs):
    rtol, atol = getTol(**kwargs)
    # Reset the option
    CFDSolver.setOption("forcesAsTractions", True)

    # Check the tractions/forces
    forces = CFDSolver.getForces()
    handler.root_print("Sum of Tractions x")
    handler.par_add_sum("Sum of Tractions x", forces[:, 0], rtol=rtol, atol=atol)
    handler.root_print("Sum of Tractions y")
    handler.par_add_sum("Sum of Tractions y", forces[:, 1], rtol=rtol, atol=atol)
    handler.root_print("Sum of Tractions z")
    handler.par_add_sum("Sum of Tractions z", forces[:, 2], rtol=rtol, atol=atol)

    # Reset the option
    CFDSolver.setOption("forcesAsTractions", False)


def assert_states_allclose(handler, CFDSolver, **kwargs):
    rtol, atol = getTol(**kwargs)
    # Get and check the states
    handler.root_print("Norm of state vector")
    states = CFDSolver.getStates()
    handler.par_add_norm("Norm of state vector", states, rtol=rtol, atol=atol)


def assert_fwd_mode_allclose(handler, CFDSolver, ap, seed=314, **kwargs):
    rtol, atol = getTol(**kwargs)
    # Now for the most fun part. Checking the derivatives. These are
    # generally the most important things to check. However, since the
    # checking requires random seeds, it is quite tricky to ensure that we
    # are actually doing the same thing in parallel as in serial. In order
    # to account for this we have to set the random set of numbers to
    # correspond to the full CGNS mesh volume ordering and then scatter
    # back to the correct locations. We have a special routine built into
    # adflow specifically for this purpose.

    # Our "full" jacobian looks like the following:
    #                       residuals  objectives   forces
    #                     +----------+------------+--------+
    #     state Variables |          |            |        |
    #                     +----------+----------- +--------+
    #        volume Nodes |          |            |        |
    #                     +----------+----------- +--------+
    #   "extra" Variables |          |            |        |
    #                     +----------+------------+--------+
    #
    # This defines everything that goes into adflow and everything we care
    # about coming back out. We will check all the derivatives using
    # forward mode AD, reverse mode AD as well as with the dot product
    # test.
    #
    handler.root_print("# ---------------------------------------------------#")
    handler.root_print("#             Forward mode testing                   #")
    handler.root_print("# ---------------------------------------------------#")

    handler.root_print("-> Derivatives with respect to states. wDot, ")
    wDot = CFDSolver.getStatePerturbation(314)

    resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True
    )

    handler.root_print("||dR/dw * wDot||")
    handler.par_add_norm("||dR/dw * wDot||", resDot, rtol=rtol, atol=atol)

    handler.root_print("dFuncs/dw * wDot")
    handler.root_add_dict("dFuncs/dw * wDot", funcsDot, rtol=rtol, atol=atol)

    handler.root_print("||dF/dw * wDot||")
    handler.par_add_norm("||dF/dw * wDot||", fDot, rtol=rtol, atol=atol)

    handler.root_print("-> Derivatives with respect to nodes")
    xVDot = CFDSolver.getSpatialPerturbation(314)

    resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True
    )

    handler.root_print("||dR/dXv * xVDot||")
    handler.par_add_norm("||dR/dXv * xVDot||", resDot, rtol=rtol, atol=atol)

    # These can be finiky sometimes so a bigger tolerance.
    handler.root_print("dFuncs/dXv * xVDot")
    handler.root_add_dict("dFuncs/dXv * xVDot", funcsDot, rtol=rtol * 10, atol=atol * 10)

    handler.root_print("||dF/dXv * xVDot||")
    handler.par_add_norm("||dF/dXv * xVDot||", fDot, rtol=rtol, atol=atol)

    handler.root_print("-> Derivatives with respect to extra variables")

    for aeroDV in ap.DVs.values():
        key = aeroDV.key
        handler.root_print("  -> %s" % key)
        xDvDot = {key: 1.0}

        resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True
        )

        handler.root_print("||dR/d%s||" % key)
        handler.par_add_norm("||dR/d%s||" % key, resDot, rtol=rtol, atol=atol)

        handler.root_print("dFuncs/d%s" % key)
        handler.root_add_dict("dFuncs/d%s" % key, funcsDot, rtol=rtol, atol=atol)

        handler.root_print("||dF/d%s||" % key)
        handler.par_add_norm("||dF/d%s||" % key, fDot, rtol=rtol, atol=atol)


def assert_bwd_mode_allclose(handler, CFDSolver, ap, seed=314, **kwargs):
    rtol, atol = getTol(**kwargs)
    handler.root_print("# ---------------------------------------------------#")
    handler.root_print("#             Reverse mode testing                   #")
    handler.root_print("# ---------------------------------------------------#")

    handler.root_print("-> Res bar Seed")
    dwBar = CFDSolver.getStatePerturbation(314)

    wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
        resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True
    )

    handler.root_print("||dwBar^T * dR/dw||")
    handler.par_add_norm("||dwBar^T * dR/dw||", wBar, rtol=rtol, atol=atol)

    handler.root_print("||dwBar^T * dR/dXv||")
    norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
    handler.root_add_val("||dwBar^T * dR/dXv||", norm, rtol=rtol, atol=atol)

    handler.root_print("||dwBar^T * dR/xDv||")
    handler.root_add_dict("||dwBar^T * dR/xDv||", xDvBar, rtol=rtol, atol=atol)

    handler.root_print("-> F Bar Seed")
    fBar = CFDSolver.getSurfacePerturbation(314)

    wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
        fBar=fBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True
    )

    handler.root_print("||FBar^T * dF/dw||")
    handler.par_add_norm("||FBar^T * dF/dw||", wBar, rtol=rtol, atol=atol)

    handler.root_print("||FBar^T * dF/dXv||")
    norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
    handler.root_add_val("||FBar^T * dF/dXv||", norm, rtol=rtol, atol=atol)

    handler.root_print("||FBar^T * dF/xDv||")
    handler.root_add_dict("||FBar^T * dF/xDv||", xDvBar, rtol=rtol, atol=atol)

    handler.root_print(" -> Objective Seeds")

    for key in ap.evalFuncs:
        handler.root_print("  -> %s" % key)
        funcsBar = {key: 1.0}

        wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
            funcsBar=funcsBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True
        )

        handler.root_print("||d%s/dw||" % key)
        handler.par_add_norm("||d%s/dw||" % key, wBar, rtol=rtol, atol=atol)

        handler.root_print("||d%s/dXv||" % key)
        norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
        handler.root_add_val("||d%s/dXv||" % key, norm, rtol=rtol, atol=atol)

        handler.root_print("||d%s/dXdv||" % key)
        handler.root_add_dict("||d%s/dXdv||" % key, xDvBar, rtol=rtol, atol=atol)


def assert_dot_products_allclose(handler, CFDSolver, seed=314, **kwargs):
    rtol, atol = getTol(**kwargs)
    handler.root_print("# ---------------------------------------------------#")
    handler.root_print("#                 Dot product Tests                  #")
    handler.root_print("# ---------------------------------------------------#")

    # Create a common set of seeds
    wDot = CFDSolver.getStatePerturbation(314)
    xVDot = CFDSolver.getSpatialPerturbation(314)
    dwBar = CFDSolver.getStatePerturbation(314)
    fBar = CFDSolver.getSurfacePerturbation(314)

    handler.root_print("Dot product test for w -> R")

    dwDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
    wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, wDeriv=True)

    dotLocal1 = numpy.sum(dwDot * dwBar)
    dotLocal2 = numpy.sum(wDot * wBar)

    handler.par_add_sum("Dot product test for w -> R", dotLocal1, rtol=rtol, atol=atol)
    handler.par_add_sum("Dot product test for w -> R", dotLocal2, rtol=rtol, atol=atol, check=True)

    handler.root_print("Dot product test for Xv -> R")
    dwDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, residualDeriv=True)
    xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, xVDeriv=True)

    dotLocal1 = numpy.sum(dwDot * dwBar)
    dotLocal2 = numpy.sum(xVDot * xVBar)

    # For laminar/Rans the DP test for nodes is generally appears a
    # little less accurate. This is ok. It has to do with the very
    # small offwall spacing on laminar/RANS meshes.
    handler.par_add_sum("Dot product test for Xv -> R", dotLocal1, rtol=rtol * 10, atol=atol * 10)
    handler.par_add_sum("Dot product test for Xv -> R", dotLocal2, rtol=rtol * 10, atol=atol * 10, check=True)

    handler.root_print("Dot product test for w -> F")
    wBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar, wDeriv=True)
    fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, fDeriv=True)

    dotLocal1 = numpy.sum(fDot.flatten() * fBar.flatten())
    dotLocal2 = numpy.sum(wDot * wBar)

    handler.par_add_sum("Dot product test for w -> F", dotLocal1, rtol=rtol, atol=atol)
    handler.par_add_sum("Dot product test for w -> F", dotLocal2, rtol=rtol, atol=atol, check=True)

    handler.root_print("Dot product test for xV -> F")

    fDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, fDeriv=True)
    xVBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar, xVDeriv=True)

    dotLocal1 = numpy.sum(fDot.flatten() * fBar.flatten())
    dotLocal2 = numpy.sum(xVDot * xVBar)

    handler.par_add_sum("Dot product test for xV -> F", dotLocal1, rtol=rtol, atol=atol)
    handler.par_add_sum("Dot product test for xV -> F", dotLocal2, rtol=rtol, atol=atol, check=True)

    handler.root_print("Dot product test for (w, xV) -> (dw, F)")

    dwDot, fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, xVDot=xVDot, residualDeriv=True, fDeriv=True)
    wBar, xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, fBar=fBar, wDeriv=True, xVDeriv=True)

    dotLocal1 = numpy.sum(dwDot * dwBar) + numpy.sum(fDot.flatten() * fBar.flatten())
    dotLocal2 = numpy.sum(wDot * wBar) + numpy.sum(xVDot * xVBar)

    handler.par_add_sum("Dot product test for (w, xV) -> (dw, F)", dotLocal1, rtol=rtol, atol=atol)
    handler.par_add_sum("Dot product test for (w, xV) -> (dw, F)", dotLocal2, rtol=rtol, atol=atol, check=True)
