import numpy
import pprint
import sys
import os
from collections import deque
import json
from mpi4py import MPI

import unittest
# =============================================================================
#                         Assert Statements 
# =============================================================================


def assert_adjoint_sens_allclose(handler, CFDSolver, ap, evalFuncs=None, rtol=1e-10, atol=1e-10):
    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens, evalFuncs=None)
    handler.root_print('Eval Functions Sens:')
    handler.root_add_dict(funcsSens, 'Eval Functions Sens:', rtol=rtol, atol=atol)

def assert_problem_size_equal(handler, CFDSolver, rtol=1e-10, atol=1e-10):
    # Now a few simple checks
    handler.root_print('Total number of state DOF')
    handler.par_add_sum(CFDSolver.getStateSize(), 'Total number of state DOF')

    handler.root_print('Total number of adjoint state DOF')
    handler.par_add_sum(CFDSolver.getAdjointStateSize(), 'Total number of adjoint state DOF')

    handler.root_print('Total number of spatial DOF')
    handler.par_add_sum(CFDSolver.getSpatialSize(), 'Total number of spatial DOF')

def assert_functions_allclose(handler, CFDSolver, ap, rtol=1e-9, atol=1e-9):
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    handler.root_print('Eval Functions:')
    handler.root_add_dict(funcs, 'Eval Functions:', rtol=rtol, atol=atol)

def assert_residuals_allclose(handler, CFDSolver, ap, rtol=1e-10, atol=1e-10):
    # Check the residual
    res = CFDSolver.getResidual(ap)
    totalR0 = CFDSolver.getFreeStreamResidual(ap)
    res /= totalR0
    handler.root_print('Norm of residual')
    handler.par_add_norm(res, 'Norm of residual', rtol=rtol, atol=atol)

# def assert_residuals_lessthan(handler, CFDSolver, ap, rtol=1e-10, atol=1e-10):
#     # Check the residual
#     res = CFDSolver.getResidual(ap)
#     totalR0 = CFDSolver.getFreeStreamResidual(ap)
#     res /= totalR0
#     unittest.TestCase.assertLessEqual(totalR0, atol)
#     unittest.TestCase.assertLessEqual(res, rtol)
#     # handler.root_print('Norm of residual')
#     # handler.par_add_norm(res, 'Norm of residual', rtol=rtol, atol=atol)

def assert_forces_allclose(handler, CFDSolver, rtol=1e-10, atol=1e-10):
    CFDSolver.setOption('forcesAsTractions', False)

    forces = CFDSolver.getForces()
    handler.root_print('Sum of Forces x')
    handler.par_add_sum(forces[:, 0], 'Sum of Forces x', rtol=rtol, atol=atol)
    handler.root_print('Sum of Forces y')
    handler.par_add_sum(forces[:, 1], 'Sum of Forces y', rtol=rtol, atol=atol)
    handler.root_print('Sum of Forces z')
    handler.par_add_sum(forces[:, 2], 'Sum of Forces z', rtol=rtol, atol=atol)

def assert_tractions_allclose(handler, CFDSolver, rtol=1e-10, atol=1e-10):
        # Reset the option
    CFDSolver.setOption('forcesAsTractions', True)


    # Check the tractions/forces
    forces = CFDSolver.getForces()
    handler.root_print('Sum of Tractions x')
    handler.par_add_sum(forces[:, 0], 'Sum of Tractions x', rtol=rtol, atol=atol)
    handler.root_print('Sum of Tractions y')
    handler.par_add_sum(forces[:, 1], 'Sum of Tractions y', rtol=rtol, atol=atol)
    handler.root_print('Sum of Tractions z')
    handler.par_add_sum(forces[:, 2], 'Sum of Tractions z', rtol=rtol, atol=atol)


    # Reset the option
    CFDSolver.setOption('forcesAsTractions', False)



def assert_states_allclose(handler, CFDSolver, rtol=1e-10, atol=1e-10):
    # Get and check the states
    handler.root_print('Norm of state vector')
    states = CFDSolver.getStates()
    handler.par_add_norm(states, 'Norm of state vector', rtol=rtol, atol=atol)

def assert_fwd_mode_allclose(handler, CFDSolver, ap, rtol=1e-10, atol=1e-10, seed=314):
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
    handler.root_print('# ---------------------------------------------------#')
    handler.root_print('#             Forward mode testing                   #')
    handler.root_print('# ---------------------------------------------------#')

    handler.root_print('-> Derivatives with respect to states. wDot, ')
    wDot = CFDSolver.getStatePerturbation(314)

    resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

    handler.root_print('||dR/dw * wDot||')
    handler.par_add_norm(resDot, '||dR/dw * wDot||', rtol=rtol, atol=atol)

    handler.root_print('dFuncs/dw * wDot')
    handler.root_add_dict(funcsDot, 'dFuncs/dw * wDot', rtol=rtol, atol=atol)

    handler.root_print('||dF/dw * wDot||')
    handler.par_add_norm(fDot, '||dF/dw * wDot||', rtol=rtol, atol=atol)

    handler.root_print('-> Derivatives with respect to nodes')
    xVDot = CFDSolver.getSpatialPerturbation(314)

    resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

    handler.root_print('||dR/dXv * xVDot||')
    handler.par_add_norm(resDot, '||dR/dXv * xVDot||', rtol=rtol, atol=atol)

    # These can be finiky sometimes so a bigger tolerance.
    handler.root_print('dFuncs/dXv * xVDot')
    handler.root_add_dict(funcsDot, 'dFuncs/dXv * xVDot', rtol=rtol*10, atol=atol*10)

    handler.root_print('||dF/dXv * xVDot||')
    handler.par_add_norm(fDot, '||dF/dXv * xVDot||', rtol=rtol, atol=atol)

    handler.root_print('-> Derivatives with respect to extra variables')

    
    for aeroDV in ap.DVs.values():
        key = aeroDV.key
        handler.root_print('  -> %s'%key)
        xDvDot = {key:1.0}

        resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

        handler.root_print('||dR/d%s||'%key)
        handler.par_add_norm(resDot, '||dR/d%s||'%key, rtol=rtol, atol=atol)

        handler.root_print('dFuncs/d%s'%key)
        handler.root_add_dict(funcsDot, 'dFuncs/d%s'%key, rtol=rtol, atol=atol)

        handler.root_print('||dF/d%s||'%key)
        handler.par_add_norm(fDot, '||dF/d%s||'%key, rtol=rtol, atol=atol)


def assert_bwd_mode_allclose(handler, CFDSolver, ap, rtol=1e-10, atol=1e-10, seed=314):
    handler.root_print('# ---------------------------------------------------#')
    handler.root_print('#             Reverse mode testing                   #')
    handler.root_print('# ---------------------------------------------------#')

    handler.root_print('-> Res bar Seed')
    dwBar = CFDSolver.getStatePerturbation(314)

    wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
        resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

    handler.root_print('||dwBar^T * dR/dw||')
    handler.par_add_norm(wBar, '||dwBar^T * dR/dw||', rtol=rtol, atol=atol)

    handler.root_print('||dwBar^T * dR/dXv||')
    norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
    handler.root_add_val(norm, '||dwBar^T * dR/dXv||', rtol=rtol, atol=atol)

    handler.root_print('||dwBar^T * dR/xDv||')
    handler.root_add_dict(xDvBar, '||dwBar^T * dR/xDv||', rtol=rtol, atol=atol)

    handler.root_print('-> F Bar Seed')
    fBar = CFDSolver.getSurfacePerturbation(314)

    wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
        fBar=fBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

    handler.root_print('||FBar^T * dF/dw||')
    handler.par_add_norm(wBar, '||FBar^T * dF/dw||', rtol=rtol, atol=atol)

    handler.root_print('||FBar^T * dF/dXv||')
    norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
    handler.root_add_val(norm, '||FBar^T * dF/dXv||', rtol=rtol, atol=atol)

    handler.root_print('||FBar^T * dF/xDv||')
    handler.root_add_dict(xDvBar, '||FBar^T * dF/xDv||', rtol=rtol, atol=atol)

    handler.root_print(' -> Objective Seeds')

    for key in ap.evalFuncs:
        handler.root_print('  -> %s'%key)
        funcsBar = {key:1.0}

        wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
            funcsBar=funcsBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

        handler.root_print('||d%s/dw||'%key)
        handler.par_add_norm(wBar, '||d%s/dw||'%key, rtol=rtol, atol=atol)

        handler.root_print('||d%s/dXv||'%key)
        norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
        handler.root_add_val(norm, '||d%s/dXv||'%key, rtol=rtol, atol=atol)

        handler.root_print('||d%s/dXdv||'%key)
        handler.root_add_dict(xDvBar,'||d%s/dXdv||'%key, rtol=rtol, atol=atol)


def assert_dot_products_allclose(handler, CFDSolver, rtol=1e-10, atol=1e-10, seed=314):
    handler.root_print('# ---------------------------------------------------#')
    handler.root_print('#                 Dot product Tests                  #')
    handler.root_print('# ---------------------------------------------------#')

    # Create a common set of seeds
    wDot = CFDSolver.getStatePerturbation(314)
    xVDot= CFDSolver.getSpatialPerturbation(314)
    dwBar = CFDSolver.getStatePerturbation(314)
    fBar = CFDSolver.getSurfacePerturbation(314)

    handler.root_print ('Dot product test for w -> R')

    dwDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
    wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar,  wDeriv=True)

    dotLocal1 = numpy.sum(dwDot*dwBar)
    dotLocal2= numpy.sum(wDot*wBar)

    handler.par_add_sum(dotLocal2, 'Dot product test for w -> R', rtol=rtol, atol=atol)
    handler.par_add_sum(dotLocal1, 'Dot product test for w -> R', rtol=rtol, atol=atol)

    handler.root_print ('Dot product test for Xv -> R')
    dwDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, residualDeriv=True)
    xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar,  xVDeriv=True)

    dotLocal1 = numpy.sum(dwDot*dwBar)
    dotLocal2 = numpy.sum(xVDot*xVBar)

    # For laminar/Rans the DP test for nodes is generally appears a
    # little less accurate. This is ok. It has to do with the very
    # small offwall spacing on laminar/RANS meshes.
    handler.par_add_sum(dotLocal1, 'Dot product test for Xv -> R', rtol=rtol*10, atol=atol*10)
    handler.par_add_sum(dotLocal2, 'Dot product test for Xv -> R', rtol=rtol*10, atol=atol*10)

    handler.root_print ('Dot product test for w -> F')
    wBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar,  wDeriv=True)
    fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, fDeriv=True)

    dotLocal1 = numpy.sum(fDot.flatten()*fBar.flatten())
    dotLocal2 = numpy.sum(wDot*wBar)

    handler.par_add_sum(dotLocal1, 'Dot product test for w -> F', rtol=rtol, atol=atol)
    handler.par_add_sum(dotLocal2, 'Dot product test for w -> F', rtol=rtol, atol=atol)

    handler.root_print ('Dot product test for xV -> F')

    fDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, fDeriv=True)
    xVBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar,  xVDeriv=True)

    dotLocal1 = numpy.sum(fDot.flatten()*fBar.flatten())
    dotLocal2 = numpy.sum(xVDot*xVBar)

    handler.par_add_sum(dotLocal1, 'Dot product test for xV -> F', rtol=rtol, atol=atol)
    handler.par_add_sum(dotLocal2, 'Dot product test for xV -> F', rtol=rtol, atol=atol)


    handler.root_print ('Dot product test for (w, xV) -> (dw, F)')

    dwDot, fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, xVDot=xVDot,
                                                            residualDeriv=True, fDeriv=True)
    wBar, xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, fBar=fBar,
                                                            wDeriv=True, xVDeriv=True)

    dotLocal1 = numpy.sum(dwDot*dwBar) + numpy.sum(fDot.flatten()*fBar.flatten())
    dotLocal2= numpy.sum(wDot*wBar) + numpy.sum(xVDot*xVBar)

    handler.par_add_sum(dotLocal1, 'Dot product test for (w, xV) -> (dw, F)', rtol=rtol, atol=atol)
    handler.par_add_sum(dotLocal2, 'Dot product test for (w, xV) -> (dw, F)', rtol=rtol, atol=atol)
