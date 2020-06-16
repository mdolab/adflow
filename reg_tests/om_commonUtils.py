# This file defines testing methods that match those in commonUtils.py, 
# but that execute using the OpenMDAO wrapper instead of bare ADflow
import numpy as np 

from mpi4py import MPI
from mdo_regression_helper import *

from commonUtils import defaultFuncList, defaultAeroDVs, adflowDefOpts, \
                        pyWarpDefOpts, IDWarpDefOpts, printHeader, parPrint

from openmdao.api import Problem 
from openmdao.utils.assert_utils import assert_rel_error


def assert_funcs_equal(test, ap, funcs, prob, tolerance):
    evalFuncs = sorted(ap.evalFuncs)

    for f_name in evalFuncs:
        adflow_name = '{}_{}'.format(ap.name, f_name)
        om_name = 'functionals.{}'.format(f_name)
        # print('func compare', f_name, funcs[adflow_name])
        assert_rel_error(test, prob[om_name], funcs[adflow_name], tolerance=tolerance)


def assert_funcsSens_equal(test, ap, funcsSens, om_funcsSens, tolerance):
     
    for f_name in ap.evalFuncs: 
        ap_f_name = '{}_{}'.format(ap.name, f_name)
        jac = funcsSens[ap_f_name]

        om_jac = om_funcsSens['functionals.{}'.format(f_name)]

        for dv_n in jac.keys(): 
            # print(f_name, dv_n,  om_jac[dv_n], jac[dv_n])
            if jac[dv_n].shape == (): # stupid trick to deal with scalars
                jac[dv_n] = np.array([[jac[dv_n]]])
            tol = tolerance
            if np.linalg.norm(jac[dv_n]) < 1e-15: 
                tol = 3
            assert_rel_error(test, om_jac[dv_n], jac[dv_n], tolerance=tol)
