from __future__ import print_function
############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW with OpenMDAO
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy

import unittest

from mpi4py import MPI
from baseclasses import AeroProblem
from mdo_regression_helper import *
from commonUtils import *

from openmdao.api import Problem, IndepVarComp
from openmdao.devtools.testutil import assert_rel_error

sys.path.append(os.path.abspath('../../'))
from python.pyADflow import ADFLOW

sys.path.append(os.path.abspath('../'))
from python.om_adflow import OM_ADFLOW

from om_commonUtils import assert_funcs_equal


class SolutionCompareTests(unittest.TestCase):

    def run_compare(self, aero_options, gridFile, ap): 

        for dv in defaultAeroDVs:
            ap.addDV(dv)


                          

        ##########################################
        # Run the strait ADflow code
        ##########################################
        # # this call is needed to initialize the state and resid vectors
        CFDSolver = ADFLOW(options=aero_options, debug=False)

        CFDSolver(ap)

        # Check the residual
        res = CFDSolver.getResidual(ap)
        
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, defaultFuncList)

        # # Get and check the states
        states = CFDSolver.getStates()

        #Some kind of fortran memory cross-contamination 
        del(CFDSolver)


        ##########################################
        # Run things through OpenMDAO
        ##########################################
        prob = Problem()
        prob.model = OM_ADFLOW(ap=ap, aero_options=aero_options, debug=False, owns_indeps=True)

        prob.setup()
        prob.run_model()


        assert_rel_error(self, res, prob.model._residuals['states.states'], tolerance=1e-8)    
        assert_rel_error(self, prob['states.states'], states, tolerance=1e-10)
        assert_funcs_equal(self, ap, prob, funcs, tolerance=1e-10)  

    def test9(self): 
        from tests.test9 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)


    def test11(self): 
        from tests.test11 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)


    def test12(self): 
        from tests.test12 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)



if __name__ == "__main__": 
    unittest.main()




