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
from adflow import ADFLOW, OM_ADFLOW

sys.path.append(os.path.abspath('../../'))

sys.path.append(os.path.abspath('../'))

from om_commonUtils import assert_funcs_equal


class Tests(unittest.TestCase):

    def run_compare(self, aero_options, gridFile, ap): 

        solve = True
        if 'solve' not in sys.argv:
            aero_options['restartfile'] = gridFile
            solve = False


        CFDSolver = ADFLOW(options=aero_options, debug=False)

        ##########################################
        # Run the strait ADflow code
        ##########################################
        # # this call is needed to initialize the state and resid vectors
        CFDSolver.getStateSize()

        for dv in defaultAeroDVs:
            ap.addDV(dv)

        if solve:
            # We are told that we must first solve the problem, most likely
            # for a training run. 
            CFDSolver(ap)

        # Check the residual
        res = CFDSolver.getResidual(ap)
        # res /= totalR0
        
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, defaultFuncList)

        # # Get and check the states
        states = CFDSolver.getStates()

        #force it to restart the flow
        if CFDSolver.getOption('restartFile') is not None: 
            CFDSolver.adflow.inputiteration.mgstartlevel = 1
            CFDSolver.adflow.initializeflow.initflowrestart()


        ##########################################
        # Run things through OpenMDAO
        ##########################################
        prob = Problem()
        prob.model = OM_ADFLOW(ap=ap, aero_options=aero_options, 
                               debug=True, owns_indeps=True, solver=CFDSolver)

        prob.setup()
        prob.model.states._do_solve = solve
        prob.model.functionals._do_solve = solve    
        prob.final_setup()
        prob.run_model()
        
        assert_rel_error(self, res, prob.model._residuals['states.states'], tolerance=1e-7)    
        assert_rel_error(self, states, prob['states.states'], tolerance=1e-10)
        assert_funcs_equal(self, ap, funcs, prob, tolerance=1e-10)         

    def test1(self): 
        from tests.test1 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test2(self): 
        from tests.test2 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test3(self): 
        from tests.test3 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test4(self): 
        from tests.test4 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test5(self): 
        from tests.test5 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test6(self): 
        from tests.test6 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test7(self): 
        from tests.test7 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test8(self): 
        from tests.test8 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)


if __name__ == "__main__": 
    unittest.main()




