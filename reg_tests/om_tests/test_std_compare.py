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
from openmdao.utils.assert_utils import assert_rel_error
from adflow import ADFLOW, OM_ADFLOW

sys.path.append(os.path.abspath('../../'))

sys.path.append(os.path.abspath('../'))

from om_commonUtils import assert_funcs_equal


class Tests(unittest.TestCase):

    def run_compare(self, setup_cb, ap, gridFile): 

        CFDSolver, _, _, _ = setup_cb(MPI.COMM_WORLD)

        solve = True
        if 'solve' not in sys.argv:
            CFDSolver.setOption('restartfile', gridFile)
            solve = False

        for dv in defaultAeroDVs:
            ap.addDV(dv)

        ##########################################
        # Run the strait ADflow code
        ##########################################
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
        prob.model = OM_ADFLOW(ap=ap, setup_cb=setup_cb, debug=True, owns_indeps=True )

        prob.setup(mode='rev')
        if not solve:
            prob.model.solver.resetFlow(ap)
            prob.model.solver.setOption('restartfile', gridFile)
            prob.model.solver.adflow.inputiteration.mgstartlevel = 1
            prob.model.solver.adflow.initializeflow.initflowrestart()
            prob.model.solver.getResidual(ap) # this does some kind of memory allocation that is needed

        prob.model.states._do_solve = solve
        prob.model.functionals._do_solve = solve    
        prob.final_setup()
        prob.run_model()
        
        assert_rel_error(self, res, prob.model._residuals['states.states'], tolerance=1e-7)    
        assert_rel_error(self, states, prob['states.states'], tolerance=1e-10)
        assert_funcs_equal(self, ap, funcs, prob, tolerance=1e-10)         

    def test1(self): 
        from tests.test1 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    def test2(self): 
        from tests.test2 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    def test3(self): 
        from tests.test3 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    def test4(self): 
        from tests.test4 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    def test5(self): 
        from tests.test5 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    def test6(self): 
        from tests.test6 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    def test7(self): 
        from tests.test7 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    def test8(self): 
        from tests.test8 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)


if __name__ == "__main__": 
    unittest.main()




