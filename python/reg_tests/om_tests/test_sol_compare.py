############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW with OpenMDAO
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
import gc

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

    def run_compare(self, setup_cb, ap, gridFile, dvs=None, defaultFuncs=True, ): 

        if dvs is None: 
            for dv in defaultAeroDVs:
                ap.addDV(dv)
        else: 
            for dv in dvs: 
                ap.addDV(dv)

        evalFuncs = defaultFuncList
        if not defaultFuncs: 
            evalFuncs = ap.evalFuncs

        ##########################################
        # Run the strait ADflow code
        ##########################################
        CFDSolver, _, _, _, = setup_cb(MPI.COMM_WORLD)

        solve = True
        if 'solve' not in sys.argv:
            CFDSolver.setOption('restartfile', gridFile)
            solve = False

        CFDSolver(ap)

        # Check the residual
        res = CFDSolver.getResidual(ap)
        
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, evalFuncs)

        # # Get and check the states
        states = CFDSolver.getStates()

        ##########################################
        # Run things through OpenMDAO
        ##########################################

        prob = Problem()
        prob.model = OM_ADFLOW(ap=ap, setup_cb=setup_cb, 
                               debug=False, owns_indeps=True)

        prob.setup(mode='rev')
        if not solve:
            prob.model.solver.resetFlow(ap)
            prob.model.solver.setOption('restartfile', gridFile)
            prob.model.solver.adflow.inputiteration.mgstartlevel = 1
            prob.model.solver.adflow.initializeflow.initflowrestart()
            prob.model.solver.getResidual(ap) # this does some kind of memory allocation that is needed

        prob.run_model()


        assert_rel_error(self, res, prob.model._residuals['states.states'], tolerance=5e-7)    
        assert_rel_error(self, states, prob['states.states'], tolerance=3e-10)

        assert_funcs_equal(self, ap, funcs, prob, tolerance=3e-10)  

    def test9(self): 
        from tests.test9 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)


    def test11(self): 
        from tests.test11 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)


    def test12(self): 
        from tests.test12 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile)

    # def test15(self):

    #     self.skipTest('time accurate stuff not working in OM yet')
    #     # This test case runs, but OpeNMDAO runs to max iterations without converging
    #     # OpenMDAO seems to start out at a different residual as well. 
    #     # there is probably some restart stuff in the AP the needs to be cleared out

    #     from tests.test15 import options as aero_options, gridFile, ap
    #     self.run_compare(aero_options, gridFile, ap, 
    #                      dvs=['alpha', 'chordRef', 'xRef', 'reynolds', 'xRot', 'areaRef', 'T', 'reynoldsLength', 'mach'])

    #     os.system('rm  0012pitching*')

    # def test16(self): 
    #     from tests.test16 import setup_cb, gridFile, ap
        
    #     self.run_compare(setup_cb, ap, gridFile, 
    #                      dvs=['areaRef', 'T', 'chordRef', 'P', 'alpha', 'mach'], 
    #                      defaultFuncs=False)
    
    def test17(self): 
        from tests.test17 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile,  
                         dvs=['areaRef', 'chordRef', 'alpha', 'altitude', 'mach'], 
                         defaultFuncs=False) 

    def test18(self): 
        from tests.test18 import setup_cb, gridFile, ap

        self.run_compare(setup_cb, ap, gridFile, 
                         dvs=['areaRef', 'chordRef', 'alpha', 'altitude', 'mach'], 
                         defaultFuncs=False) 

if __name__ == "__main__": 
    unittest.main()




