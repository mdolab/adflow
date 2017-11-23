from __future__ import print_function
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
from openmdao.devtools.testutil import assert_rel_error
from adflow import ADFLOW

sys.path.append(os.path.abspath('../../'))
sys.path.append(os.path.abspath('../'))
from python.om_adflow import OM_ADFLOW

from om_commonUtils import assert_funcs_equal


class Tests(unittest.TestCase):

    def run_compare(self, aero_options, gridFile, ap, dvs=None, defaultFuncs=True, setup_cb=None): 

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
        # # this call is needed to initialize the state and resid vectors
        CFDSolver = ADFLOW(options=aero_options, debug=False)
        if setup_cb is not None: 
            setup_cb(CFDSolver)

        CFDSolver(ap)

        # Check the residual
        res = CFDSolver.getResidual(ap)
        
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs, evalFuncs)

        # # Get and check the states
        states = CFDSolver.getStates()

            # #Some kind of fortran memory cross-contamination 
            # del(CFDSolver)
            # gc.collect()

        ##########################################
        # Run things through OpenMDAO
        ##########################################

        #force it to restart the flow
        if CFDSolver.getOption('restartFile') is not None: 
            CFDSolver.adflow.inputiteration.mgstartlevel = 1
            CFDSolver.adflow.initializeflow.initflowrestart()

        prob = Problem()
        prob.model = OM_ADFLOW(ap=ap, aero_options=aero_options, 
                               debug=False, owns_indeps=True, #adflow_setup_cb=setup_cb, 
                               solver=CFDSolver)

        prob.setup()
        prob.run_model()


        assert_rel_error(self, res, prob.model._residuals['states.states'], tolerance=5e-7)    
        assert_rel_error(self, states, prob['states.states'], tolerance=1e-10)

        assert_funcs_equal(self, ap, funcs, prob, tolerance=1e-10)  

    def test9(self): 
        from tests.test9 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)


    def test11(self): 
        from tests.test11 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)


    def test12(self): 
        from tests.test12 import options as aero_options, gridFile, ap

        self.run_compare(aero_options, gridFile, ap)

    def test15(self):

        self.skipTest('time accurate stuff not working in OM yet')
        # This test case runs, but OpeNMDAO runs to max iterations without converging
        # OpenMDAO seems to start out at a different residual as well. 
        # there is probably some restart stuff in the AP the needs to be cleared out

        from tests.test15 import options as aero_options, gridFile, ap
        self.run_compare(aero_options, gridFile, ap, 
                         dvs=['alpha', 'chordRef', 'xRef', 'reynolds', 'xRot', 'areaRef', 'T', 'reynoldsLength', 'mach'])

        os.system('rm  0012pitching*')

    def test16(self): 

        from tests.test16 import options as aero_options, gridFile, ap, setupADFlow
        self.run_compare(aero_options, gridFile, ap, 
                         dvs=['areaRef', 'T', 'chordRef', 'P', 'alpha', 'mach'], 
                         defaultFuncs=False, 
                         setup_cb=setupADFlow)
    
    def test17(self): 

        from tests.test17 import options as aero_options, gridFile, ap, setupADFlow
        self.run_compare(aero_options, gridFile, ap, 
                         dvs=['areaRef', 'chordRef', 'alpha', 'altitude', 'mach'], 
                         defaultFuncs=False, 
                         setup_cb=setupADFlow) 

    def test18(self): 


        from tests.test18 import options as aero_options, gridFile, ap, setupADFlow
        self.run_compare(aero_options, gridFile, ap, 
                         dvs=['areaRef', 'chordRef', 'alpha', 'altitude', 'mach'], 
                         defaultFuncs=False, 
                         setup_cb=setupADFlow) 
if __name__ == "__main__": 
    unittest.main()




