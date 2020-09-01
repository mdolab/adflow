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
from pygeo import DVGeometry
from idwarp import USMesh
from pyspline import Curve

sys.path.append(os.path.abspath('../../'))

sys.path.append(os.path.abspath('../'))
from adflow.om_utils import get_dvs_and_cons

from om_commonUtils import assert_funcs_equal, assert_funcsSens_equal


class Tests(unittest.TestCase):

    def run_compare(self, setup_cb, ap, gridFile, ffdFile): 
       
        for dv in defaultAeroDVs:
            ap.addDV(dv)
        ap.addDV('alpha')
        ap.addDV('mach')
        
        # Create the solver
        CFDSolver, mesh, DVGeo, _ = setup_cb(MPI.COMM_WORLD)

        solve = True
        if 'solve' not in sys.argv:
            CFDSolver.setOption('restartfile', gridFile)
            solve = False

        if solve:
            # We are told that we must first solve the problem, most likely
            # for a training run. 
            CFDSolver(ap)

        res = CFDSolver.getResidual(ap)
        
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        ##########################################
        # Run things through OpenMDAO
        ##########################################
        prob = Problem()
      
        prob.model = OM_ADFLOW(ap=ap, setup_cb=setup_cb, 
                               debug=True, owns_indeps=True)

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

        dvs, _ = get_dvs_and_cons(ap=ap, geo=DVGeo) 
        dv_names = [args[0] for args,kwargs in dvs]
        om_funcsSens = prob.compute_totals(of=['functionals.{}'.format(s) for s in ap.evalFuncs], 
                                          wrt=dv_names, return_format='dict')

        # om_funcsSens = prob.compute_totals(of=['functionals.fx'], 
        #                                   wrt=dv_names, return_format='dict')

        # print(om_funcsSens)
        # print()
        # print(funcsSens)

        assert_rel_error(self, res, prob.model._residuals['states.states'], tolerance=1e-5)    

        assert_funcsSens_equal(self, ap, funcsSens, om_funcsSens, tolerance=1e-9)         

    def test13(self): 
        from tests.test13 import setup_cb, gridFile, ffdFile, ap

        self.run_compare(setup_cb, ap, gridFile, ffdFile)

    def test14(self): 
        from tests.test14 import setup_cb, gridFile, ffdFile, ap

        self.run_compare(setup_cb, ap, gridFile, ffdFile)



    

if __name__ == "__main__": 
    unittest.main()




