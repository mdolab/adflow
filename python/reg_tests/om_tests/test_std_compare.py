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

from om_commonUtils import standardCompare




class StandardCompareTests(unittest.TestCase):

    def setUp(self): 
        solve = True
        if 'solve' not in sys.argv:
            aero_options['restartfile'] = gridFile
            solve = False
        self.solve = solve

    def test1(self): 
        from tests.test1 import options as aero_options, gridFile, ap
        
        # Create the solver
        CFDSolver = CFDSolver = ADFLOW(options=aero_options, debug=True)        
        # Create the OM wrapped solver
        OM_adflow = OM_ADFLOW(ap=ap, aero_options=aero_options, debug=True, owns_indeps=True)
        
        standardCompare(self, CFDSolver, ap, OM_adflow, self.solve)

    def test2(self): 
        from tests.test2 import options as aero_options, gridFile, ap
        
        # Create the solver
        CFDSolver = CFDSolver = ADFLOW(options=aero_options, debug=True)        
        # Create the OM wrapped solver
        OM_adflow = OM_ADFLOW(ap=ap, aero_options=aero_options, debug=True, owns_indeps=True)
        
        standardCompare(self, CFDSolver, ap, OM_adflow, self.solve)


if __name__ == "__main__": 
    unittest.main()




