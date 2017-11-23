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

sys.path.append(os.path.abspath('../../'))
from python.pyADflow import ADFLOW
from pygeo import DVGeometry
from pywarpustruct import USMesh
from pyspline import Curve

sys.path.append(os.path.abspath('../'))
from python.om_adflow import OM_ADFLOW
from python.om_utils import get_dvs_and_cons

from om_commonUtils import assert_funcs_equal, assert_funcsSens_equal


class StandardCompareTests(unittest.TestCase):

    def run_compare(self, aeroOptions, meshOptions, gridFile, ffdFile, ap): 
        solve = True
        if 'solve' not in sys.argv:
            aeroOptions['restartfile'] = gridFile
            solve = False

        for dv in defaultAeroDVs:
            ap.addDV(dv)
        ap.addDV('alpha')
        ap.addDV('mach')
        
        # Create the solver
        CFDSolver = ADFLOW(options=aeroOptions, debug=False)

        if solve:
            # We are told that we must first solve the problem, most likely
            # for a training run. 
            CFDSolver(ap)

        # # Setup geometry/mesh
        DVGeo = DVGeometry(ffdFile)
        nTwist = 6
        DVGeo.addRefAxis('wing', Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                       y=numpy.zeros(nTwist),
                                       z=numpy.linspace(0,14, nTwist), k=2))
        def twist(val, geo):
            for i in range(nTwist):
                geo.rot_z['wing'].coef[i] = val[i]

        DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
        DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
        mesh = USMesh(options=meshOptions)
        CFDSolver.setMesh(mesh)
        CFDSolver.setDVGeo(DVGeo)

        res = CFDSolver.getResidual(ap)
        
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        #Some kind of fortran memory cross-contamination 
        del(CFDSolver)
        del(DVGeo)
        del(mesh)
        gc.collect()

        ##########################################
        # Run things through OpenMDAO
        ##########################################
        prob = Problem()
        DVGeo = DVGeometry(ffdFile)
        nTwist = 6
        DVGeo.addRefAxis('wing', Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                       y=numpy.zeros(nTwist),
                                       z=numpy.linspace(0,14, nTwist), k=2))
        def twist(val, geo):
            for i in range(nTwist):
                geo.rot_z['wing'].coef[i] = val[i]

        DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
        DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
        mesh = USMesh(options=meshOptions)
        prob.model = OM_ADFLOW(ap=ap, aero_options=aeroOptions, mesh_options=meshOptions, 
                               dvgeo=DVGeo, 
                               debug=True, owns_indeps=True)

        prob.setup()
        # from openmdao.api import view_model 
        # view_model(prob)
        # exit()
        prob.model.states._do_solve = solve
        prob.model.functionals._do_solve = solve    
        prob.final_setup()
        prob.run_model()
        # prob.model.run_linearize()

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
        from tests.test13 import options as aeroOptions, meshOptions, gridFile, ffdFile, ap

        self.run_compare(aeroOptions, meshOptions, gridFile, ffdFile, ap)

    def test14(self): 
        from tests.test14 import options as aeroOptions, meshOptions, gridFile, ffdFile, ap

        self.run_compare(aeroOptions, meshOptions, gridFile, ffdFile, ap)



    

if __name__ == "__main__": 
    unittest.main()




