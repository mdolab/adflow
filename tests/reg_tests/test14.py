from __future__ import print_function
import unittest
import numpy
from baseclasses import BaseRegTest

class RegTest14(unittest.TestCase):
    '''
    Test 14: MDO tutorial -- Rans -- Adjoint Test
    '''
    N_PROCS = 4

    def setUp(self):
        self.ref_file = 'reg_tests/ref/test14.ref'

    def train(self):
        with BaseRegTest(self.ref_file, train=True) as handler:
            self.regression_test(handler)

    def test(self):
        with BaseRegTest(self.ref_file, train=False) as handler:
            self.regression_test(handler)

    def regression_test(self, handler):
        '''
        This is where the actual testing happens.
        '''
        import copy
        from mpi4py import MPI
        from baseclasses import AeroProblem
        from reg_tests.commonUtils import adjoint_test, adflowDefOpts
        from adflow import ADFLOW
        from pyspline import Curve
        from pygeo import DVGeometry
        from idwarp import USMesh
        comm = MPI.COMM_WORLD

        gridFile = 'input_files/mdo_tutorial_rans_scalar_jst.cgns'
        ffdFile = 'input_files/mdo_tutorial_ffd.fmt'

        options = copy.copy(adflowDefOpts)
        options.update({
            'gridfile': gridFile,
            'restartfile': gridFile,
            'mgcycle':'2w',
            'equationtype':'RANS',
            'smoother':'dadi',
            'cfl':1.5,
            'cflcoarse':1.25,
            'resaveraging':'noresaveraging',
            'nsubiter':3,
            'nsubiterturb':3,
            'ncyclescoarse':100,
            'ncycles':1000,
            'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cmz','yplus','totalr'],
            'usenksolver':True,
            'l2convergence':1e-14,
            'l2convergencecoarse':1e-4,
            'nkswitchtol':1e-3,
            'adjointl2convergence': 1e-14,
            'frozenturbulence':False,
            'blockSplitting':False,
        })

        # Setup aeroproblem, cfdsolver, mesh and geometry.
        ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, P=20000.0,
                    T=220.0, areaRef=45.5, chordRef=3.25, beta=0.0, R=287.87,
                    xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=['fx', 'mz'])

        # Create the solver
        CFDSolver = ADFLOW(options=options, comm=comm, debug=False)

        # Create mesh warper
        meshOptions={'gridFile':gridFile}
        mesh = USMesh(options=meshOptions, comm=comm)
        CFDSolver.setMesh(mesh)

        # Setup geometry/mesh
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
        CFDSolver.setDVGeo(DVGeo)

        # Run adjoint test
        adjoint_test(handler, CFDSolver, ap)
