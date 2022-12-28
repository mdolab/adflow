# built-ins
import unittest
import os
import copy
import numpy

# MACH testing class
from adflow import ADFLOW
from idwarp import USMesh

# import the testing utilities
import reg_test_utils as utils

from reg_default_options import adflowDefOpts

from reg_aeroproblems import ap_naca64A010_time_spectral
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))


class TestSolve(reg_test_classes.RegTest):
    """
    Tests for time spectral case for an airfoil.
    """

    N_PROCS = 2
    ref_file = "solve_euler_time_spectral_naca64A010.json"

    def setUp(self):

        # Introduce a transfer class for displacement transfer from struct to aero.
        class Transfer:

            # simplified transfer class
            # converting csd displacement to cfd surface nodes

            def __init__(self, alpha, xRot, aeroSolver):

                # takes in displacement history

                self.alpha = alpha
                self.ntimeintervalsspectral = len(alpha)

                self.aeroSolver = aeroSolver  # a shallow copy of CFD solver

                self.xRot = xRot

            def getUndeformedSurfaceNodes(self):

                self.MDGroup = self.aeroSolver.allWallsGroup

                self.cfdPts0 = self.aeroSolver.getSurfaceCoordinates(self.MDGroup, includeZipper=False)

            def setDisplacements(self):

                xRot = self.xRot
                ntimeintervalsspectral = self.ntimeintervalsspectral
                alpha = self.alpha  # notice a shallow copy introduced here; dont change the underlying obj!
                cfdPoints_init = self.cfdPts0  # notice a shallow copy introduced here; dont change the underlying obj!

                N_pts = cfdPoints_init.shape[0]

                self.cfdPts = []

                for sps in range(ntimeintervalsspectral):

                    cfdPoints_deformed = numpy.zeros((N_pts, 3))

                    ptch_loc = alpha[sps]

                    cc = numpy.cos(ptch_loc)
                    ss = numpy.sin(ptch_loc)

                    for j in range(N_pts):

                        cfdPoints_deformed[j, 0] = cc * (cfdPoints_init[j, 0] - xRot) + ss * cfdPoints_init[j, 1] + xRot
                        cfdPoints_deformed[j, 1] = -ss * (cfdPoints_init[j, 0] - xRot) + cc * cfdPoints_init[j, 1]
                        cfdPoints_deformed[j, 2] = cfdPoints_init[j, 2]

                    self.cfdPts.append(cfdPoints_deformed)

            def setVolumeMesh(self):

                ntimeintervalsspectral = self.ntimeintervalsspectral

                for sps in range(ntimeintervalsspectral):

                    self.aeroSolver.mesh.setSurfaceCoordinates(self.cfdPts[sps])
                    self.aeroSolver.mesh.warpMesh()
                    m = self.aeroSolver.mesh.getSolverGrid()
                    self.aeroSolver.adflow.warping.setgridforoneinstance(m, sps=sps + 1)

                self.aeroSolver._updateGeomInfo = True
                self.aeroSolver.updateGeometryInfo()

        # Set up the test for the parent RegTest class.
        super().setUp()

        gridFile = os.path.join(baseDir, "../../input_files/naca64A010_euler-L2.cgns")

        ntimeintervalsspectral = 3
        options = copy.copy(adflowDefOpts)
        options.update(
            {
                "gridfile": gridFile,
                "outputDirectory": os.path.join(baseDir, "../output_files"),
                "writeVolumeSolution": False,
                "writeSurfaceSolution": False,
                "blocksplitting": True,
                "useblockettes": False,
                "equationtype": "Euler",
                "equationmode": "time spectral",
                "mgcycle": "sg",
                "l2convergence": 1e-15,
                "ncycles": 200000,
                "monitorvariables": ["resrho", "cl"],
                "usenksolver": True,
                "nkswitchtol": 1e-4,
                "NKSubSpaceSize": 400,
                "applypcsubspacesize": 400,
                "useanksolver": True,
                "ankswitchtol": 1e-2,
                "anksubspacesize": 50,
                "alphafollowing": False,
                "timeintervals": ntimeintervalsspectral,
                "useexternaldynamicmesh": True,
                "usetsinterpolatedgridvelocity": True,
            }
        )

        # Grid option
        meshOptions = {
            "gridFile": gridFile,
        }

        # Setup aeroproblem
        self.ap = copy.copy(ap_naca64A010_time_spectral)

        # Motion history
        alpha_0 = 1.01
        deltaAlpha = -alpha_0 * numpy.pi / 180.0
        alpha = numpy.linspace(0.0, 2.0 * numpy.pi, ntimeintervalsspectral + 1)[:-1]
        alpha = -numpy.sin(alpha)
        alpha *= deltaAlpha

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)
        # self.CFDSolver.addSlices("z", [0.5])

        # Deform the mesh
        mesh = USMesh(options=meshOptions)
        self.CFDSolver.setMesh(mesh)

        # deformation
        xRot = 0.25  # Hard copied from the reference file.
        TSTransfer = Transfer(alpha, xRot, self.CFDSolver)
        TSTransfer.getUndeformedSurfaceNodes()
        TSTransfer.setDisplacements()
        TSTransfer.setVolumeMesh()

    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check if the solution failed
        self.assert_solution_failure()

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-8)


if __name__ == "__main__":
    unittest.main()
