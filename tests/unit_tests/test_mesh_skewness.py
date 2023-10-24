import unittest
import os
from adflow import ADFLOW
from idwarp import USMesh
import numpy as np
import sys

baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir, "../reg_tests"))

from reg_aeroproblems import ap_tutorial_wing  # noqa E402


class BasicTests(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        gridFile = "input_files/cube_1x1x1.cgns"
        self.options = {
            "gridfile": os.path.join(baseDir, "../../", gridFile),
            "useSkewnessCheck": True,
            "meshMaxSkewness": 0.5,
            "printBadlySkewedCells": True,
        }
        self.idwarp_options = {
            "gridfile": os.path.join(baseDir, "../../", gridFile),
        }
        eps = 1e-12
        self.warp_no_fail = 1 - eps
        self.warp_fail = 1 + eps

        # only here for debugging
        self.write_mesh = False

    def test_skew_x(self):
        hold, move = 1, 0

        # warp the mesh slightly
        CFDSolver = self.create_solver()
        self.warp_mesh(CFDSolver, hold, move, self.warp_no_fail)
        self.assertFalse(CFDSolver.adflow.killsignals.fatalfail)

        # warp it more
        CFDSolver = self.create_solver()
        self.warp_mesh(CFDSolver, hold, move, self.warp_fail)
        self.assertTrue(CFDSolver.adflow.killsignals.fatalfail)

    def test_skew_y(self):
        hold, move = 0, 1

        # warp the mesh slightly
        CFDSolver = self.create_solver()
        self.warp_mesh(CFDSolver, hold, move, self.warp_no_fail)
        self.assertFalse(CFDSolver.adflow.killsignals.fatalfail)

        # warp it more
        CFDSolver = self.create_solver()
        self.warp_mesh(CFDSolver, hold, move, self.warp_fail)
        self.assertTrue(CFDSolver.adflow.killsignals.fatalfail)

    def test_skew_z(self):
        hold, move = 0, 2

        # warp the mesh slightly
        CFDSolver = self.create_solver()
        self.warp_mesh(CFDSolver, hold, move, self.warp_no_fail)
        self.assertFalse(CFDSolver.adflow.killsignals.fatalfail)

        # warp it more
        CFDSolver = self.create_solver()
        self.warp_mesh(CFDSolver, hold, move, self.warp_fail)
        self.assertTrue(CFDSolver.adflow.killsignals.fatalfail)

    def create_solver(self):
        CFDSolver = ADFLOW(options=self.options, debug=False)
        mesh = USMesh(options=self.idwarp_options)
        CFDSolver.setMesh(mesh)
        return CFDSolver

    def warp_mesh(self, CFDSolver, hold, move, amount):
        # setup warping
        coords = CFDSolver.mesh.getSurfaceCoordinates().copy()
        index = np.argwhere(coords[:, hold] == 0)
        coords[index, move] += amount

        # warp
        CFDSolver.mesh.setSurfaceCoordinates(coords)
        CFDSolver.setAeroProblem(ap_tutorial_wing)  # dummy AP
        if self.write_mesh:
            CFDSolver.writeMeshFile(f"skewed_{hold}_{move}_{amount}.cgns")


if __name__ == "__main__":
    unittest.main()
