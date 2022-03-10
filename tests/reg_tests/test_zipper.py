import unittest
import os
import copy
from adflow import ADFLOW
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_simple_cart_cube


baseDir = os.path.dirname(os.path.abspath(__file__))

# Tests for zipper meshes using a mesh with two overlapping cubes

gridFile = os.path.join(baseDir, "../../input_files/cube_overset.cgns")
commonTestOptions = {
    "gridFile": gridFile,
    "equationType": "RANS",
    "writeTecplotSurfaceSolution": True,
    "monitorVariables": ["cpu", "resrho", "resturb", "cd"],
    "volumeVariables": ["resrho", "resturb", "cp", "mach", "blank"],
    "mgcycle": "sg",
    "L2Convergence": 1e-13,
    "nCycles": 500,
    "useANKSolver": True,
    "useNKSolver": True,
    "nearWallDist": 1.0,
}


class TestCubeOverset(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):

        super().setUp()

        # Start with the default testing options dictionary
        options = copy.copy(adflowDefOpts)

        # Set the output directory
        options["outputDirectory"] = os.path.join(baseDir, options["outputDirectory"])

        # These are the modified options common to these tests
        options.update(commonTestOptions)

        # Define the AeroProblem
        self.ap = copy.deepcopy(ap_simple_cart_cube)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

    def test_convergence(self):

        # Solve the flow
        self.CFDSolver(self.ap)

        # Check that the flow converged
        funcs = {}
        self.CFDSolver.checkSolutionFailure(self.ap, funcs)
        self.assertFalse(funcs["fail"])


if __name__ == "__main__":
    unittest.main()
