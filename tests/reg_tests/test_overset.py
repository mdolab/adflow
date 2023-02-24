import unittest
import os
import copy
from parameterized import parameterized_class
from adflow import ADFLOW
from reg_default_options import adflowDefOpts, IDWarpDefOpts
from reg_aeroproblems import ap_simple_cart_cube
from idwarp import USMesh
from pygeo import DVGeometry

baseDir = os.path.dirname(os.path.abspath(__file__))

# Tests for overset and zipper meshes using a mesh with two overlapping cubes


# This does not affect the hole cutting for this case but tests that cutCallback is working
def cutCallback(xCen, CGNSZoneNameIDs, cellIDs, flags):
    # Blank cells in the negative y-axis
    flags[xCen[:, 1] < 0] = 1


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
    "cutCallback": cutCallback,
}

test_params = [
    {
        "name": "frozen",
        "options": {
            "oversetUpdateMode": "frozen",
        },
    },
    {
        "name": "fast",
        "options": {
            "oversetUpdateMode": "fast",
        },
    },
    {
        "name": "full",
        "options": {
            "oversetUpdateMode": "full",
        },
    },
]


@parameterized_class(test_params)
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

        # Add the parameterized options
        options.update(self.options)

        # Define the AeroProblem
        self.ap = copy.deepcopy(ap_simple_cart_cube)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

        # Create mesh warping object
        meshOptions = copy.copy(IDWarpDefOpts)
        meshOptions.update({"gridFile": options["gridfile"]})
        self.CFDSolver.setMesh(USMesh(options=meshOptions))

        # Create geometry object
        ffdFile = os.path.join(baseDir, "../../input_files/cube_overset_ffd.xyz")
        DVGeo = DVGeometry(ffdFile)
        DVGeo.addLocalDV("shape", lower=-0.5, upper=0.5, axis="z")
        self.CFDSolver.setDVGeo(DVGeo)

        # Apply shape change
        dvDict = self.CFDSolver.DVGeo.getValues()
        dvDict["shape"][4] = 0.1
        self.CFDSolver.setAeroProblem(self.ap)
        self.CFDSolver.DVGeo.setDesignVars(dvDict)

    def test_convergence(self):
        # Solve the flow
        self.CFDSolver(self.ap)

        # Check that the flow converged
        funcs = {}
        self.CFDSolver.checkSolutionFailure(self.ap, funcs)
        self.assertFalse(funcs["fail"])


if __name__ == "__main__":
    unittest.main()
