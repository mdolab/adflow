import copy
import os
import sys
import unittest
from parameterized import parameterized_class
from adflow import ADFLOW


baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir, "../reg_tests"))

from reg_aeroproblems import ap_tutorial_wing  # noqa E402


@parameterized_class(
    [
        {
            "name": "three_digits",
            "writeSolutionDigits": 3,
        },
        {
            "name": "five_digits",
            "writeSolutionDigits": 5,
        },
    ]
)
class TestSolutionFileNames(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.options = {
            "outputDirectory": os.path.join(baseDir, "../output_files"),
            "writeSolutionDigits": self.writeSolutionDigits,
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "writeSurfaceSolution": False,
        }

    def test_failed_mesh(self):
        # Create the solver
        CFDSolver = ADFLOW(options=self.options)
        ap = copy.copy(ap_tutorial_wing)

        # Pretend mesh warping failed
        CFDSolver.adflow.killsignals.fatalfail = True

        # Solve
        CFDSolver(ap)

        if self.writeSolutionDigits == 3:
            self.mesh_file = os.path.join(self.options["outputDirectory"], f"failed_mesh_{ap.name}_000.cgns")
        elif self.writeSolutionDigits == 5:
            self.mesh_file = os.path.join(self.options["outputDirectory"], f"failed_mesh_{ap.name}_00000.cgns")

        # Check that a mesh file with this name exists
        self.assertTrue(os.path.isfile(self.mesh_file))

    def test_volume_solution(self):
        # Create the solver
        CFDSolver = ADFLOW(options=self.options)

        # Set the AeroProblem and call counter
        ap = copy.copy(ap_tutorial_wing)
        CFDSolver.setAeroProblem(ap)
        CFDSolver.curAP.adflowData.callCounter = 5

        CFDSolver.writeSolution()

        if self.writeSolutionDigits == 3:
            self.mesh_file = os.path.join(self.options["outputDirectory"], f"{ap.name}_005_vol.cgns")
        elif self.writeSolutionDigits == 5:
            self.mesh_file = os.path.join(self.options["outputDirectory"], f"{ap.name}_00005_vol.cgns")

        # Check that a mesh file with this name exists
        self.assertTrue(os.path.isfile(self.mesh_file))

    def tearDown(self):
        if os.path.isfile(self.mesh_file):
            os.remove(self.mesh_file)


@parameterized_class(
    [
        {
            "name": "unsteady_three_digits",
            "writeSolutionDigits": 3,
        },
    ]
)
class TestSolutionFileNamesUnsteady(unittest.TestCase):
    N_PROCS = 1

    def setUp(self):
        self.options = {
            "outputDirectory": os.path.join(baseDir, "../output_files"),
            "writeSolutionDigits": self.writeSolutionDigits,
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "writeSurfaceSolution": True,
            "writeTecplotSurfaceSolution": True,
            "writeVolumeSolution": True,
            "nCycles": 1,
            "nTimeStepsFine": 2,
            "equationMode": "unsteady",
        }

    def _cleanup_output_files(self, output_dir):
        if not os.path.isdir(output_dir):
            return
        for fname in os.listdir(output_dir):
            if fname == "README.md":
                continue
            if not (fname.endswith(".cgns") or fname.endswith(".plt")):
                continue
            full_path = os.path.join(output_dir, fname)
            if os.path.isfile(full_path):
                os.remove(full_path)

    def test_unsteady_file_names(self):
        output_dir = self.options["outputDirectory"]
        before = set(os.listdir(output_dir)) if os.path.isdir(output_dir) else set()

        try:
            CFDSolver = ADFLOW(options=self.options)
            ap = copy.copy(ap_tutorial_wing)
            CFDSolver.setAeroProblem(ap)
            CFDSolver(ap)

            after = set(os.listdir(output_dir)) if os.path.isdir(output_dir) else set()
            new_files = sorted(after - before)
            for fname in new_files:
                if fname == "README.md":
                    continue
                full_path = os.path.join(output_dir, fname)
                if not os.path.isfile(full_path):
                    continue
                self.assertTrue(
                    fname.endswith(".cgns") or fname.endswith(".plt"),
                    msg=f"Unexpected output file extension: {fname}",
                )
        finally:
            self._cleanup_output_files(output_dir)


if __name__ == "__main__":
    unittest.main()
