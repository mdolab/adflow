import copy
import os
import sys
import unittest
from parameterized import parameterized_class
from adflow import ADFLOW


baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir, "../reg_tests"))

from reg_aeroproblems import ap_tutorial_wing  # noqa E402

# Set to True to keep test output files on disk for post-run inspection.
KEEP_OUTPUT_FILES = False


class _OutputFilesTestBase(unittest.TestCase):
    """Base class that isolates each test's outputs in its own subdirectory of output_files/."""

    def setUp(self):
        self.output_dir = os.path.join(baseDir, "../output_files", self.id().replace(".", "_"))
        os.makedirs(self.output_dir, exist_ok=True)
        self._remove_output_files()

    def tearDown(self):
        if KEEP_OUTPUT_FILES:
            return
        self._remove_output_files()
        try:
            os.rmdir(self.output_dir)
        except OSError:
            pass

    def _remove_output_files(self):
        if not os.path.isdir(self.output_dir):
            return
        for fname in os.listdir(self.output_dir):
            full = os.path.join(self.output_dir, fname)
            if os.path.isfile(full):
                os.remove(full)

    def _files_in_output(self):
        if not os.path.isdir(self.output_dir):
            return set()
        return {f for f in os.listdir(self.output_dir) if os.path.isfile(os.path.join(self.output_dir, f))}

    def assertOutputFiles(self, expected):
        self.assertSetEqual(self._files_in_output(), set(expected))


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
class TestSolutionFileNames(_OutputFilesTestBase):
    N_PROCS = 1

    def setUp(self):
        super().setUp()
        self.options = {
            "outputDirectory": self.output_dir,
            "writeSolutionDigits": self.writeSolutionDigits,
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "writeSurfaceSolution": False,
        }

    def test_failed_mesh(self):
        CFDSolver = ADFLOW(options=self.options)
        ap = copy.copy(ap_tutorial_wing)

        # Pretend mesh warping failed
        CFDSolver.adflow.killsignals.fatalfail = True
        CFDSolver(ap)

        digits = self.writeSolutionDigits
        self.assertOutputFiles({f"failed_mesh_{ap.name}_{'0' * digits}.cgns"})

    def test_volume_solution(self):
        CFDSolver = ADFLOW(options=self.options)
        ap = copy.copy(ap_tutorial_wing)
        CFDSolver.setAeroProblem(ap)
        CFDSolver.curAP.adflowData.callCounter = 5

        CFDSolver.writeSolution()

        digits = self.writeSolutionDigits
        self.assertOutputFiles({f"{ap.name}_{5:0{digits}d}_vol.cgns"})


@parameterized_class(
    [
        {
            "name": "unsteady_three_digits",
            "writeSolutionDigits": 3,
        },
    ]
)
class TestSolutionFileNamesUnsteady(_OutputFilesTestBase):
    N_PROCS = 1

    def setUp(self):
        super().setUp()
        self.options = {
            "outputDirectory": self.output_dir,
            "writeSolutionDigits": self.writeSolutionDigits,
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "writeSurfaceSolution": True,
            "writeTecplotSurfaceSolution": True,
            "writeVolumeSolution": True,
            "nCycles": 1,
            "nTimeStepsFine": 2,
            "equationMode": "unsteady",
        }

    def _expected_output_files(self, ap_name):
        base_name = f"{ap_name}_000"
        expected = {
            f"{base_name}_surf_Timestep0001.cgns",
            f"{base_name}_surf_Timestep0001.plt",
            f"{base_name}_surf_Timestep0002.cgns",
            f"{base_name}_surf_Timestep0002.plt",
        }
        for timestep in range(self.options["nTimeStepsFine"] + 1):
            expected.add(f"{base_name}_vol_Timestep{timestep:04d}.cgns")
        return expected

    def _expected_large_timestep_output_files(self, ap_name):
        base_name = f"{ap_name}_000"
        return {
            f"{base_name}_surf_Timestep10000.cgns",
            f"{base_name}_surf_Timestep10000.plt",
            f"{base_name}_vol_Timestep10000.cgns",
            f"{base_name}_vol_Timestep9999.cgns",
        }

    def test_unsteady_file_names(self):
        CFDSolver = ADFLOW(options=self.options)
        ap = copy.copy(ap_tutorial_wing)
        CFDSolver.setAeroProblem(ap)
        CFDSolver(ap)

        self.assertOutputFiles(self._expected_output_files(ap.name))

    def test_unsteady_file_names_large_timestep(self):
        CFDSolver = ADFLOW(options=self.options)
        ap = copy.copy(ap_tutorial_wing)
        CFDSolver.setAeroProblem(ap)
        CFDSolver.curAP.adflowData.callCounter = 0
        CFDSolver.adflow.monitor.timestepunsteady = 10000

        CFDSolver.writeSolution()

        self.assertOutputFiles(self._expected_large_timestep_output_files(ap.name))


if __name__ == "__main__":
    unittest.main()
