import copy
import glob
import os
import shutil
import sys
import unittest

from mpi4py import MPI
from adflow import ADFLOW


baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir, "../reg_tests"))

from reg_aeroproblems import ap_tutorial_wing  # noqa: E402


class SolutionFileNamesMixin:
    writeSolutionDigits = None
    N_PROCS = 1

    def _comm(self):
        return MPI.COMM_WORLD

    def _rank(self):
        return self._comm().Get_rank()

    def _barrier(self):
        self._comm().Barrier()

    def _bcast(self, value, root=0):
        return self._comm().bcast(value, root=root)

    def _assert_collective_true(self, local_flag, msg):
        ok = self._comm().allreduce(bool(local_flag), op=MPI.LAND)
        self.assertTrue(ok, msg)

    def setUp(self):
        comm = self._comm()
        rank = self._rank()

        self.output_dir = os.path.join(baseDir, "../output_files")
        self.unsteady_output_dir = os.path.join(
            baseDir, f"../output_files/unsteady_name_check_{self.writeSolutionDigits}"
        )

        self.mesh_file = None

        if rank == 0:
            os.makedirs(self.output_dir, exist_ok=True)
            if os.path.isdir(self.unsteady_output_dir):
                shutil.rmtree(self.unsteady_output_dir, ignore_errors=True)
            os.makedirs(self.unsteady_output_dir, exist_ok=True)
        comm.Barrier()

        self.options = {
            "outputDirectory": self.output_dir,
            "writeSolutionDigits": self.writeSolutionDigits,
            "gridFile": os.path.join(
                baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"
            ),
            "writeSurfaceSolution": False,
        }

        self.unsteady_options = {
            "gridFile": os.path.join(
                baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"
            ),
            "outputDirectory": self.unsteady_output_dir,
            "writeSolutionDigits": self.writeSolutionDigits,
            "equationType": "Euler",
            "equationMode": "unsteady",
            "timeIntegrationScheme": "BDF",
            "nTimeStepsFine": 3,
            "deltaT": 1.0e-2,
            "nCycles": 2,
            "MGCycle": "sg",
            "writeVolumeSolution": True,
            "writeSurfaceSolution": True,
            "monitorVariables": ["cpu", "resrho"],
            "printIterations": False,
        }

    def test_failed_mesh(self):
        solver = ADFLOW(options=self.options)
        ap = copy.copy(ap_tutorial_wing)

        solver.adflow.killsignals.fatalfail = True
        solver(ap)

        self._barrier()

        expected_name = f"failed_mesh_{ap.name}_{0:0{self.writeSolutionDigits}d}.cgns"
        self.mesh_file = os.path.join(self.output_dir, expected_name)

        exists = os.path.isfile(self.mesh_file) if self._rank() == 0 else True
        exists = self._bcast(exists, root=0)

        self._assert_collective_true(
            exists,
            f"Expected failed mesh file not found: {self.mesh_file}",
        )

    def test_volume_solution(self):
        solver = ADFLOW(options=self.options)

        ap = copy.copy(ap_tutorial_wing)
        solver.setAeroProblem(ap)
        solver.curAP.adflowData.callCounter = 5
        solver.writeSolution()

        self._barrier()

        expected_name = f"{ap.name}_{5:0{self.writeSolutionDigits}d}_vol.cgns"
        self.mesh_file = os.path.join(self.output_dir, expected_name)

        exists = os.path.isfile(self.mesh_file) if self._rank() == 0 else True
        exists = self._bcast(exists, root=0)

        self._assert_collective_true(
            exists,
            f"Expected volume solution file not found: {self.mesh_file}",
        )

    def test_unsteady_file_naming(self):
        solver = ADFLOW(options=self.unsteady_options, debug=False)
        ap = copy.copy(ap_tutorial_wing)
        solver(ap)

        self._barrier()

        if self._rank() == 0:
            files = sorted(glob.glob(os.path.join(self.unsteady_output_dir, "*Timestep*")))
            basenames = [os.path.basename(f) for f in files]

            ok_files_exist = len(files) > 0
            bad_old_names = [n for n in basenames if ".cgnsTimestep" in n]
            ok_no_old_names = len(bad_old_names) == 0
            bad_non_cgns = [n for n in basenames if not n.endswith(".cgns")]
            ok_all_cgns = len(bad_non_cgns) == 0

            status = {
                "ok_files_exist": ok_files_exist,
                "ok_no_old_names": ok_no_old_names,
                "ok_all_cgns": ok_all_cgns,
                "bad_old_names": bad_old_names,
                "bad_non_cgns": bad_non_cgns,
                "files": basenames,
            }
        else:
            status = None

        status = self._bcast(status, root=0)

        self._assert_collective_true(
            status["ok_files_exist"],
            "No unsteady timestep files were written.",
        )
        self._assert_collective_true(
            status["ok_no_old_names"],
            f"Found old-style names (.cgnsTimestep): {status['bad_old_names']}",
        )
        self._assert_collective_true(
            status["ok_all_cgns"],
            f"Timestep files must end with .cgns. Found: {status['bad_non_cgns']}",
        )

    def tearDown(self):
        comm = self._comm()
        rank = self._rank()

        comm.Barrier()

        if rank == 0:
            if self.mesh_file is not None and os.path.isfile(self.mesh_file):
                os.remove(self.mesh_file)

            if os.path.isdir(self.unsteady_output_dir):
                shutil.rmtree(self.unsteady_output_dir, ignore_errors=True)

        comm.Barrier()


class TestSolutionFileNamesThreeDigits(SolutionFileNamesMixin, unittest.TestCase):
    writeSolutionDigits = 3


class TestSolutionFileNamesFiveDigits(SolutionFileNamesMixin, unittest.TestCase):
    writeSolutionDigits = 5


if __name__ == "__main__":
    unittest.main()