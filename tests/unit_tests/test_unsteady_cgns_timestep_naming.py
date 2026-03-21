import copy
import glob
import os
import shutil
import sys
import unittest

from adflow import ADFLOW
from petsc4py import PETSc


baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir, "../reg_tests"))

from reg_aeroproblems import ap_tutorial_wing  # noqa E402


class TestUnsteadyCGNSTimestepNaming(unittest.TestCase):
    """
    Check that unsteady CGNS outputs use the new naming:
      <name>_TimestepXXXX.cgns
    and do not use the old (broken) naming:
      <name>.cgnsTimestepXXXX
    """

    N_PROCS = 1

    def setUp(self):
        self.output_dir = os.path.join(baseDir, "../output_files/unsteady_name_check")
        os.makedirs(self.output_dir, exist_ok=True)

        self.options = {
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "outputDirectory": self.output_dir,
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

    def test_unsteady_file_naming(self):
        solver = ADFLOW(options=self.options, debug=False)
        ap = copy.copy(ap_tutorial_wing)
        solver(ap)

        files = sorted(glob.glob(os.path.join(self.output_dir, "*Timestep*")))
        comm = PETSc.COMM_WORLD
        rank = comm.getRank()
        basenames = [os.path.basename(f) for f in files]

        if rank == 0:
            self.assertTrue(files, "No unsteady timestep files were written.")

        if rank == 0:
            bad_names = [n for n in basenames if ".cgnsTimestep" in n]
            self.assertFalse(
                bad_names,
                msg=f"Found old-style names (.cgnsTimestep): {bad_names}",
            )
            non_cgns = [n for n in basenames if not n.endswith(".cgns")]
            self.assertFalse(
                non_cgns,
                msg=f"Timestep files must end with .cgns. Found: {non_cgns}",
            )

    def tearDown(self):
        comm = PETSc.COMM_WORLD
        rank = comm.getRank()
        comm.Barrier()
        if rank == 0:
            if os.path.isdir(self.output_dir):
                shutil.rmtree(self.output_dir, ignore_errors=True)
        comm.Barrier()


if __name__ == "__main__":
    unittest.main()
