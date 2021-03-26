import unittest
import os
from adflow import ADFLOW


baseDir = os.path.dirname(os.path.abspath(__file__))


class BasicTests(unittest.TestCase):
    N_PROCS = 1

    def test_import(self):
        gridFile = "input_files/mdo_tutorial_euler.cgns"
        options = {"gridfile": os.path.join(baseDir, "../../", gridFile)}
        ADFLOW(options=options, debug=False)
