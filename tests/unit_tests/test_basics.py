from __future__ import print_function
import unittest
import numpy

import os
import sys
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from adflow import ADFLOW


class BasicTests(unittest.TestCase):
    N_PROCS = 1

    def test_import(self):
        gridFile = 'input_files/mdo_tutorial_euler.cgns'
        options = {'gridfile': gridFile}
        CFDSolver = ADFLOW(options=options, debug=False)
