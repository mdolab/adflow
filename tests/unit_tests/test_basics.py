from __future__ import print_function
import unittest
import numpy

class BasicTests(unittest.TestCase):

    def test_import(self):
        from adflow import ADFLOW
        gridFile = 'input_files/mdo_tutorial_euler.cgns'
        options = {'gridfile': gridFile}
        CFDSolver = ADFLOW(options=options)
