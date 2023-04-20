import unittest
import os
import numpy as np
from adflow import ADFLOW
from adflow import ADFLOW_C


baseDir = os.path.dirname(os.path.abspath(__file__))


class BasicTests(unittest.TestCase):
    N_PROCS = 1

    def test_import(self):
        gridFile = "input_files/mdo_tutorial_euler.cgns"
        options = {"gridfile": os.path.join(baseDir, "../../", gridFile)}
        ADFLOW(options=options, debug=False)

    def verifyOptions(self, CFDSolver):
        pythonOptions, deprecatedOptions, specialOptions = CFDSolver._getSpecialOptionLists()

        for optionName in CFDSolver.defaultOptions:
            optionName = optionName.lower()

            if optionName in pythonOptions or optionName in deprecatedOptions or optionName in specialOptions:
                # Ignore for now
                continue

            # Get the python value from the dict
            value = CFDSolver.getOption(optionName)

            # Treat all other options
            if isinstance(CFDSolver.optionMap[optionName], dict):
                # For values that are set in fortran using constant from fortran
                # we just overwrite the python value. An example is "solutionprecision" which is
                # a str value in python that corresponds to an integer in fortran, which sets a module variable precisionsol
                module = CFDSolver.moduleMap[CFDSolver.optionMap[optionName]["location"][0]]
                variable = CFDSolver.optionMap[optionName]["location"][1]
                value = CFDSolver.optionMap[optionName][value.lower()]
            else:
                module = CFDSolver.moduleMap[CFDSolver.optionMap[optionName][0]]
                variable = CFDSolver.optionMap[optionName][1]

            # Fetch the variable value form the module
            fortranValue = getattr(module, variable)

            # Now check if the set parameter is correct in fortran
            if type(value) is bool:
                fortranValue = bool(fortranValue)

            elif type(value) is str:
                # To not mess string arrays too much the following converts string byte array to regular array and strip whitespace
                fortranValue = str(np.array(fortranValue, dtype=str)).strip()

            elif optionName == "mgstartlevel" and (value == -1 or value == 0):
                # mgstartlevel is set in Fortran but the default value (-1) defaults to the coarsest level
                # value of -1 or 0 should default to coarsest level, which is extracted based on the number of MG levels strategy

                # Just overwrite value with what is set in fortran
                value = CFDSolver.adflow.inputiteration.nmglevels

            # Note/todo: When running in complex mode the real part is compared to the default value. Options should in most cases always be real, but are complex in fortran. This should be updated at some point.
            self.assertEqual(value, fortranValue)

    def cmplx_test_options(self):
        """
        Test that options are set properly from python to Fortran in complex mode
        """
        gridFile = "input_files/cube_4x4x4.cgns"
        options = {"gridfile": os.path.join(baseDir, "../../", gridFile)}
        CFDSolver = ADFLOW_C(options=options, debug=False)
        self.verifyOptions(CFDSolver)

    def test_options(self):
        """
        Test that options are set properly from python to Fortran
        """
        gridFile = "input_files/cube_4x4x4.cgns"
        options = {"gridfile": os.path.join(baseDir, "../../", gridFile)}
        CFDSolver = ADFLOW(options=options, debug=False)
        self.verifyOptions(CFDSolver)


if __name__ == "__main__":
    unittest.main()
