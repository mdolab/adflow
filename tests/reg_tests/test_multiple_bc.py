import copy
import os
import unittest

import numpy as np
import reg_test_classes
import reg_test_utils as utils
from reg_aeroproblems import ap_multi_bc
from reg_default_options import adflowDefOpts

from adflow import ADFLOW, ADFLOW_C

baseDir = os.path.dirname(os.path.abspath(__file__))


class MultiBCTests(reg_test_classes.RegTest):
    """
    Tests for the multi_bc problem
    """

    N_PROCS = 2
    ref_file = "multi_bc_tests.json"

    def setUp(self):
        super().setUp()

        self.options = {
            "gridFile": os.path.join(baseDir, "../../input_files/multi_bc.cgns"),
            "writevolumesolution": False,
            "writesurfacesolution": False,
            "writetecplotsurfacesolution": False,
            "mgcycle": "sg",
            "ncycles": 2000,
            "useanksolver": True,
            "usenksolver": False,
            "anksecondordswitchtol": 1e-4,
            "ankcoupledswitchtol": 1e-4,
            "ankinnerpreconits": 3,
            "ankouterpreconits": 3,
            "ankpcilufill": 2,
            "ankasmoverlap": 2,
            "anklinresmax": 0.1,
            "volumevariables": ["temp", "mach", "resrho"],
            "surfacevariables": ["temp", "vx", "vy", "vz", "p", "ptloss", "mach", "rho"],
            "equationType": "RANS",
            "l2convergence": 1e-15,
            "L2ConvergenceRel": 1e-15,
            "adjointl2convergence": 1e-15,
            "adjointL2ConvergenceRel": 1e-15,
            "ankcfllimit": 1e5,
            "adpc": True,
            "innerpreconits": 2,
            "outerpreconits": 3,
            "asmoverlap": 2,
            "adjointsubspacesize": 100,
            "adjointmaxiter": 100,
            "skipafterfailedadjoint": False,
            "smoother": "DADI",
            "solutionPrecision": "double",
        }

        self.options["outputdirectory"] = adflowDefOpts["outputdirectory"]

        CFDSolver = ADFLOW(options=self.options)

        families = ["outflow", "inflow1", "inflow2"]
        funcs = [
            "mdot",
            "aavgptot",
            "mavgttot",
            "aavgps",
            "area",
            "mavgvx",
            "forcezpressure",
            "forcezmomentum",
            "mavgvi",
        ]

        evalFuncs = []
        for family in families:
            for func in funcs:
                CFDSolver.addFunction(func, family, name=f"{func}_{family}")
                evalFuncs.append(f"{func}_{family}")

        self.CFDSolver = CFDSolver

        # This is imported from reg_aeroproblems utility script
        self.ap = ap_multi_bc

        self.ap.setBCVar("Pressure", 75000.0, "outflow")
        self.ap.setBCVar("PressureStagnation", 107400.0, "inflow1")
        self.ap.setBCVar("PressureStagnation", 105400.0, "inflow2")
        self.ap.setBCVar("TemperatureStagnation", 330.0, "inflow1")
        self.ap.setBCVar("TemperatureStagnation", 310.0, "inflow2")

        self.ap.addDV("Pressure", value=75000.0, family="outflow", name="ps1", units="Pa")
        self.ap.addDV("PressureStagnation", value=107400.0, family="inflow1", name="pt1", units="Pa")
        self.ap.addDV("PressureStagnation", value=105400.0, family="inflow2", name="pt2", units="Pa")
        self.ap.addDV("TemperatureStagnation", value=330.0, family="inflow1", name="tt1", units="K")
        self.ap.addDV("TemperatureStagnation", value=310.0, family="inflow2", name="tt2", units="K")

    def test_bc_functions(self):
        """
        Test the BC functions
        """
        self.CFDSolver(self.ap)
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-12)

    def test_bc_adjoint(self):
        """
        Test the adjoint sensitivities
        """
        self.ap.setDesignVars(
            {
                "ps1": 75000.0,
                "pt1": 107400.0,
                "pt2": 105400.0,
                "tt1": 330.0,
                "tt2": 310.0,
            }
        )
        self.CFDSolver(self.ap)

        funcs = {}
        self.CFDSolver.checkSolutionFailure(self.ap, funcs)

        # Check if solution failed
        self.assert_solution_failure()

        funcs = {}
        funcsSens = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSens)

        for key in funcsSens:
            self.handler.root_add_dict(f"{key} sens", funcsSens[key], rtol=1e-11, atol=1e-11)


class MultiBCCmplxTests(reg_test_classes.CmplxRegTest):
    """
    Complex step tests for the multi bc problem.
    """

    N_PROCS = 2
    ref_file = "multi_bc_tests.json"
    h = 1e-40

    def setUp(self):
        super().setUp()

        self.options = {
            "gridFile": os.path.join(baseDir, "../../input_files/multi_bc.cgns"),
            "writevolumesolution": False,
            "writesurfacesolution": False,
            "writetecplotsurfacesolution": False,
            "mgcycle": "sg",
            "ncycles": 2000,
            "useanksolver": True,
            "usenksolver": False,
            "anksecondordswitchtol": 1e-4,
            "ankcoupledswitchtol": 1e-4,
            "ankinnerpreconits": 3,
            "ankouterpreconits": 3,
            "ankpcilufill": 2,
            "ankasmoverlap": 2,
            "anklinresmax": 0.1,
            "volumevariables": ["temp", "mach", "resrho"],
            "surfacevariables": ["temp", "vx", "vy", "vz", "p", "ptloss", "mach", "rho"],
            "equationType": "RANS",
            "l2convergence": 1e-13,
            "L2ConvergenceRel": 1e-13,
            "ankcfllimit": 1e9,
            "smoother": "DADI",
            "solutionPrecision": "double",
        }

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        CFDSolver = ADFLOW_C(options=options)

        families = ["outflow", "inflow1", "inflow2"]
        funcs = [
            "mdot",
            "aavgptot",
            "mavgttot",
            "aavgps",
            "area",
            "mavgvx",
            "forcezpressure",
            "forcezmomentum",
            "mavgvi",
        ]

        evalFuncs = []
        for family in families:
            for func in funcs:
                CFDSolver.addFunction(func, family, name=f"{func}_{family}")
                evalFuncs.append(f"{func}_{family}")

        self.CFDSolver = CFDSolver

        # This is imported from reg_aeroproblems utility script
        self.ap = ap_multi_bc

        self.ap.setBCVar("Pressure", 75000.0, "outflow")
        self.ap.setBCVar("PressureStagnation", 107400.0, "inflow1")
        self.ap.setBCVar("PressureStagnation", 105400.0, "inflow2")
        self.ap.setBCVar("TemperatureStagnation", 330.0, "inflow1")
        self.ap.setBCVar("TemperatureStagnation", 310.0, "inflow2")

        self.ap.addDV("PressureStagnation", value=107400.0, family="inflow1", name="pt1", units="Pa")
        self.ap.addDV("PressureStagnation", value=105400.0, family="inflow2", name="pt2", units="Pa")

    def cmplx_test_bc_functions(self):
        """
        Test if the adjoint sensitivities are correct.
        """

        families = ["outflow", "inflow1", "inflow2"]
        funcs = [
            "mdot",
            "aavgptot",
            "mavgttot",
            "aavgps",
            "area",
            "mavgvx",
            "forcezpressure",
            "forcezmomentum",
            "mavgvi",
        ]

        evalFuncs = []
        for family in families:
            for func in funcs:
                evalFuncs.append(f"{func}_{family}")

        aDV = {"pt1": 107400.0, "pt2": 105400.0}
        self.ap.setDesignVars(aDV)

        self.CFDSolver(self.ap)

        self.assert_solution_failure()

        funcs = {}
        funcsSensCS = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=evalFuncs)

        funcs_plus = {}

        for dv in ["pt1", "pt2"]:
            # Save the old dv
            dvsave = aDV[dv]

            # Perturb the design variable
            aDV[dv] += self.h * 1j
            self.ap.setDesignVars(aDV)

            # Call the solver again
            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver(self.ap)

            # Check if solution failed
            self.assert_solution_failure()

            # Save the new funcs
            funcs_plus[dv] = {}
            self.CFDSolver.evalFunctions(self.ap, funcs_plus[dv], evalFuncs=evalFuncs)

            # Compute the sens
            funcsSensCS[dv] = {}
            for key in funcs:
                funcsSensCS[dv][key] = np.imag(funcs_plus[dv][key]) / self.h

            # reset the DV
            aDV[dv] = dvsave

        for dv in ["pt1", "pt2"]:
            for key in funcsSensCS[dv]:
                ref_val = self.handler.db[f"{key} sens"][dv]
                np.testing.assert_allclose(funcsSensCS[dv][key], ref_val, atol=1e-9, rtol=1e-9)


if __name__ == "__main__":
    unittest.main()
