# built-ins
import unittest
import os
import copy
import numpy
from parameterized import parameterized_class

# MACH classes
from adflow import ADFLOW
from adflow import ADFLOW_C

# import the testing utilities that live a few directories up
import reg_test_utils as utils
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_simple_cart_cube
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))

# tests for different solver combinations using a simple cartesian block.
# we are interested in these tests converging without seg-faulting
# and switching between the solvers correctly. We test
# Euler + laminar NS + RANS, with smoother, ANK, SANK, CSANK, NK transition order
# with turbDADI and turbKSP. So 3 eqn types and 2 turbulence solvers, we have 6 tests
# for both real and complex cases. These do not cover every possible combination, but
# represent the most common use cases. Testing turbulence solver options with euler
# and laminar NS makes sure the solver does not segfault etc.

solver_combo_params = [
    # Euler Tests
    {
        "name": "euler_smoother_ank_sank_csank_nk_turbdadi",
        "options": {
            "equationtype": "Euler",
            # euler fails to get to machine zero in parallel...
            "L2Convergence": 1e-8,
            "ankuseturbdadi": True,
        },
    },

    {
        "name": "euler_smoother_ank_sank_csank_nk_turbksp",
        "options": {
            "equationtype": "Euler",
            # euler fails to get to machine zero in parallel...
            "L2Convergence": 1e-8,
            "ankuseturbdadi": False,
        },
    },

    # laminar NS Tests
    {
        "name": "laminar_ns_smoother_ank_sank_csank_nk_turbdadi",
        "options": {
            "equationtype": "laminar NS",
            "ankuseturbdadi": True,
        },
    },

    {
        "name": "laminar_ns_smoother_ank_sank_csank_nk_turbksp",
        "options": {
            "equationtype": "laminar NS",
            "ankuseturbdadi": False,
        },
    },

    # RANS tests
    {
        "name": "rans_smoother_ank_sank_csank_nk_turbdadi",
        "options": {
            "equationtype": "RANS",
            "ankuseturbdadi": True,
        },
    },

    {
        "name": "rans_smoother_ank_sank_csank_nk_turbksp",
        "options": {
            "equationtype": "RANS",
            "ankuseturbdadi": False,
        },
    },

]

# common options dict for both real and complex
gridFile = os.path.join(baseDir, '../../input_files/cube.cgns')
commonTestOptions = {
    'gridfile': gridFile,
    "writevolumesolution": False,
    'writesurfacesolution': False,
    'monitorvariables': ['cpu', 'resrho', 'resturb', 'cd'],
    "mgcycle": "SG",
    "mgstartlevel": -1,
    'l2convergence': 1e-12,
    'ncycles': 500,
    'useblockettes': False,

    # solver flags
    'useanksolver': True,
    "usenksolver": True,

    # switch tolerances
    "ankswitchtol": 0.5,
    "anksecondordswitchtol": 1e-1,
    "ankcoupledswitchtol": 1e-3,
    "nkswitchtol": 1e-5,
}

@parameterized_class(solver_combo_params)
class TestSolverCombos(reg_test_classes.RegTest):
    N_PROCS = 2
    ref_file = 'test_solver_combos.json'

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        super().setUp()

        # start with the default options dictionary
        options = copy.copy(adflowDefOpts)

        # set the output directory
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])

        # these are the modified options common to these tests
        options.update(commonTestOptions)

        # finally, bring in the specific options for each parameterized test
        options.update(self.options)

        self.ap = copy.deepcopy(ap_simple_cart_cube)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options, debug=False)

    def test_convergence(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # do the solve
        self.CFDSolver(self.ap)

        # get residual norms
        r0, _, rfinal = self.CFDSolver.getResNorms()

        # get the target
        l2conv = self.CFDSolver.getOption("L2Convergence")

        # we should get 12 orders of magnitude relative convergence
        numpy.testing.assert_array_less(rfinal / r0, l2conv)

# we do the same tests as the previous one but with complex mode here
@parameterized_class(solver_combo_params)
class TestCmplxSolverCombos(reg_test_classes.RegTest):

    # TODO add a convergence test with a complex perturbed DV also

    N_PROCS = 2
    ref_file = 'test_solver_combos.json'

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        super().setUp()

        # start with the default options dictionary
        options = copy.copy(adflowDefOpts)

        # set the output directory
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])

        # these are the modified options common to these tests
        options.update(commonTestOptions)

        # finally, bring in the specific options for each parameterized test
        options.update(self.options)

        self.ap = copy.deepcopy(ap_simple_cart_cube)

        # Create the solver
        self.CFDSolver = ADFLOW_C(options=options, debug=False)

    def cmplx_test_convergence(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # do the solve
        self.CFDSolver(self.ap)

        # get residual norms
        r0, _, rfinal = self.CFDSolver.getResNorms()

        # get the target
        l2conv = self.CFDSolver.getOption("L2Convergence")

        # we should get 12 orders of magnitude relative convergence
        numpy.testing.assert_array_less(rfinal / r0, l2conv)

if __name__ == "__main__":
    unittest.main()
