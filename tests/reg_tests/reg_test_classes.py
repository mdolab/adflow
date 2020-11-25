import unittest
import os

import reg_test_utils as utils
from baseclasses import BaseRegTest

# this is based on
# https://stackoverflow.com/questions/1323455/python-unit-test-with-base-and-sub-class
# another options is to use a class factory like this answer
# https://stackoverflow.com/questions/48234672/how-to-use-same-unit-test-for-different-implementations-in-python
# or yet another is to define a new classed like a parameterized class
# https://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases
baseDir = os.path.dirname(os.path.abspath(__file__))

refDir = os.path.join(baseDir, "refs")


class test_objects:
    class RegTest(unittest.TestCase):
        def setUp(self):
            ref_file = os.path.join(refDir, self.ref_file)
            self.handler = BaseRegTest(ref_file)

        def train(self):
            if not hasattr(self, "name"):
                # return immediately when the train method is being called on the based class and NOT the
                # classes created using parametrized
                # this will happen when testing, but will hopefully be fixed down the line
                return

            self.handler.train = True
            self.handler.setRef({})

            # get all of the testing methods
            # all of the tests written in this framework must start with "test_"
            tests = [x for x in dir(self) if x.startswith("test_")]
            for test in tests:
                test_func = getattr(self, test)
                test_func()

            self.handler.writeRef(os.path.join(refDir, self.ref_file))


class RegTest(object):
    def train(self):
        self.handler.setMode(train=True)
        self.handler.setRef({})

        # get all of the testing methods
        # all of the tests written in this framework must start with "test_"
        tests = [x for x in dir(self) if x.startswith("test_")]
        for test in tests:
            test_func = getattr(self, test)
            test_func()

        self.handler.writeRef(os.path.join(refDir, self.ref_file))


class SolveRegTest(RegTest):
    """this is a base class for use with in the mdolab regression tests

    in the parent class you must set
        self.CFDSolver
        self.ap
        self.handler

    """

    def test_solve(self):

        # do the solve
        self.CFDSolver(self.ap)

        # check its accuracy
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-9)
        utils.assert_states_allclose(self.handler, self.CFDSolver, tol=1e-10)
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)


class FunctionalsRegTest(RegTest):
    """this a a base class for use with in the mdolab regression tests

    in the parent class you must set
        self.CFDSolver
        self.ap
        self.handler

    AND call getResidual()

    """

    def test_restart_read(self):
        utils.assert_problem_size_equal(self.handler, self.CFDSolver, tol=1e-10)
        utils.assert_states_allclose(self.handler, self.CFDSolver, tol=1e-10)

    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_functions(self):
        utils.assert_functions_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-9)

    def test_forces_and_tractions(self):
        utils.assert_forces_allclose(self.handler, self.CFDSolver, tol=1e-10)
        utils.assert_tractions_allclose(self.handler, self.CFDSolver, tol=1e-10)

        # Reset the option
        self.CFDSolver.setOption("forcesAsTractions", True)

        # Make sure we can write the force file.
        forces_file = os.path.join(self.CFDSolver.getOption("outputdirectory"), "forces.txt")
        self.CFDSolver.writeForceFile(forces_file)

    # ------------------- Derivative routine checks ----------------------------
    def test_jac_vec_prod_fwd(self):
        utils.assert_fwd_mode_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_jac_vec_prod_bwd(self):
        utils.assert_bwd_mode_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_dot_products(self):
        utils.assert_dot_products_allclose(self.handler, self.CFDSolver, tol=1e-10)


class AdjointRegTest(RegTest):
    # Standard test for solving multiple adjoints and going right back
    # to the DVs. This solves for whatever functions are in the
    # aeroProblem.
    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_adjoint(self):
        utils.assert_adjoint_sens_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)
