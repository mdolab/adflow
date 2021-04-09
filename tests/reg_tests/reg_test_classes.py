import unittest
import os
import copy
from baseclasses import BaseRegTest

# this is based on
# https://stackoverflow.com/questions/1323455/python-unit-test-with-base-and-sub-class
# another options is to use a class factory like this answer
# https://stackoverflow.com/questions/48234672/how-to-use-same-unit-test-for-different-implementations-in-python
# or yet another is to define a new classed like a parameterized class
# https://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases
baseDir = os.path.dirname(os.path.abspath(__file__))
refDir = os.path.join(baseDir, "refs")


class RegTest(unittest.TestCase):
    def setUp(self):
        ref_file = os.path.join(refDir, self.ref_file)
        train = "train" in self.id().split(".")[-1]
        self.handler = BaseRegTest(ref_file, train=train)

    def assert_solution_failure(self):
        funcs = {}
        self.CFDSolver.checkSolutionFailure(self.ap, funcs)
        self.assertFalse(funcs["fail"])

    def assert_adjoint_failure(self):
        funcsSens = {}
        self.CFDSolver.checkAdjointFailure(self.ap, funcsSens)
        self.assertFalse(funcsSens["fail"])

    def train(self):
        if not hasattr(self, "handler"):
            # return immediately when the train method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when testing, but will hopefully be fixed down the line
            return

        # get all of the testing methods
        # all of the tests written in this framework must start with "test_"
        tests = [x for x in dir(self) if x.startswith("test_")]
        for test in tests:
            test_func = getattr(self, test)
            test_func()

            # Convert paths to relative before storing the options
            opts = copy.deepcopy(self.CFDSolver.getOptions())
            for key in ["outputDirectory", "gridFile", "restartFile"]:
                if key in opts and opts[key] is not None:
                    opts[key] = os.path.relpath(opts[key])

            self.handler.add_metadata(opts)

        self.handler.writeRef()
