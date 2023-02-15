import unittest
import os
from baseclasses.testing import BaseRegTest

# this is based on
# https://stackoverflow.com/questions/1323455/python-unit-test-with-base-and-sub-class
# another options is to use a class factory like this answer
# https://stackoverflow.com/questions/48234672/how-to-use-same-unit-test-for-different-implementations-in-python
# or yet another is to define a new classed like a parameterized class
# https://eli.thegreenplace.net/2011/08/02/python-unit-testing-parametrized-test-cases
baseDir = os.path.dirname(os.path.abspath(__file__))
refDir = os.path.join(baseDir, "refs")


class RegTest(unittest.TestCase):
    # if True, do not run the inherited train() function
    no_train = False

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
        if not hasattr(self, "handler") or self.no_train:
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
            options = self.CFDSolver.getOptions()
            # this is the root adflow directory, we use this so paths are always taken
            # relative to this instead of cwd
            startDir = os.path.join(baseDir, "../../")
            for key in ["outputDirectory", "gridFile", "restartFile"]:
                if key in options and options[key] is not None:
                    options[key] = os.path.relpath(options[key], start=startDir)

            self.handler.add_metadata(options)

        self.handler.writeRef()


# all complex tests should inherit this class instead of RegTest
class CmplxRegTest(RegTest):
    # complex tests have no training method, setting this attribute to True skips the train() function
    # defined in the parent class
    no_train = True
