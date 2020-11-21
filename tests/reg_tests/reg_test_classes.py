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
            self.handler = BaseRegTest(ref_file, train=False)

        def train(self):
            if not hasattr(self, "handler"):
                # return immediately when the train method is being called on the based class and NOT the
                # classes created using parametrized
                # this will happen when testing, but will hopefully be fixed down the line
                return

            self.handler.setTrainingMode()

            # get all of the testing methods
            # all of the tests written in this framework must start with "test_"
            tests = [x for x in dir(self) if x.startswith("test_")]
            for test in tests:
                test_func = getattr(self, test)
                test_func()

            self.handler.writeRef()
