import unittest
import os
import copy
from parameterized import parameterized_class
from adflow import ADFLOW
from adflow import ADFLOW_C
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_naca0012_time_acc

baseDir = os.path.dirname(os.path.abspath(__file__))

gridFile = os.path.join(baseDir, "../../input_files/naca0012_rans-L2.cgns")
commonTestOptions = {
    "gridFile": gridFile,
    "equationType": "RANS",
    "writeVolumeSolution": False,
    "writeSurfaceSolution": False,
    "monitorVariables": ["cpu", "resrho", "resturb", "cd"],
    "MGCycle": "sg",
    "L2Convergence": 1e-12,
    "adjointL2Convergence": 1e-12,
    "nCycles": 2000,
    "useBlockettes": False,
    "useANKSolver": True,
    "useNKSolver": True,
    # Switch tolerances to test ANK, SANK, CSANK, NK
    "ANKSwitchTol": 10.0,
    "ANKSecondOrdSwitchTol": 1e-3,
    "ANKCoupledSwitchTol": 1e-6,
    "NKSwitchTol": 1e-8,
    # Use AMG for all solvers
    "ANKGlobalPreconditioner": "multigrid",
    "NKGlobalPreconditioner": "multigrid",
    "globalPreconditioner": "multigrid",
}

test_params = [
    {
        "name": "two_levels_one_iteration",
        "options": {
            # AMG levels
            "ANKAMGLevels": 2,
            "NKAMGLevels": 2,
            "adjointAMGLevels": 2,
            # Outer iterations
            "ANKOuterPreconIts": 1,
            "NKOuterPreconIts": 1,
            "outerPreconIts": 1,
        },
    },
    {
        "name": "two_levels_two_iterations",
        "options": {
            # AMG levels
            "ANKAMGLevels": 2,
            "NKAMGLevels": 2,
            "adjointAMGLevels": 2,
            # Outer iterations
            "ANKOuterPreconIts": 2,
            "NKOuterPreconIts": 2,
            "outerPreconIts": 2,
        },
    },
    {
        "name": "different_levels_and_iterations",
        "options": {
            # AMG levels
            "ANKAMGLevels": 2,
            "NKAMGLevels": 2,
            "adjointAMGLevels": 1,
            # Outer iterations
            "ANKOuterPreconIts": 1,
            "NKOuterPreconIts": 3,
            "outerPreconIts": 2,
            # Smoothing iterations
            "ANKAMGNSmooth": 2,
            "NKAMGNSmooth": 3,
            "adjointAMGNSmooth": 4,
            # Inner iterations
            "ANKInnerPreconIts": 2,
            "NKInnerPreconIts": 1,
            "innerPreconIts": 3,
        },
    },
    {
        "name": "nondefault_coarse_levels",
        "options": {
            # AMG levels
            "ANKAMGLevels": 2,
            "NKAMGLevels": 2,
            "adjointAMGLevels": 2,
            # ANK coarse level options
            "ANKASMOverlapCoarse": 2,
            "ANKPCILUFillCoarse": 2,
            "ANKInnerPreconItsCoarse": 2,
            # NK coarse level options
            "NKASMOverlapCoarse": 2,
            "NKPCILUFillCoarse": 2,
            "NKInnerPreconItsCoarse": 2,
            # Adjoint coarse level options
            "ASMOverlapCoarse": 2,
            "ILUFillCoarse": 2,
            "innerPreconItsCoarse": 2,
        },
    },
]


@parameterized_class(test_params)
class TestAMGReal(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        # Start with the default testing options dictionary
        options = copy.copy(adflowDefOpts)

        # Set the output directory
        options["outputDirectory"] = os.path.join(baseDir, options["outputDirectory"])

        # These are the modified options common to these tests
        options.update(commonTestOptions)

        # Add the parameterized options
        options.update(self.options)

        # Define the AeroProblem
        self.ap = copy.deepcopy(ap_naca0012_time_acc)

        # Create the solver
        self.CFDSolver = ADFLOW(options=options)

    def test_convergence(self):
        # Solve the flow
        self.CFDSolver(self.ap)

        # Check that the flow converged
        funcs = {}
        self.CFDSolver.checkSolutionFailure(self.ap, funcs)
        self.assertFalse(funcs["fail"])

        # Solve the adjoint
        funcsSens = {}
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSens, evalFuncs=["cd"])

        # Check that the adjoint converged
        self.CFDSolver.checkAdjointFailure(self.ap, funcsSens)
        self.assertFalse(funcsSens["fail"])


class TestAMGCoarseMatrixAssembly(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        self.options = copy.copy(adflowDefOpts)
        self.options["outputDirectory"] = os.path.join(baseDir, self.options["outputDirectory"])
        self.options.update(commonTestOptions)
        self.ap = copy.deepcopy(ap_naca0012_time_acc)

    def test_nk_coarse_matrix_is_reassembled_not_accumulated(self):
        """
        Given a flow solve that ends in the NK solver with a two-level multigrid preconditioner,
        when the NK preconditioner matrix is assembled twice in a row for the same state,
        then the coarse-level matrix has the same nonzero norm after each assembly.
        """
        solver = ADFLOW(options=self.options)
        solver(self.ap)

        norms = []
        for _ in range(2):
            solver.adflow.nksolver.formjacobiannk()
            norms.append(solver.adflow.amg.getcoarsematnorm(2))

        self.assertGreater(norms[0], 0.0)
        self.assertAlmostEqual(norms[1] / norms[0], 1.0, places=12)

    def test_ank_coarse_matrix_is_reassembled_not_accumulated(self):
        """
        Given a flow solve that ends in the ANK solver with a two-level multigrid preconditioner,
        when the ANK time-step terms, and then the full ANK preconditioner, are each assembled twice in a row,
        then the coarse-level matrix has the same nonzero norm after each repeat,
        and the full preconditioner's coarse matrix differs from the time-step terms alone.
        """
        self.options.update({"useNKSolver": False, "L2Convergence": 1e-8})
        solver = ADFLOW(options=self.options)
        solver(self.ap)
        ank = solver.adflow.anksolver

        timeStepNorms = []
        for _ in range(2):
            ank.computetimestepmat(True)
            timeStepNorms.append(solver.adflow.amg.getcoarsematnorm(2))

        fullNorms = []
        for _ in range(2):
            ank.computetimestepmat(True)
            ank.formjacobianank()
            fullNorms.append(solver.adflow.amg.getcoarsematnorm(2))

        self.assertGreater(timeStepNorms[0], 0.0)
        self.assertAlmostEqual(timeStepNorms[1] / timeStepNorms[0], 1.0, places=12)
        self.assertAlmostEqual(fullNorms[1] / fullNorms[0], 1.0, places=12)
        self.assertNotAlmostEqual(fullNorms[0] / timeStepNorms[0], 1.0, places=6)


# Run the flow convergence tests in complex mode
@parameterized_class(test_params)
class TestAMGComplex(unittest.TestCase):
    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # Return immediately when the setUp method is being called on the base class and NOT the
            # classes created using parameterized.
            # This is required for running complex tests with ``testflo -m 'cmplx*'`` but will hopefully be
            # fixed down the line.
            return

        # Start with the default testing options dictionary
        options = copy.copy(adflowDefOpts)

        # Set the output directory
        options["outputDirectory"] = os.path.join(baseDir, options["outputDirectory"])

        # These are the modified options common to these tests
        options.update(commonTestOptions)

        # Add the parameterized options
        options.update(self.options)

        # Define the AeroProblem
        self.ap = copy.deepcopy(ap_naca0012_time_acc)

        # Create the solver
        self.CFDSolver = ADFLOW_C(options=options)

    def cmplx_test_convergence(self):
        if not hasattr(self, "name"):
            # Return immediately when the setUp method is being called on the base class and NOT the
            # classes created using parameterized.
            # This is required for running complex tests with ``testflo -m 'cmplx*'`` but will hopefully be
            # fixed down the line.
            return

        # Add an imaginary perturbation to the angle of attack
        self.ap.alpha += 1e-100 * 1j

        # Solve the flow
        self.CFDSolver(self.ap)

        # Check that the flow converged
        funcs = {}
        self.CFDSolver.checkSolutionFailure(self.ap, funcs)
        self.assertFalse(funcs["fail"])


if __name__ == "__main__":
    unittest.main()
