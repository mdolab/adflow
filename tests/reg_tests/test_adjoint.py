# built-ins
import unittest
import numpy
import os
import copy
from collections import defaultdict
from parameterized import parameterized_class
from mpi4py import MPI

# MACH classes
from pygeo import DVGeometry
from pyspline import Curve
from idwarp import USMesh

# from pywarp import MBMesh

# MACH testing class
from adflow import ADFLOW

from adflow import ADFLOW_C
from idwarp import USMesh_C


import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs, IDWarpDefOpts

from reg_aeroproblems import ap_tutorial_wing, ap_tutorial_wing_laminar, ap_tutorial_wing_rotating
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))


def getDVGeo(ffdFile, isComplex=False):
    # Setup geometry/mesh
    DVGeo = DVGeometry(ffdFile, isComplex=isComplex)

    nTwist = 6
    DVGeo.addRefAxis(
        "wing",
        Curve(
            x=numpy.linspace(5.0 / 4.0, 1.5 / 4.0 + 7.5, nTwist),
            y=numpy.zeros(nTwist),
            z=numpy.linspace(0, 14, nTwist),
            k=2,
        ),
    )

    def twist(val, geo):
        for i in range(nTwist):
            geo.rot_z["wing"].coef[i] = val[i]

    def span(val, geo):
        C = geo.extractCoef("wing")
        s = geo.extractS("wing")
        for i in range(len(C)):
            C[i, 2] += s[i] * val[0]
        geo.restoreCoef(C, "wing")

    DVGeo.addGlobalDV("twist", [0] * nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGlobalDV("span", [0], span, lower=-10, upper=10, scale=1.0)
    DVGeo.addLocalDV("shape", lower=-0.5, upper=0.5, axis="y", scale=10.0)

    return DVGeo


test_params = [
    # # Tutorial scalar JST
    {
        "name": "euler_scalar_JST_tut_wing_1core",
        "options": {
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "restartFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "L2Convergence": 1e-14,
            "monitorVariables": ["cpu", "resrho", "totalr"],
            "adjointL2Convergence": 1e-14,
            "MGCycle": "2w",
            "nCycles": 1000,
            "nCyclesCoarse": 250,
            "useNKsolver": True,
            "NKSwitchTol": 2.5e-4,
            "ANKSwitchTol": 1e-2,
            "ANKSecondordSwitchTol": 1e-2,
            "useBlockettes": False,
            "blockSplitting": False,
        },
        "ref_file": "adjoint_euler_scalar_jst_tut_wing.json",
        "aero_prob": ap_tutorial_wing,
        "evalFuncs": ["cl", "cd"],
        "N_PROCS": 1,
    },
    # Tutorial scalar JST
    {
        "name": "euler_scalar_JST_tut_wing",
        "options": {
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "restartFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
            "L2Convergence": 1e-14,
            "monitorVariables": ["cpu", "resrho", "totalr"],
            "adjointL2Convergence": 1e-14,
            "MGCycle": "2w",
            "nCycles": 1000,
            "nCyclesCoarse": 250,
            "useNKsolver": True,
            "NKSwitchTol": 2.5e-4,
            "ANKSwitchtol": 1e-2,
            "ANKSecondordSwitchTol": 1e-2,
            "useBlockettes": False,
            "blockSplitting": False,
        },
        "ref_file": "adjoint_euler_scalar_jst_tut_wing.json",
        "aero_prob": ap_tutorial_wing,
        "evalFuncs": ["cl", "cd", "cmz", "lift", "drag"],
    },
    # # Tutorial wing laminar
    {
        "name": "laminar_tut_wing",
        "options": {
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_viscous_scalar_jst.cgns"),
            "restartFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_viscous_scalar_jst.cgns"),
            "L2Convergence": 1e-15,
            "L2ConvergenceCoarse": 1e-2,
            "monitorVariables": ["cpu", "resrho", "totalr"],
            "adjointL2Convergence": 1e-16,
            "nCycles": 500,
            "CFL": 1.5,
            "CFLCoarse": 1.25,
            "MGCycle": "2w",
            "nCyclesCoarse": 250,
            "useNKsolver": True,
            "ANKSwitchTol": 1e-2,
            "ANKSecondordSwitchTol": 1e-2,
            "NKSwitchTol": 1e-2,
            "equationType": "laminar NS",
            "useBlockettes": False,
        },
        "ref_file": "adjoint_laminar_tut_wing.json",
        "aero_prob": ap_tutorial_wing_laminar,
        "evalFuncs": ["cl", "cd", "cmz", "lift", "drag"],
    },
    # # Tutorial wing RANS
    {
        "name": "rans_tut_wing",
        "options": {
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
            "restartFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
            "MGCycle": "2w",
            "equationType": "RANS",
            "smoother": "DADI",
            "CFL": 1.5,
            "CFLCoarse": 1.25,
            "resAveraging": "never",
            "nSubiter": 3,
            "nSubiterTurb": 3,
            "nCyclesCoarse": 100,
            "nCycles": 1000,
            "monitorVariables": ["cpu", "resrho", "resturb", "totalr"],
            "useNKsolver": True,
            "ANKSwitchTol": 1e-2,
            "ANKSecondordSwitchTol": 1e-2,
            "L2Convergence": 1e-15,
            "NKSwitchTol": 1e-5,
            "adjointL2Convergence": 1e-16,
            "blockSplitting": True,
            "NKjacobianlag": 2,
            "computecavitation": True,
        },
        "ref_file": "adjoint_rans_tut_wing.json",
        "aero_prob": ap_tutorial_wing,
        "evalFuncs": ["fx", "mz", "cl", "cd", "cmz", "lift", "drag", "cavitation", "colx", "coly", "colz"],
    },
    # # Rotating frame test
    {
        "name": "Rotating_wing",
        "options": {
            "gridFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_rotframe.cgns"),
            "restartFile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_rotframe.cgns"),
            "nCycles": 20000,
            "equationType": "RANS",
            "smoother": "DADI",
            "useRotationSA": True,
            "useQCR": True,
            "useBlockettes": False,
            "restrictionRelaxation": 1.0,
            "CFL": 1.5,
            "NKinnerpreconits": 2,
            "NKjacobianlag": 3,
            "NKouterpreconits": 3,
            "NKsubspacesize": 100,
            "ANKmaxiter": 80,
            "ANKasmoverlap": 2,
            "ANKinnerpreconits": 2,
            "ANKouterpreconits": 2,
            "nsubiterturb": 10,
            "ANKSecondordSwitchTol": 1e-2,
            "ANKcoupledSwitchTol": 1e-8,
            "ANKCFLlimit": 1e10,
            "L2Convergence": 1e-14,
            "adjointL2Convergence": 1e-14,
            "ADPC": True,
            "adjointMaxIter": 1000,
            "adjointSubspaceSize": 500,
            "monitorVariables": ["cpu", "resrho", "resturb", "cmx"],
            "volumevariables": ["resrho"],
        },
        "ref_file": "adjoint_rans_rotating.json",
        "aero_prob": ap_tutorial_wing_rotating,
        "evalFuncs": ["fy", "my"],
        "N_PROCS": 2,
    },
]


@parameterized_class(test_params)
class TestAdjoint(reg_test_classes.RegTest):
    """
    Tests that sensitives calculated from solving an adjoint are correct.
    and jacobian vector products are accurate.

    based on old regression tests 12, and 14
    """

    N_PROCS = 2

    options = None
    ap = None
    ref_file = None

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.ffdFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_ffd.fmt")

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridFile"]})

        self.ap = copy.deepcopy(self.aero_prob)

        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

        # add the default dvs to the problem
        if self.name != "Rotating_wing":
            for dv in defaultAeroDVs:
                self.ap.addDV(dv)
        else:
            for dv in ["alpha", "beta", "mach", "T", "xRef", "yRef", "zRef"]:
                self.ap.addDV(dv)

        self.CFDSolver = ADFLOW(options=options, debug=True)

        self.CFDSolver.setMesh(USMesh(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=False), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})
        if self.name == "Rotating_wing":
            # Add rotation component to the frame
            rotRate_x = 0.5
            self.CFDSolver.setRotationRate([0, 0, 0], [0, -rotRate_x, 0])

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def test_residuals(self):
        utils.assert_residuals_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)

    def test_adjoint(self):
        utils.assert_adjoint_sens_allclose(self.handler, self.CFDSolver, self.ap, tol=1e-10)
        self.assert_adjoint_failure()


@parameterized_class(test_params)
class TestCmplxStep(reg_test_classes.CmplxRegTest):
    """
    Tests that sensitives calculated from solving an adjoint are correct.
    and jacobian vector products are accurate.

    based on old regression tests_cs 12, and 14
    """

    N_PROCS = 2

    h = 1e-40

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return
        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        self.ffdFile = os.path.join(baseDir, "../../input_files/mdo_tutorial_ffd.fmt")

        mesh_options = copy.copy(IDWarpDefOpts)
        mesh_options.update({"gridFile": options["gridFile"]})

        self.ap = copy.deepcopy(self.aero_prob)

        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

        # add the default dvs to the problem
        if self.name != "Rotating_wing":
            for dv in defaultAeroDVs:
                self.ap.addDV(dv)
        else:
            for dv in ["alpha", "beta", "mach", "T", "xRef", "yRef", "zRef"]:
                self.ap.addDV(dv)

        self.CFDSolver = ADFLOW_C(options=options, debug=True)

        self.CFDSolver.setMesh(USMesh_C(options=mesh_options))
        self.CFDSolver.setDVGeo(getDVGeo(self.ffdFile, isComplex=True), pointSetKwargs={"embTol": 1e-12, "eps": 1e-14})
        if self.name == "Rotating_wing":
            # Add rotation component to the frame
            rotRate_x = 0.5
            self.CFDSolver.setRotationRate([0, 0, 0], [0, -rotRate_x, 0])

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def cmplx_test_aero_dvs(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        if self.name == "Rotating_wing":
            rtol = 2e-6
        else:
            rtol = 1e-8
        atol = 5e-10

        for dv in ["alpha", "mach"]:  # defaultAeroDVs:
            funcsSens = defaultdict(lambda: {})
            setattr(self.ap, dv, getattr(self.ap, dv) + self.h * 1j)

            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver(self.ap, writeSolution=False)
            self.assert_solution_failure()

            funcs = {}
            self.CFDSolver.evalFunctions(self.ap, funcs)
            setattr(self.ap, dv, getattr(self.ap, dv) - self.h * 1j)

            for f in self.ap.evalFuncs:
                key = self.ap.name + "_" + f
                dv_key = dv + "_" + self.ap.name
                funcsSens[key][dv_key] = numpy.imag(funcs[key]) / self.h

        if MPI.COMM_WORLD.rank == 0:
            print("====================================")
            print(self.ap.alpha)
            print(self.ap.mach)
            print(self.name, funcsSens)
            print("====================================")

        self.handler.root_add_dict("Eval Functions Sens:", funcsSens, rtol=rtol, atol=atol)

    def cmplx_test_geom_dvs(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # redo the setup for a cmplx test
        funcsSens = defaultdict(lambda: {})

        xRef = {"twist": [0.0] * 6, "span": [0.0], "shape": numpy.zeros(72, dtype="D")}

        if self.name == "Rotating_wing":
            rtol = 4e-6
        else:
            rtol = 5e-9
        atol = 5e-9

        for dv in ["span", "twist", "shape"]:
            xRef[dv][0] += self.h * 1j

            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver.DVGeo.setDesignVars(xRef)
            self.CFDSolver(self.ap, writeSolution=False)
            self.assert_solution_failure()

            funcs = {}
            self.CFDSolver.evalFunctions(self.ap, funcs)

            xRef[dv][0] -= self.h * 1j

            for f in self.ap.evalFuncs:
                key = self.ap.name + "_" + f
                dv_key = dv
                funcsSens[key][dv_key] = numpy.imag(funcs[key]) / self.h

                err_msg = "Failed value for: {}".format(key + " " + dv_key)

                ref_val = self.handler.db["Eval Functions Sens:"][key][dv_key]
                if not isinstance(ref_val, float):
                    ref_val = ref_val.flatten()[0]

                numpy.testing.assert_allclose(funcsSens[key][dv_key], ref_val, atol=atol, rtol=rtol, err_msg=err_msg)

        if MPI.COMM_WORLD.rank == 0:
            print("====================================")
            print(self.name, funcsSens)
            print("====================================")


if __name__ == "__main__":
    unittest.main()
