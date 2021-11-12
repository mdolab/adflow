# built-ins
import unittest
import os
import copy
from parameterized import parameterized_class
import numpy as np
# MACH classes
from adflow import ADFLOW

# import the testing utilities
import reg_test_utils as utils

from reg_default_options import adflowDefOpts, defaultAeroDVs

from reg_aeroproblems import ap_tutorial_wing, ap_CRM, ap_tutorial_wing_laminar
import reg_test_classes


baseDir = os.path.dirname(os.path.abspath(__file__))

test_params = [
        # scalar JST
        {
            "name": "euler_scalar_jst_tut_wing_1core",
            "options": {
                "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
                "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
                "l2convergence": 1e-14,
                "mgcycle": "2w",
                "ncyclescoarse": 250,
                "usenksolver": True,
                "useblockettes": False,
            },
            "ref_file": "funcs_euler_scalar_jst_tut_wing.json",
            "aero_prob": copy.deepcopy(ap_tutorial_wing),
            "N_PROCS": 1,
        },
        # # scalar JST
        # {
        #     "name": "euler_scalar_jst_tut_wing",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_scalar_jst.cgns"),
        #         "l2convergence": 1e-14,
        #         "mgcycle": "2w",
        #         "ncyclescoarse": 250,
        #         "usenksolver": True,
        #         "useblockettes": False,
        #     },
        #     "ref_file": "funcs_euler_scalar_jst_tut_wing.json",
        #     "aero_prob": ap_tutorial_wing,
        # },
        # # Matrix JST
        # {
        #     "name": "euler_matrix_jst_tut_wing",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_matrix.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_matrix.cgns"),
        #         "mgcycle": "2w",
        #         "ncyclescoarse": 250,
        #         "usenksolver": True,
        #         "nkswitchtol": 1e-2,
        #         "vis4": 0.1,
        #         "discretization": "central plus matrix dissipation",
        #         "coarsediscretization": "central plus matrix dissipation",
        #         "l2convergence": 1e-14,
        #         "useblockettes": False,
        #     },
        #     "ref_file": "funcs_euler_matrix_jst_tut_wing.json",
        #     "aero_prob": ap_tutorial_wing,
        # },
        # # Upwind
        # {
        #     "name": "euler_upwind_tut_wing",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_upwind.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_euler_upwind.cgns"),
        #         "mgcycle": "2w",
        #         "ncyclescoarse": 250,
        #         "usenksolver": True,
        #         "nkswitchtol": 1e-2,
        #         "vis4": 0.1,
        #         "discretization": "upwind",
        #         "useblockettes": False,
        #         "l2convergence": 1e-14,
        #     },
        #     "ref_file": "funcs_euler_upwind_tut_wing.json",
        #     "aero_prob": ap_tutorial_wing,
        # },
        # # Tutorial wing random block order
        # {
        #     "name": "euler_scalar_jst_rand_tut_wing",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_random_euler_scalar_jst.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_random_euler_scalar_jst.cgns"),
        #         "l2convergence": 1e-14,
        #         "mgcycle": "2w",
        #         "ncyclescoarse": 250,
        #         "usenksolver": True,
        #         "useblockettes": False,
        #     },
        #     "ref_file": "funcs_euler_scalar_jst_rand_tut_wing.json",
        #     "aero_prob": ap_tutorial_wing,
        # },
        # # Tutorial wing laminar
        # {
        #     "name": "laminar_tut_wing",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_viscous_scalar_jst.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_viscous_scalar_jst.cgns"),
        #         "l2convergence": 1e-14,
        #         "l2convergencecoarse": 1e-4,
        #         "ncycles": 1000,
        #         "cfl": 1.5,
        #         "cflcoarse": 1.25,
        #         "mgcycle": "2w",
        #         "ncyclescoarse": 250,
        #         "usenksolver": True,
        #         "nkswitchtol": 1e-2,
        #         "equationtype": "laminar NS",
        #         "useblockettes": False,
        #     },
        #     "ref_file": "funcs_laminar_tut_wing.json",
        #     "aero_prob": ap_tutorial_wing_laminar,
        # },
        # # Tutorial wing RANS
        # {
        #     "name": "rans_tut_wing",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_rans_scalar_jst.cgns"),
        #         "mgcycle": "sg",
        #         "equationtype": "RANS",
        #         "smoother": "DADI",
        #         "cfl": 1.5,
        #         "cflcoarse": 1.25,
        #         "resaveraging": "never",
        #         "nsubiter": 3,
        #         "nsubiterturb": 3,
        #         "ncyclescoarse": 100,
        #         "ncycles": 1000,
        #         "monitorvariables": ["cpu", "resrho", "resturb", "cl", "cd", "cmz", "yplus", "totalr"],
        #         "usenksolver": True,
        #         "l2convergence": 1e-14,
        #         "l2convergencecoarse": 1e-4,
        #         "nkswitchtol": 1e-3,
        #         "adjointl2convergence": 1e-14,
        #         "frozenturbulence": False,
        #     },
        #     "ref_file": "funcs_rans_tut_wing.json",
        #     "aero_prob": ap_tutorial_wing,
        # },
        # # Tutorial wing random RANS
        # {
        #     "name": "rans_rand_tut_wing",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_random_rans_scalar_jst.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/mdo_tutorial_random_rans_scalar_jst.cgns"),
        #         "mgcycle": "sg",
        #         "equationtype": "RANS",
        #         "smoother": "DADI",
        #         "cfl": 1.5,
        #         "cflcoarse": 1.25,
        #         "resaveraging": "never",
        #         "nsubiter": 3,
        #         "nsubiterturb": 3,
        #         "ncyclescoarse": 100,
        #         "ncycles": 1000,
        #         "monitorvariables": ["cpu", "resrho", "resturb", "cl", "cd", "cmz", "yplus", "totalr"],
        #         "usenksolver": True,
        #         "l2convergence": 1e-14,
        #         "l2convergencecoarse": 1e-4,
        #         "nkswitchtol": 1e-3,
        #         "adjointl2convergence": 1e-14,
        #         "frozenturbulence": False,
        #     },
        #     "ref_file": "funcs_rans_rand_tut_wing.json",
        #     "aero_prob": ap_tutorial_wing,
        # },
        # # CRM WBT
        # {
        #     "name": "euler_scalar_jst_CRM_WBT",
        #     "options": {
        #         "gridfile": os.path.join(baseDir, "../../input_files/CRM_wbt_scalar_jst.cgns"),
        #         "restartfile": os.path.join(baseDir, "../../input_files/CRM_wbt_scalar_jst.cgns"),
        #         "mgcycle": "sg",
        #         "cfl": 1.5,
        #         "cflcoarse": 1.25,
        #         "resaveraging": "never",
        #         "ncycles": 1000,
        #         "monitorvariables": ["resrho", "cl", "cd", "cmy", "yplus", "totalr"],
        #         "usenksolver": True,
        #         "l2convergence": 1e-14,
        #         "l2convergencecoarse": 1e-4,
        #         "nkswitchtol": 1e-1,
        #         "adjointl2convergence": 1e-14,
        #         "liftindex": 3,
        #         "useblockettes": False,
        #     },
        #     "ref_file": "funcs_euler_scalar_jst_CRM_WBT.json",
        #     "aero_prob": ap_CRM,
        # },
    ]
@parameterized_class(test_params)
class TestJacVecFwd(reg_test_classes.RegTest):
    """
    Tests that given a flow state the residuals, function, forces/tractions,
    and jacobian vector products are accurate.

    """

    N_PROCS = 2

    def setUp(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when testing, but will hopefully be fixed down the line
            return

        super().setUp()

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        # Create the solver
        self.CFDSolver = ADFLOW(options=copy.deepcopy(options), debug=True)

        self.ap = copy.deepcopy(self.aero_prob)
        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)


    # ------------------- Derivative routine checks ----------------------------
    def test_wDot(self):
        # perturb each input and check that the outputs match the FD to with in reason 
        wDot = self.CFDSolver.getStatePerturbation(321)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True
        )
        resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode='FD', h=1e-8
        )
        
        np.testing.assert_allclose(resDot_FD, resDot, rtol=5e-4, err_msg='residual')
        
        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func],  rtol=1e-5, err_msg=func)

        np.testing.assert_allclose(fDot_FD, fDot, rtol=5e-4, err_msg='forces')
    
    def test_xVDot(self):
        # perturb each input and check that the outputs match the FD to with in reason 
        xVDot = self.CFDSolver.getSpatialPerturbation(314)

        resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True
        )
        
        resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode='FD', h=1e-8
        )
        
        idx_max = np.argmax((resDot_FD -resDot)/resDot)
        print(resDot[idx_max], resDot_FD[idx_max])
        
        np.testing.assert_allclose(resDot_FD, resDot, atol=5e-4, err_msg='residual')
        
        for func in funcsDot:
            np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func],  rtol=5e-6, err_msg=func)

        np.testing.assert_allclose(fDot_FD, fDot, rtol=5e-4, err_msg='forces')
    
    def test_xDvDot(self):
        # perturb each input and check that the outputs match the FD to with in reason 
        step_size = {
            "alpha":1e-4,
            "beta":1e-5,
            "mach":1e-5,
            "P": 1e-1,
            "T": 1e-4,
            "xRef": 1e-5,
            "yRef": 1e-5,
            "zRef": 1e-5,
            
        }
        
        for aeroDV in self.ap.DVs.values():
            key = aeroDV.key
            xDvDot = {key: 1.0}

            resDot, funcsDot, fDot = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True
            )
            
            resDot_FD, funcsDot_FD, fDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, mode='FD', h=step_size[key]
            )
            
            # the tolerances here are loose becuase different ouputs have different optimal steps
            np.testing.assert_allclose(resDot_FD, resDot, atol=5e-5, err_msg=f'residual wrt {key}')
            
            for func in funcsDot:
                if np.abs(funcsDot[func]) <= 1e-16:
                    np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func],  atol=5e-5, err_msg=f'{func} wrt {key}')
                else:   
                    np.testing.assert_allclose(funcsDot_FD[func], funcsDot[func],  rtol=1e-3, err_msg=f'{func} wrt {key}')

            np.testing.assert_allclose(fDot_FD, fDot, atol=5e-7, err_msg=f'forces wrt {key}')


@parameterized_class(test_params)
class TestJacVecFwdCS(reg_test_classes.CmplxRegTest):
    """
    Tests jacobian vector products against CS.

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
        mesh_options.update({"gridFile": options["gridfile"]})

        self.ap = copy.deepcopy(self.aero_prob)

        # Setup aeroproblem
        self.ap.evalFuncs = self.evalFuncs

        # add the default dvs to the problem
        for dv in defaultAeroDVs:
            self.ap.addDV(dv)

        self.CFDSolver = ADFLOW_C(options=options, debug=True)

        self.CFDSolver.setMesh(USMesh_C(options=mesh_options))
        self.CFDSolver.setDVGeo(setDVGeo(self.ffdFile, cmplx=True))

        # propagates the values from the restart file throughout the code
        self.CFDSolver.getResidual(self.ap)

    def cmplx_test_aero_dvs(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

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

        self.handler.root_add_dict("Eval Functions Sens:", funcsSens, rtol=1e-8, atol=5e-10)

    def cmplx_test_geom_dvs(self):
        if not hasattr(self, "name"):
            # return immediately when the setup method is being called on the based class and NOT the
            # classes created using parametrized
            # this will happen when training, but will hopefully be fixed down the line
            return

        # redo the setup for a cmplx test
        funcsSens = defaultdict(lambda: {})

        xRef = {"twist": [0.0] * 6, "span": [0.0], "shape": numpy.zeros(72, dtype="D")}

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
                ref_val = ref_val.flatten()[0]

                numpy.testing.assert_allclose(funcsSens[key][dv_key], ref_val, atol=5e-9, rtol=5e-9, err_msg=err_msg)

        if MPI.COMM_WORLD.rank == 0:
            print("====================================")
            print(self.name, funcsSens)
            print("====================================")



if __name__ == "__main__":
    unittest.main()
