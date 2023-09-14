# built-ins
import unittest
import numpy as np
import os
import copy

# need to import mpi4py for dot product tests
from mpi4py import MPI

# MACH classes
from adflow import ADFLOW
from idwarp import USMesh
from pygeo import DVGeometry
from adflow import ADFLOW_C

# import the testing utilities
import reg_test_utils as utils
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_actuator_pipe
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))


class ActuatorBasicTests(reg_test_classes.RegTest):
    """
    Tests for the actuator zone.
    """

    N_PROCS = 2
    ref_file = "actuator_tests.json"

    def setUp(self):
        super().setUp()

        self.options = {
            "gridfile": os.path.join(baseDir, "../../input_files/actuator_test_pipe.cgns"),
            # the restart file was ran with thrust = 600 N and heat = 1e5 W
            "restartfile": os.path.join(baseDir, "../../input_files/actuator_test_pipe.cgns"),
            "writevolumesolution": False,
            "writesurfacesolution": False,
            "writetecplotsurfacesolution": False,
            "mgcycle": "sg",
            "ncycles": 1000,
            "useanksolver": True,
            "usenksolver": True,
            "anksecondordswitchtol": 1e-2,
            "nkswitchtol": 1e-4,
            "volumevariables": ["temp", "mach", "resrho"],
            "surfacevariables": ["temp", "vx", "vy", "vz", "p", "ptloss", "mach", "rho"],
            "equationType": "Euler",
            "l2convergence": 1e-13,
            "adjointl2convergence": 1e-13,
        }

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        CFDSolver = ADFLOW(options=options)

        CFDSolver.addFunction("mdot", "inlet", name="mdot_in")
        CFDSolver.addFunction("mdot", "outlet", name="mdot_out")

        CFDSolver.addFunction("aavgptot", "outlet", name="aavgptot_out")
        CFDSolver.addFunction("aavgptot", "inlet", name="aavgptot_in")

        CFDSolver.addFunction("mavgttot", "outlet", name="mavgttot_out")
        CFDSolver.addFunction("mavgttot", "inlet", name="mavgttot_in")

        CFDSolver.addFunction("aavgps", "outlet", name="aavgps_out")
        CFDSolver.addFunction("aavgps", "inlet", name="aavgps_in")

        CFDSolver.addFunction("area", "inlet", name="area_in")
        CFDSolver.addFunction("area", "outlet", name="area_out")

        CFDSolver.addFunction("mavgvx", "inlet", name="mavgvx_in")
        CFDSolver.addFunction("mavgvx", "outlet", name="mavgvx_out")

        CFDSolver.addFunction("forcexpressure", "inlet", name="forcexpressure_in")
        CFDSolver.addFunction("forcexpressure", "outlet", name="forcexpressure_out")

        CFDSolver.addFunction("forcexmomentum", "inlet", name="forcexmomentum_in")
        CFDSolver.addFunction("forcexmomentum", "outlet", name="forcexmomentum_out")

        CFDSolver.addFunction("mavgvi", "inlet", name="mavgvi_in")
        CFDSolver.addFunction("mavgvi", "outlet", name="mavgvi_out")

        self.CFDSolver = CFDSolver

        # this is imported from reg_aeroproblems utility script
        self.ap = ap_actuator_pipe

        actuatorFile = os.path.join(baseDir, "../../input_files/actuator_test_disk.xyz")
        self.CFDSolver.addActuatorRegion(
            actuatorFile,
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            "actuator_region",
            # we will set these individually in the tests below
            thrust=0.0,
            torque=0.0,
            heat=0.0,
        )

        # add thrust and heat as AP DVs
        self.ap.setBCVar("Thrust", 0.0, "actuator_region")
        self.ap.addDV("Thrust", family="actuator_region", units="N", name="thrust")

        self.ap.setBCVar("Heat", 0.0, "actuator_region")
        self.ap.addDV("Heat", family="actuator_region", units="J/s", name="heat")

        # also add flowpower as an AZ function
        CFDSolver.addFunction("flowpower", "actuator_region", name="flowpower_az")

    def test_actuator_thrust(self):
        "Tests if the correct amount of momentum is added to the flow by the actuator"

        # set the az force
        az_force = 600.0
        # need to set all dvs because training may re-use leftover dvs from a previous test
        self.ap.setDesignVars({"thrust": az_force, "heat": 0.0})

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        # negate mdot out because of the normal, mdot in is already positive
        mdot_i = funcs[self.ap.name + "_mdot_in"]
        mdot_o = -funcs[self.ap.name + "_mdot_out"]

        vx_i = funcs[self.ap.name + "_mavgvx_in"]
        vx_o = funcs[self.ap.name + "_mavgvx_out"]

        area_i = funcs[self.ap.name + "_area_in"]
        area_o = funcs[self.ap.name + "_area_out"]

        aavgps_i = funcs[self.ap.name + "_aavgps_in"]
        aavgps_o = funcs[self.ap.name + "_aavgps_out"]

        ttot_i = funcs[self.ap.name + "_mavgttot_in"]
        ttot_o = funcs[self.ap.name + "_mavgttot_out"]

        # also get the pressure and momentum forces directly from CFD
        fp_i = funcs[self.ap.name + "_forcexpressure_in"]
        fp_o = funcs[self.ap.name + "_forcexpressure_out"]

        fm_i = funcs[self.ap.name + "_forcexmomentum_in"]
        fm_o = funcs[self.ap.name + "_forcexmomentum_out"]

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # this is the analytical force based on primitive values (like mdot, ps etc)
        my_force = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)

        # this is the force computed by the momentum integration directly from CFD
        # just sum these up, the forces contain the correct normals from CFD
        cfd_force = fp_o + fm_o + fp_i + fm_i

        # The low accuracy is because the integrated quantities don't have a lot of precision
        np.testing.assert_allclose(my_force, az_force, rtol=1e-3)
        np.testing.assert_allclose(cfd_force, az_force, rtol=1e-3)

        ##################
        # TEST POWER ADDED
        ##################

        # this is the integration done in the AZ
        az_power = funcs[self.ap.name + "_flowpower_az"]

        # this is the energy balance of the control volume.
        # Cp of air is taken as 1004.5 J/kg
        my_power = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)

        # the tolerance is slightly worse but not terrible
        np.testing.assert_allclose(my_power, az_power, rtol=1.5e-3)

    def test_actuator_heat(self):
        "Tests if the correct amount of heat is added to the flow by the actuator"

        # set the az heat
        az_heat = 1e5
        # need to set all dvs because training may re-use leftover dvs from a previous test
        self.ap.setDesignVars({"thrust": 0.0, "heat": az_heat})

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        # negate mdot out because of the normal, mdot in is already positive
        mdot_i = funcs[self.ap.name + "_mdot_in"]
        mdot_o = -funcs[self.ap.name + "_mdot_out"]

        vx_i = funcs[self.ap.name + "_mavgvx_in"]
        vx_o = funcs[self.ap.name + "_mavgvx_out"]

        area_i = funcs[self.ap.name + "_area_in"]
        area_o = funcs[self.ap.name + "_area_out"]

        aavgps_i = funcs[self.ap.name + "_aavgps_in"]
        aavgps_o = funcs[self.ap.name + "_aavgps_out"]

        ttot_i = funcs[self.ap.name + "_mavgttot_in"]
        ttot_o = funcs[self.ap.name + "_mavgttot_out"]

        # also get the pressure and momentum forces directly from CFD
        fp_i = funcs[self.ap.name + "_forcexpressure_in"]
        fp_o = funcs[self.ap.name + "_forcexpressure_out"]

        fm_i = funcs[self.ap.name + "_forcexmomentum_in"]
        fm_o = funcs[self.ap.name + "_forcexmomentum_out"]

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # this is the analytical force based on primitive values (like mdot, ps etc)
        my_force = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)

        # this is the force computed by the momentum integration directly from CFD
        # just sum these up, the forces contain the correct normals from CFD
        cfd_force = fp_o + fm_o + fp_i + fm_i

        # the true value is 0. However, we get an error around 0.07964251
        # because of accumulated integration/convergence/precision errors.
        np.testing.assert_allclose(my_force, 0.0, atol=0.08)
        np.testing.assert_allclose(cfd_force, 0.0, atol=0.08)

        ##################
        # TEST POWER ADDED
        ##################

        # this is the total energy rise based on total temperature
        cfd_heat = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)

        # the tolerance is slightly worse but not terrible
        np.testing.assert_allclose(cfd_heat, az_heat, rtol=1.5e-3)

    def test_actuator_thrust_and_heat(self):
        "Tests if the correct amount of momentum and heat is added to the flow by the actuator"

        # set the az force
        az_force = 600.0
        az_heat = 1e5
        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        # negate mdot out because of the normal, mdot in is already positive
        mdot_i = funcs[self.ap.name + "_mdot_in"]
        mdot_o = -funcs[self.ap.name + "_mdot_out"]

        vx_i = funcs[self.ap.name + "_mavgvx_in"]
        vx_o = funcs[self.ap.name + "_mavgvx_out"]

        area_i = funcs[self.ap.name + "_area_in"]
        area_o = funcs[self.ap.name + "_area_out"]

        aavgps_i = funcs[self.ap.name + "_aavgps_in"]
        aavgps_o = funcs[self.ap.name + "_aavgps_out"]

        ttot_i = funcs[self.ap.name + "_mavgttot_in"]
        ttot_o = funcs[self.ap.name + "_mavgttot_out"]

        # also get the pressure and momentum forces directly from CFD
        fp_i = funcs[self.ap.name + "_forcexpressure_in"]
        fp_o = funcs[self.ap.name + "_forcexpressure_out"]

        fm_i = funcs[self.ap.name + "_forcexmomentum_in"]
        fm_o = funcs[self.ap.name + "_forcexmomentum_out"]

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # this is the analytical force based on primitive values (like mdot, ps etc)
        my_force = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)

        # this is the force computed by the momentum integration directly from CFD
        # just sum these up, the forces contain the correct normals from CFD
        cfd_force = fp_o + fm_o + fp_i + fm_i

        # The low accuracy is because the integrated quantities don't have a lot of precision
        np.testing.assert_allclose(my_force, az_force, rtol=1e-3)
        np.testing.assert_allclose(cfd_force, az_force, rtol=1e-3)

        ##################
        # TEST POWER ADDED
        ##################

        # this is the integration done in the AZ plus the heat added by us
        az_power = funcs[self.ap.name + "_flowpower_az"] + az_heat

        # this is the energy balance of the control volume.
        # Cp of air is taken as 1004.5 J/kg
        my_power = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)

        # the tolerance is slightly worse but not terrible
        np.testing.assert_allclose(my_power, az_power, rtol=1.5e-3)

        #################
        # TEST STATE NORM
        #################
        utils.assert_states_allclose(self.handler, self.CFDSolver)

        ################
        # TEST ALL FUNCS
        ################
        self.handler.root_print("all functionals")
        self.handler.root_add_dict("all functionals", funcs, rtol=1e-12, atol=1e-12)

    def test_actuator_adjoint(self):
        "Tests if the adjoint sensitivities are correct for the AZ DVs"

        # define user functions
        my_force_functions = [
            "mdot_in",
            "mdot_out",
            "mavgvx_in",
            "mavgvx_out",
            "area_in",
            "area_out",
            "aavgps_in",
            "aavgps_out",
        ]
        cfd_force_functions = [
            "forcexpressure_in",
            "forcexpressure_out",
            "forcexmomentum_in",
            "forcexmomentum_out",
        ]
        my_power_functions = [
            "mdot_in",
            "mdot_out",
            "mavgttot_in",
            "mavgttot_out",
        ]

        def f_my_force(funcs):
            mdot_i = funcs["mdot_in"]
            mdot_o = -funcs["mdot_out"]

            vx_i = funcs["mavgvx_in"]
            vx_o = funcs["mavgvx_out"]

            area_i = funcs["area_in"]
            area_o = funcs["area_out"]

            aavgps_i = funcs["aavgps_in"]
            aavgps_o = funcs["aavgps_out"]

            funcs["my_force"] = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)
            return funcs

        def f_cfd_force(funcs):
            fp_i = funcs["forcexpressure_in"]
            fp_o = funcs["forcexpressure_out"]

            fm_i = funcs["forcexmomentum_in"]
            fm_o = funcs["forcexmomentum_out"]

            funcs["cfd_force"] = fp_o + fm_o + fp_i + fm_i
            return funcs

        def f_my_power(funcs):
            mdot_i = funcs["mdot_in"]
            mdot_o = -funcs["mdot_out"]

            ttot_i = funcs["mavgttot_in"]
            ttot_o = funcs["mavgttot_out"]

            funcs["my_power"] = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)
            return funcs

        self.CFDSolver.addUserFunction("my_force", my_force_functions, f_my_force)
        self.CFDSolver.addUserFunction("cfd_force", cfd_force_functions, f_cfd_force)
        self.CFDSolver.addUserFunction("my_power", my_power_functions, f_my_power)

        # set the az force
        az_force = 0.0
        az_heat = 1e5
        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        funcs = {}
        funcsSens = {}
        # self.CFDSolver.evalFunctions(self.ap, funcs)
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["my_force", "cfd_force", "my_power", "mavgvi_out"])
        self.CFDSolver.evalFunctionsSens(
            self.ap, funcsSens, evalFuncs=["my_force", "cfd_force", "my_power", "mavgvi_out"]
        )
        # check if adjoint failed
        self.assert_adjoint_failure()

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # we test these w.r.t. thrust and heat analytically
        # The low accuracy is because the integrated quantities don't have a lot of precision

        np.testing.assert_allclose(funcsSens["actuator_pipe_my_force"]["thrust"], 1, rtol=1e-3)
        np.testing.assert_allclose(funcsSens["actuator_pipe_cfd_force"]["thrust"], 1, rtol=1e-3)

        # heat addition should not affect these
        np.testing.assert_allclose(funcsSens["actuator_pipe_my_force"]["heat"] + 100, 100, rtol=1e-3)
        np.testing.assert_allclose(funcsSens["actuator_pipe_cfd_force"]["heat"] + 100, 100, rtol=1e-3)

        # also test the actual values from the ref file
        self.handler.root_print("my_force sens")
        self.handler.root_add_dict("my_force sens", funcsSens["actuator_pipe_my_force"], rtol=1e-12, atol=1e-12)

        self.handler.root_print("cfd_force sens")
        self.handler.root_add_dict("cfd_force sens", funcsSens["actuator_pipe_cfd_force"], rtol=1e-12, atol=1e-12)

        ##################
        # TEST POWER ADDED
        ##################

        # analytically test heat addition
        # this should be equal to one because we are not adding any thrust in this test.
        # if we also added thrust, addition of heat would affect flow power integration
        # due to the changes in the flowfield and as a result the derivative would not be one.
        # flowpower is more complicated here so we just check with json reference
        np.testing.assert_allclose(funcsSens["actuator_pipe_my_power"]["heat"], 1, rtol=1e-3)

        # test values with the ref file
        self.handler.root_print("my_power sens")
        self.handler.root_add_dict("my_power sens", funcsSens["actuator_pipe_my_power"], rtol=1e-12, atol=1e-12)

        # test mavgvi
        self.handler.root_print("mavgvi sens")
        self.handler.root_add_dict("mavgvi sens", funcsSens["actuator_pipe_mavgvi_out"], rtol=1e-12, atol=1e-12)

    def test_actuator_flowpower_adjoint(self):
        "we test this adjoint separately because we need to have a finite thrust for this to actually test"

        # set the az force
        az_force = 600.0
        az_heat = 1e5
        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        funcs = {}
        funcsSens = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)
        self.CFDSolver.evalFunctionsSens(self.ap, funcsSens, evalFuncs=["flowpower_az"])

        #############################
        # TEST FLOW POWER INTEGRATION
        #############################

        self.handler.root_print("flowpower sens")
        self.handler.root_add_dict("flowpower sens", funcsSens["actuator_pipe_flowpower_az"], rtol=1e-12, atol=1e-12)

    def test_actuator_partials(self):
        az_force = 600.0
        az_heat = 1e5
        # need to set all dvs because training may re-use leftover dvs from a previous test
        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        #############
        # TEST FWD AD
        #############

        # save these for the dot product tests
        resDot = {}

        for DV in self.ap.DVs.values():
            # regular aero DVs do not need the family name, but the BC names do require the fam names
            xDvDot = {DV.key + "_actuator_region": 1}

            resDot[DV.key] = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot,
                # actuator changes only affect residuals
                residualDeriv=True,
            )

            self.handler.root_print("||dR/d%s||" % DV.key)
            self.handler.par_add_norm("||dR/d%s||" % DV.key, resDot[DV.key], rtol=1e-12, atol=1e-12)

        #############
        # TEST BWD AD
        #############

        # we are only interested in checking the xDvBar result here
        dwBar = self.CFDSolver.getStatePerturbation(123)
        xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, xDvDerivAero=True)

        for DV in self.ap.DVs.values():
            self.handler.root_print("||dR/d%s^T||" % DV.key)
            self.handler.root_add_val("||dR/d%s^T||" % DV.key, xDvBar[DV.key.lower()], rtol=1e-12, atol=1e-12)

        ###############
        # TEST DOT PROD
        ###############

        for DV in self.ap.DVs.values():
            # this product is already the same on all procs because adflow handles comm for xDvBar
            first = np.dot(xDvBar[DV.key.lower()], 1.0)
            # we need to reduce the resDot dot dwBar ourselves
            secondLocal = np.dot(resDot[DV.key], dwBar)
            second = self.CFDSolver.comm.allreduce(secondLocal, op=MPI.SUM)

            # compare the final products
            np.testing.assert_array_almost_equal(first, second, decimal=14)

    def test_actuator_partials_warp(self):
        az_force = 600.0
        az_heat = 1e5

        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        ffd_file = os.path.join(baseDir, "../../input_files/actuator_test_ffd.xyz")
        DVGeo = DVGeometry(ffd_file)
        DVGeo.addLocalDV("xdir", lower=-0.5, upper=0.5, axis="x", scale=1.0)
        DVGeo.addLocalDV("ydir", lower=-0.5, upper=0.5, axis="y", scale=1.0)
        DVGeo.addLocalDV("zdir", lower=-0.5, upper=0.5, axis="z", scale=1.0)
        self.CFDSolver.setDVGeo(DVGeo)

        dv_keys = ["xdir", "ydir", "zdir"]

        grid_file = os.path.join(baseDir, "../../input_files/actuator_test_pipe.cgns")
        mesh = USMesh({"gridFile": grid_file})
        self.CFDSolver.setMesh(mesh)

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        #############
        # TEST FWD AD
        #############

        # save these for the dot product tests
        resDot = {}

        for key in dv_keys:
            xDvDot = {key: 1}

            resDot[key] = self.CFDSolver.computeJacobianVectorProductFwd(
                xDvDot=xDvDot,
                # actuator changes only affect residuals
                residualDeriv=True,
            )

            self.handler.root_print("||dR/d%s||" % key)
            self.handler.par_add_norm("||dR/d%s||" % key, resDot[key], rtol=1e-12, atol=1e-12)

        #############
        # TEST BWD AD
        #############

        # we are only interested in checking the xDvBar result here
        dwBar = self.CFDSolver.getStatePerturbation(123)
        xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, xDvDeriv=True)

        for key in dv_keys:
            self.handler.root_print("||dR/d%s^T||" % key)
            self.handler.root_add_val("||dR/d%s^T||" % key, xDvBar[key.lower()], rtol=1e-12, atol=1e-12)

        ###############
        # TEST DOT PROD
        ###############

        for key in dv_keys:
            # this product is already the same on all procs because adflow handles comm for xDvBar
            xDvBar_tmp = xDvBar[key.lower()]
            first = np.dot(xDvBar_tmp, np.ones(xDvBar_tmp.shape[1]))
            # we need to reduce the resDot dot dwBar ourselves
            secondLocal = np.dot(resDot[key], dwBar)
            second = self.CFDSolver.comm.allreduce(secondLocal, op=MPI.SUM)

            # compare the final products
            np.testing.assert_array_almost_equal(first, second, decimal=11)


class ActuatorCmplxTests(reg_test_classes.CmplxRegTest):
    """
    Complex step tests for the actuator zone.
    """

    N_PROCS = 2
    ref_file = "actuator_tests.json"
    h = 1e-40

    def setUp(self):
        super().setUp()

        self.options = {
            "gridfile": os.path.join(baseDir, "../../input_files/actuator_test_pipe.cgns"),
            # the restart file was ran with thrust = 600 N and heat = 1e5 W
            # "restartfile": os.path.join(baseDir, "../../input_files/actuator_test_pipe.cgns"),
            "writevolumesolution": False,
            "writesurfacesolution": False,
            "writetecplotsurfacesolution": False,
            "mgcycle": "sg",
            "ncycles": 1000,
            "useanksolver": True,
            "usenksolver": True,
            "ankswitchtol": 10.0,
            "anksecondordswitchtol": 1e-2,
            "ankinnerpreconits": 1,
            "monitorvariables": ["cpu", "resrho", "resrhoe"],
            "volumevariables": ["temp", "mach", "resrho"],
            "surfacevariables": ["temp", "vx", "vy", "vz", "p", "ptloss", "mach", "rho"],
            "equationType": "Euler",
            "l2convergence": 1e-13,
            "adjointl2convergence": 1e-13,
        }

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        # CFDSolver = ADFLOW(options=options)
        CFDSolver = ADFLOW_C(options=options, debug=True)

        CFDSolver.addFunction("mdot", "inlet", name="mdot_in")
        CFDSolver.addFunction("mdot", "outlet", name="mdot_out")

        CFDSolver.addFunction("aavgptot", "outlet", name="aavgptot_out")
        CFDSolver.addFunction("aavgptot", "inlet", name="aavgptot_in")

        CFDSolver.addFunction("mavgttot", "outlet", name="mavgttot_out")
        CFDSolver.addFunction("mavgttot", "inlet", name="mavgttot_in")

        CFDSolver.addFunction("aavgps", "outlet", name="aavgps_out")
        CFDSolver.addFunction("aavgps", "inlet", name="aavgps_in")

        CFDSolver.addFunction("area", "inlet", name="area_in")
        CFDSolver.addFunction("area", "outlet", name="area_out")

        CFDSolver.addFunction("mavgvx", "inlet", name="mavgvx_in")
        CFDSolver.addFunction("mavgvx", "outlet", name="mavgvx_out")

        # we just test mavgvi on the outlet plane for derivatives
        CFDSolver.addFunction("mavgvi", "outlet", name="mavgvi_out")

        CFDSolver.addFunction("forcexpressure", "inlet", name="forcexpressure_in")
        CFDSolver.addFunction("forcexpressure", "outlet", name="forcexpressure_out")

        CFDSolver.addFunction("forcexmomentum", "inlet", name="forcexmomentum_in")
        CFDSolver.addFunction("forcexmomentum", "outlet", name="forcexmomentum_out")

        self.CFDSolver = CFDSolver

        # this is imported from reg_aeroproblems utility script
        self.ap = ap_actuator_pipe

        actuatorFile = os.path.join(baseDir, "../../input_files/actuator_test_disk.xyz")
        self.CFDSolver.addActuatorRegion(
            actuatorFile,
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            "actuator_region",
            # we will set these individually in the tests below
            thrust=0.0,
            torque=0.0,
            heat=0.0,
        )

        # add thrust and heat as AP DVs
        self.ap.setBCVar("Thrust", 0.0, "actuator_region")
        self.ap.addDV("Thrust", family="actuator_region", units="N", name="thrust")

        self.ap.setBCVar("Heat", 0.0, "actuator_region")
        self.ap.addDV("Heat", family="actuator_region", units="J/s", name="heat")

        # also add flowpower as an AZ function
        CFDSolver.addFunction("flowpower", "actuator_region", name="flowpower_az")

    def cmplx_test_actuator_adjoint(self):
        "Tests if the adjoint sensitivities are correct for the AZ DVs"

        # define user functions
        my_force_functions = [
            "mdot_in",
            "mdot_out",
            "mavgvx_in",
            "mavgvx_out",
            "area_in",
            "area_out",
            "aavgps_in",
            "aavgps_out",
        ]
        cfd_force_functions = [
            "forcexpressure_in",
            "forcexpressure_out",
            "forcexmomentum_in",
            "forcexmomentum_out",
        ]
        my_power_functions = [
            "mdot_in",
            "mdot_out",
            "mavgttot_in",
            "mavgttot_out",
        ]

        def f_my_force(funcs):
            mdot_i = funcs["mdot_in"]
            mdot_o = -funcs["mdot_out"]

            vx_i = funcs["mavgvx_in"]
            vx_o = funcs["mavgvx_out"]

            area_i = funcs["area_in"]
            area_o = funcs["area_out"]

            aavgps_i = funcs["aavgps_in"]
            aavgps_o = funcs["aavgps_out"]

            funcs["my_force"] = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)
            return funcs

        def f_cfd_force(funcs):
            fp_i = funcs["forcexpressure_in"]
            fp_o = funcs["forcexpressure_out"]

            fm_i = funcs["forcexmomentum_in"]
            fm_o = funcs["forcexmomentum_out"]

            funcs["cfd_force"] = fp_o + fm_o + fp_i + fm_i
            return funcs

        def f_my_power(funcs):
            mdot_i = funcs["mdot_in"]
            mdot_o = -funcs["mdot_out"]

            ttot_i = funcs["mavgttot_in"]
            ttot_o = funcs["mavgttot_out"]

            funcs["my_power"] = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)
            return funcs

        self.CFDSolver.addUserFunction("my_force", my_force_functions, f_my_force)
        self.CFDSolver.addUserFunction("cfd_force", cfd_force_functions, f_cfd_force)
        self.CFDSolver.addUserFunction("my_power", my_power_functions, f_my_power)

        # set the az force
        az_force = 0.0
        az_heat = 1e5
        aDV = {"thrust": az_force, "heat": az_heat}
        self.ap.setDesignVars(aDV)

        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        funcs = {}
        funcsSensCS = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["my_force", "cfd_force", "my_power", "mavgvi_out"])

        funcs_plus = {}
        for dv in ["thrust", "heat"]:
            # save the old dv
            dvsave = aDV[dv]

            # perturb
            aDV[dv] += self.h * 1j
            self.ap.setDesignVars(aDV)

            # call solver again
            # we can also not re-set the flow and call the solver a few times until complex parts converge
            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver(self.ap)

            # check if solution failed
            self.assert_solution_failure()

            # save the new funcs in a dict
            funcs_plus[dv] = {}

            # eval functions
            self.CFDSolver.evalFunctions(
                self.ap, funcs_plus[dv], evalFuncs=["my_force", "cfd_force", "my_power", "mavgvi_out"]
            )

            # compute the sens
            funcsSensCS[dv] = {}
            for f in ["my_force", "cfd_force", "my_power", "mavgvi_out"]:
                fname = "actuator_pipe_" + f
                funcsSensCS[dv][f] = np.imag(funcs_plus[dv][fname]) / self.h

            # reset the DV
            aDV[dv] = dvsave

        ##################
        # TEST DERIVATIVES
        ##################

        # here we compare the CS derivatives to the adjoint values computed by the real test
        # we treat the CS value as the truth, so if this test passes,
        # we assume the adjoint sensitivities are also true

        # my force

        # thrust
        ref_val = self.handler.db["my_force sens"]["thrust"]
        np.testing.assert_allclose(funcsSensCS["thrust"]["my_force"], ref_val, atol=1e-10, rtol=1e-10)

        # heat
        ref_val = self.handler.db["my_force sens"]["heat"]
        np.testing.assert_allclose(funcsSensCS["heat"]["my_force"], ref_val, atol=1e-10, rtol=1e-10)

        # cfd force

        # thrust
        ref_val = self.handler.db["cfd_force sens"]["thrust"]
        np.testing.assert_allclose(funcsSensCS["thrust"]["cfd_force"], ref_val, atol=1e-10, rtol=1e-10)

        # heat
        ref_val = self.handler.db["cfd_force sens"]["heat"]
        np.testing.assert_allclose(funcsSensCS["heat"]["cfd_force"], ref_val, atol=1e-10, rtol=1e-10)

        # my power

        # thrust
        ref_val = self.handler.db["my_power sens"]["thrust"]
        np.testing.assert_allclose(funcsSensCS["thrust"]["my_power"], ref_val, atol=1e-10, rtol=1e-10)

        # heat
        ref_val = self.handler.db["my_power sens"]["heat"]
        np.testing.assert_allclose(funcsSensCS["heat"]["my_power"], ref_val, atol=1e-10, rtol=1e-10)

        # mavgvi
        # thrust
        ref_val = self.handler.db["mavgvi sens"]["thrust"]
        np.testing.assert_allclose(funcsSensCS["thrust"]["mavgvi_out"], ref_val, atol=1e-10, rtol=1e-10)

        # heat
        ref_val = self.handler.db["mavgvi sens"]["heat"]
        np.testing.assert_allclose(funcsSensCS["heat"]["mavgvi_out"], ref_val, atol=1e-10, rtol=1e-10)

    def cmplx_test_actuator_flowpower_adjoint(self):
        "we test this adjoint separately because we need to have a finite thrust for this to actually test"

        # set the az force
        az_force = 600.0
        az_heat = 1e5
        aDV = {"thrust": az_force, "heat": az_heat}
        self.ap.setDesignVars(aDV)

        # We dont need to rerun, restart file has the correct state. just run a residual
        self.CFDSolver(self.ap)

        # check if solution failed
        self.assert_solution_failure()

        funcs = {}
        funcsSensCS = {}
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=["flowpower_az"])

        funcs_plus = {}
        for dv in ["thrust", "heat"]:
            # save the old dv
            dvsave = aDV[dv]

            # perturb
            aDV[dv] += self.h * 1j
            self.ap.setDesignVars(aDV)

            # call solver again
            # we can also not re-set the flow and call the solver a few times until complex parts converge
            self.CFDSolver.resetFlow(self.ap)
            self.CFDSolver(self.ap)

            # check if solution failed
            self.assert_solution_failure()

            # save the new funcs in a dict
            funcs_plus[dv] = {}

            # eval functions
            self.CFDSolver.evalFunctions(self.ap, funcs_plus[dv], evalFuncs=["flowpower_az"])

            # compute the sens
            funcsSensCS[dv] = np.imag(funcs_plus[dv]["actuator_pipe_flowpower_az"]) / self.h

            # reset the DV
            aDV[dv] = dvsave

        ##################
        # TEST DERIVATIVES
        ##################

        # here we compare the CS derivatives to the adjoint values computed by the real test
        # we treat the CS value as the truth, so if this test passes,
        # we assume the adjoint sensitivities are also true

        # thrust
        ref_val = self.handler.db["flowpower sens"]["thrust"]
        np.testing.assert_allclose(funcsSensCS["thrust"], ref_val, atol=1e-10, rtol=1e-10)

        # heat
        ref_val = self.handler.db["flowpower sens"]["heat"]
        np.testing.assert_allclose(funcsSensCS["heat"], ref_val, atol=1e-10, rtol=1e-10)


if __name__ == "__main__":
    unittest.main()
