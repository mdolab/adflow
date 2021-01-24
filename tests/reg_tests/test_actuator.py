# built-ins
import unittest
import numpy as np
import os
import copy
from pprint import pprint as pp

# need to import mpi4py for dot product tests
from mpi4py import MPI

# MACH classes
from adflow import ADFLOW

# import the testing utilities that live a few directories up
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

    # TODO tests:
    # CS test; test flowpower and user functions compare with adjoint

    def setUp(self):

        super().setUp()

        self.options = {'gridfile': os.path.join(baseDir, '../../inputFiles/actuator_test_pipe.cgns'),
                        # the restart file was ran with thrust = 600 N and heat = 1e5 W
                        'restartfile': os.path.join(baseDir, '../../inputFiles/actuator_test_pipe.cgns'),
                        "writevolumesolution": False,
                        "writesurfacesolution": False,
                        "writetecplotsurfacesolution": False,
                        'mgcycle':'sg',
                        'ncycles':1000,
                        'useanksolver':True,
                        'usenksolver':True,
                        'ankswitchtol':1.0,
                        'anksubspacesize': -1,
                        'anklinearsolvetol': 0.05,
                        'ankjacobianlag': 10,
                        'ankinnerpreconits':1,
                        'volumevariables': ['temp', 'mach', 'resrho' ],
                        'surfacevariables':['temp', 'vx', 'vy', 'vz', 'p', 'ptloss', 'mach', 'rho'],
                        'equationType':'Euler',
                        'l2convergence': 1e-13,
                        'adjointl2convergence': 1e-13,
                        }

        options = copy.copy(adflowDefOpts)
        options["outputdirectory"] = os.path.join(baseDir, options["outputdirectory"])
        options.update(self.options)

        CFDSolver = ADFLOW(options=options)

        CFDSolver.addFunction('mdot', 'inlet', name="mdot_in")
        CFDSolver.addFunction('mdot', 'outlet', name="mdot_out")

        CFDSolver.addFunction('aavgptot', 'outlet', name="aavgptot_out")
        CFDSolver.addFunction('aavgptot', 'inlet', name="aavgptot_in")

        CFDSolver.addFunction('mavgttot', 'outlet', name="mavgttot_out")
        CFDSolver.addFunction('mavgttot', 'inlet', name="mavgttot_in")

        CFDSolver.addFunction('aavgps', 'outlet', name="aavgps_out")
        CFDSolver.addFunction('aavgps', 'inlet', name="aavgps_in")

        CFDSolver.addFunction('area', 'inlet', name="area_in")
        CFDSolver.addFunction('area', 'outlet', name="area_out")

        CFDSolver.addFunction('mavgvx', 'inlet', name="mavgvx_in")
        CFDSolver.addFunction('mavgvx', 'outlet', name="mavgvx_out")

        CFDSolver.addFunction('forcexpressure', 'inlet', name="forcexpressure_in")
        CFDSolver.addFunction('forcexpressure', 'outlet', name="forcexpressure_out")

        CFDSolver.addFunction('forcexmomentum', 'inlet', name="forcexmomentum_in")
        CFDSolver.addFunction('forcexmomentum', 'outlet', name="forcexmomentum_out")

        self.CFDSolver = CFDSolver

        # this is imported from reg_aeroproblems utility script
        self.ap = ap_actuator_pipe

        actuatorFile = os.path.join(baseDir, '../../inputFiles/actuator_test_disk.xyz')
        self.CFDSolver.addActuatorRegion(
            actuatorFile,
            np.array([0, 0, 0]),
            np.array([1, 0, 0]),
            'actuator_region',
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

    def test_actuator_momentum_flux(self):
        "Tests if the correct amount of momentum is added to the flow by the actuator"

        # set the az force
        az_force = 600.
        # need to set all dvs because training may re-use leftover dvs from a previous test
        self.ap.setDesignVars({"thrust": az_force, "heat": 0.0})

        # we have to rerun because restart file was ran with a different heat value
        self.CFDSolver(self.ap)
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        # negate mdot out because of the normal, mdot in is already positive
        mdot_i = funcs[self.ap.name + '_mdot_in']
        mdot_o = -funcs[self.ap.name + '_mdot_out']

        vx_i = funcs[self.ap.name + '_mavgvx_in']
        vx_o = funcs[self.ap.name + '_mavgvx_out']

        area_i = funcs[self.ap.name +'_area_in']
        area_o = funcs[self.ap.name +'_area_out']

        aavgps_i = funcs[self.ap.name + '_aavgps_in']
        aavgps_o = funcs[self.ap.name + '_aavgps_out']

        ttot_i = funcs[self.ap.name + '_mavgttot_in']
        ttot_o = funcs[self.ap.name + '_mavgttot_out']

        # also get the pressure and momentum forces directly from CFD
        fp_i = funcs[self.ap.name + '_forcexpressure_in']
        fp_o = funcs[self.ap.name + '_forcexpressure_out']

        fm_i = funcs[self.ap.name + '_forcexmomentum_in']
        fm_o = funcs[self.ap.name + '_forcexmomentum_out']

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # this is the analytical force based on primitive values (like mdot, ps etc)
        my_force = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)

        # this is the force computed by the momentum integration directly from CFD
        # just sum these up, the forces contain the correct normals from CFD
        cfd_force = fp_o + fm_o + fp_i + fm_i

        # The low accuracy is because the intgrated quantities don't have a lot of precision
        np.testing.assert_allclose(my_force, az_force, rtol=1e-3)
        np.testing.assert_allclose(cfd_force, az_force, rtol=1e-3)

        ##################
        # TEST POWER ADDED
        ##################

        # this is the integration done in the AZ
        az_power = funcs[self.ap.name + '_flowpower_az']

        # this is the energy balance of the control volume.
        # Cp of air is taken as 1004.5 J/kg
        my_power = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)

        # the tolerance is slightly worse but not terrible
        np.testing.assert_allclose(my_power, az_power, rtol=1.5e-3)

    def test_actuator_heat_flux(self):
        "Tests if the correct amount of heat is added to the flow by the actuator"

        # set the az heat
        az_heat = 1e5
        # need to set all dvs because training may re-use leftover dvs from a previous test
        self.ap.setDesignVars({"thrust": 0.0, "heat": az_heat})

        # we have to rerun because restart file was ran with a different thrust value
        self.CFDSolver(self.ap)
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        # negate mdot out because of the normal, mdot in is already positive
        mdot_i = funcs[self.ap.name + '_mdot_in']
        mdot_o = -funcs[self.ap.name + '_mdot_out']

        vx_i = funcs[self.ap.name + '_mavgvx_in']
        vx_o = funcs[self.ap.name + '_mavgvx_out']

        area_i = funcs[self.ap.name +'_area_in']
        area_o = funcs[self.ap.name +'_area_out']

        aavgps_i = funcs[self.ap.name + '_aavgps_in']
        aavgps_o = funcs[self.ap.name + '_aavgps_out']

        ttot_i = funcs[self.ap.name + '_mavgttot_in']
        ttot_o = funcs[self.ap.name + '_mavgttot_out']

        # also get the pressure and momentum forces directly from CFD
        fp_i = funcs[self.ap.name + '_forcexpressure_in']
        fp_o = funcs[self.ap.name + '_forcexpressure_out']

        fm_i = funcs[self.ap.name + '_forcexmomentum_in']
        fm_o = funcs[self.ap.name + '_forcexmomentum_out']

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # this is the analytical force based on primitive values (like mdot, ps etc)
        my_force = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)

        # this is the force computed by the momentum integration directly from CFD
        # just sum these up, the forces contain the correct normals from CFD
        cfd_force = fp_o + fm_o + fp_i + fm_i

        # the true value is 0, so we add 100 on both sides to have a reasonable comparison
        np.testing.assert_allclose(100 + my_force, 100.0, rtol=1e-3)
        np.testing.assert_allclose(100 + cfd_force, 100.0, rtol=1e-3)

        ##################
        # TEST POWER ADDED
        ##################

        # this is the total energy rise based on total temperature
        cfd_heat = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)

        # the tolerance is slightly worse but not terrible
        np.testing.assert_allclose(cfd_heat, az_heat, rtol=1.5e-3)

    def test_actuator_both_fluxes(self):
        "Tests if the correct amount of momentum and heat is added to the flow by the actuator"

        # set the az force
        az_force = 600.
        az_heat = 1e5
        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        # We dont need to rerun, restart file has the correct state. just run a residual
        # self.CFDSolver(self.ap)
        self.CFDSolver.getResidual(self.ap)

        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        # negate mdot out because of the normal, mdot in is already positive
        mdot_i = funcs[self.ap.name + '_mdot_in']
        mdot_o = -funcs[self.ap.name + '_mdot_out']

        vx_i = funcs[self.ap.name + '_mavgvx_in']
        vx_o = funcs[self.ap.name + '_mavgvx_out']

        area_i = funcs[self.ap.name +'_area_in']
        area_o = funcs[self.ap.name +'_area_out']

        aavgps_i = funcs[self.ap.name + '_aavgps_in']
        aavgps_o = funcs[self.ap.name + '_aavgps_out']

        ttot_i = funcs[self.ap.name + '_mavgttot_in']
        ttot_o = funcs[self.ap.name + '_mavgttot_out']

        # also get the pressure and momentum forces directly from CFD
        fp_i = funcs[self.ap.name + '_forcexpressure_in']
        fp_o = funcs[self.ap.name + '_forcexpressure_out']

        fm_i = funcs[self.ap.name + '_forcexmomentum_in']
        fm_o = funcs[self.ap.name + '_forcexmomentum_out']

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # this is the analytical force based on primitive values (like mdot, ps etc)
        my_force = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)

        # this is the force computed by the momentum integration directly from CFD
        # just sum these up, the forces contain the correct normals from CFD
        cfd_force = fp_o + fm_o + fp_i + fm_i

        # The low accuracy is because the intgrated quantities don't have a lot of precision
        np.testing.assert_allclose(my_force, az_force, rtol=1e-3)
        np.testing.assert_allclose(cfd_force, az_force, rtol=1e-3)

        ##################
        # TEST POWER ADDED
        ##################

        # this is the integration done in the AZ plus the heat added by us
        az_power = funcs[self.ap.name + '_flowpower_az'] + az_heat

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
            'mdot_in',
            'mdot_out',
            'mavgvx_in',
            'mavgvx_out',
            'area_in',
            'area_out',
            'aavgps_in',
            'aavgps_out',
        ]
        cfd_force_functions = [
            'forcexpressure_in',
            'forcexpressure_out',
            'forcexmomentum_in',
            'forcexmomentum_out',
        ]
        my_power_functions = [
            'mdot_in',
            'mdot_out',
            'mavgttot_in',
            'mavgttot_out',
        ]

        def f_my_force(funcs):
            mdot_i = funcs['mdot_in']
            mdot_o = -funcs['mdot_out']

            vx_i = funcs['mavgvx_in']
            vx_o = funcs['mavgvx_out']

            area_i = funcs['area_in']
            area_o = funcs['area_out']

            aavgps_i = funcs['aavgps_in']
            aavgps_o = funcs['aavgps_out']

            funcs['my_force'] = mdot_o * vx_o + aavgps_o * area_o - (mdot_i * vx_i + aavgps_i * area_i)
            return funcs

        def f_cfd_force(funcs):
            fp_i = funcs['forcexpressure_in']
            fp_o = funcs['forcexpressure_out']

            fm_i = funcs['forcexmomentum_in']
            fm_o = funcs['forcexmomentum_out']

            funcs['cfd_force'] = fp_o + fm_o + fp_i + fm_i
            return funcs

        def f_my_power(funcs):
            mdot_i = funcs['mdot_in']
            mdot_o = -funcs['mdot_out']

            ttot_i = funcs['mavgttot_in']
            ttot_o = funcs['mavgttot_out']

            funcs['my_power'] = 1004.5 * (mdot_o * ttot_o - mdot_i * ttot_i)
            return funcs

        self.CFDSolver.addUserFunction('my_force', my_force_functions, f_my_force)
        self.CFDSolver.addUserFunction('cfd_force', cfd_force_functions, f_cfd_force)
        self.CFDSolver.addUserFunction('my_power', my_power_functions, f_my_power)

        # set the az force
        az_force = 0.
        az_heat = 1e5
        self.ap.setDesignVars({"thrust": az_force, "heat":az_heat})

        self.CFDSolver(self.ap)
        funcs = {}
        funcsSens = {}
        # self.CFDSolver.evalFunctions(self.ap, funcs)
        self.CFDSolver.evalFunctions(self.ap, funcs, evalFuncs=['my_force', 'cfd_force', 'my_power'])
        self.CFDSolver.evalFunctionsSens(
            self.ap,
            funcsSens,
            evalFuncs=['my_force', 'cfd_force', 'my_power']
        )

        #####################
        # TEST MOMENTUM ADDED
        #####################

        # we test these w.r.t. thrust and heat analytically
        # The low accuracy is because the intgrated quantities don't have a lot of precision

        np.testing.assert_allclose(funcsSens['actuator_pipe_my_force']['thrust'], 1, rtol=1e-3)
        np.testing.assert_allclose(funcsSens['actuator_pipe_cfd_force']['thrust'], 1, rtol=1e-3)

        # heat addition should not affect these
        np.testing.assert_allclose(funcsSens['actuator_pipe_my_force']['heat'] + 100, 100, rtol=1e-3)
        np.testing.assert_allclose(funcsSens['actuator_pipe_cfd_force']['heat'] + 100, 100, rtol=1e-3)

        # also test the actual values from the ref file
        self.handler.root_print("my_force sens")
        self.handler.root_add_dict("my_force sens", funcsSens['actuator_pipe_my_force'], rtol=1e-12, atol=1e-12)

        self.handler.root_print("cfd_force sens")
        self.handler.root_add_dict("cfd_force sens", funcsSens['actuator_pipe_cfd_force'], rtol=1e-12, atol=1e-12)

        ##################
        # TEST POWER ADDED
        ##################

        # analytically test heat addition
        # this should be equal to one because we are not adding any thrust in this test.
        # if we also added thrust, addition of heat would affect flow power integration
        # due to the changes in the flowfield and as a result the derivative would not be one.
        # flowpower is more complicated here so we just check with json reference
        np.testing.assert_allclose(funcsSens['actuator_pipe_my_power']['heat'], 1, rtol=1e-3)

        # test values with the ref file
        self.handler.root_print("my_power sens")
        self.handler.root_add_dict("my_power sens", funcsSens['actuator_pipe_my_power'], rtol=1e-12, atol=1e-12)

    def test_actuator_flowpower_adjoint(self):
        "we test this adjoint separately because we need to have a finite thrust for this to actually test"

        # set the az force
        az_force = 600.
        az_heat = 1e5
        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        # We dont need to rerun, restart file has the correct state. just run a residual
        # self.CFDSolver(self.ap)
        self.CFDSolver.getResidual(self.ap)

        funcs = {}
        funcsSens = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)
        self.CFDSolver.evalFunctionsSens(
            self.ap,
            funcsSens,
            evalFuncs=['flowpower_az']
        )

        #############################
        # TEST FLOW POWER INTEGRATION
        #############################

        self.handler.root_print("flowpower sens")
        self.handler.root_add_dict("flowpower sens", funcsSens['actuator_pipe_flowpower_az'], rtol=1e-12, atol=1e-12)

    def test_actuator_partials(self):

        az_force = 600.
        az_heat = 1e5
        # need to set all dvs because training may re-use leftover dvs from a previous test
        self.ap.setDesignVars({"thrust": az_force, "heat": az_heat})

        # We dont need to rerun, restart file has the correct state. just run a residual
        # self.CFDSolver(self.ap)
        self.CFDSolver.getResidual(self.ap)

        #############
        # TEST FWD AD
        #############

        # save these for the dot product tests
        resDot = {}

        for DV in self.ap.DVs.values():
            # regular aero DVs do not need the family name, but the BC names do require the fam names
            xDvDot = {DV.key + '_actuator_region': 1}

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

class ActuatorDerivTests(unittest.TestCase):

    def setUp(self):
        gridFile = os.path.join(baseDir, '../../inputFiles/actuator_test_pipe.cgns')
        self.options = {'gridfile': gridFile,
                        "outputdirectory": os.path.join(baseDir, "../output_files"),
                        'mgcycle':'sg',
                        'ncycles':1,
                        'useanksolver':True,
                        # 'usenksolver':True,
                        'volumevariables': ['temp', 'mach', 'resrho' ],
                        'surfacevariables':['temp', 'vx', 'vy', 'vz', 'p', 'ptloss', 'mach', 'rho'],
                        'equationType':'Euler',
                        'l2convergence': 1e-13}

        CFDSolver = ADFLOW(options=self.options)

        CFDSolver.addFunction('mdot', 'inlet', name="mdot_in")
        CFDSolver.addFunction('mdot', 'outlet', name="mdot_out")

        CFDSolver.addFunction('mavgptot', 'outlet', name="mavgptot_out")
        CFDSolver.addFunction('mavgptot', 'inlet', name="mavgptot_in")

        CFDSolver.addFunction('aavgptot', 'outlet', name="aavgptot_out")
        CFDSolver.addFunction('aavgptot', 'inlet', name="aavgptot_in")

        CFDSolver.addFunction('mavgttot', 'outlet', name="mavgttot_out")
        CFDSolver.addFunction('mavgttot', 'inlet', name="mavgttot_in")

        CFDSolver.addFunction('mavgps', 'outlet', name="mavgps_out")
        CFDSolver.addFunction('mavgps', 'inlet', name="mavgps_in")

        CFDSolver.addFunction('aavgps', 'outlet', name="aavgps_out")
        CFDSolver.addFunction('aavgps', 'inlet', name="aavgps_in")

        CFDSolver.addFunction('area', 'inlet', name="area_in")
        CFDSolver.addFunction('area', 'outlet', name="area_out")

        CFDSolver.addFunction('mavgvx', 'inlet', name="vx_in")
        CFDSolver.addFunction('mavgvx', 'outlet', name="vx_out")

        CFDSolver.addFunction('mavgps', 'inlet', name="ps_in")
        CFDSolver.addFunction('mavgps', 'outlet', name="ps_out")

        self.CFDSolver = CFDSolver

        self.ap = AeroProblem(name='actuator_in_pipe', alpha=00, mach=0.6, altitude=0,
                        areaRef=1.0, chordRef=1.0,
                    evalFuncs=['mdot_in', 'mdot_out',
                               'mavgptot_in', 'mavgptot_out',
                               'mavgttot_in', 'mavgttot_out',
                               'mavgps_in', 'mavgps_out',
                               'area_in', 'area_out',
                               'aavgps_in', 'aavgps_out',
                               'aavgptot_in', 'aavgptot_out',
                               'ps_in', 'ps_out',
                               'vx_in', 'vx_out'] )

        self.ap.setActuatorVar('Thrust',  599.0, 'actuator')
        self.ap.addDV('Thrust', familyGroup='actuator', name='actuator_thrust')
        self.ap.setActuatorVar('Torque',  599.0, 'actuator')
        self.ap.addDV('Torque', familyGroup='actuator', name='actuator_torque')

        actuatorFile = os.path.join(baseDir, '../../inputFiles/actuator_test_disk.xyz')
        self.CFDSolver.addActuatorRegion(actuatorFile, np.array([0,0,0]),np.array([1,0,0]), 'actuator', thrust=600 )

    # def test_fwd_CS(self):
    #     # arch=self.arch['complex']
    #     import petsc4py
    #     petsc4py.init(arch='complex-debug')
    #     from petsc4py import PETSc
    #     from python.pyADflow_C import ADFLOW_C

    #     self.options['useanksolver'] = False
    #     CFDSolver = ADFLOW_C(options=self.options)

    #     CFDSolver.addFamilyGroup('upstream',['inlet'])
    #     CFDSolver.addFamilyGroup('downstream',['outlet'])
    #     CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
    #     CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

    #     CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
    #     CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")

    #     CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
    #     CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")

    #     CFDSolver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
    #     CFDSolver.addFunction('aavgptot', 'upstream', name="aavgptot_up")

    #     CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
    #     CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")

    #     CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
    #     CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")

    #     CFDSolver.addFunction('aavgps', 'downstream', name="aavgps_down")
    #     CFDSolver.addFunction('aavgps', 'upstream', name="aavgps_up")

    #     actuatorFile = os.path.join(baseDir, '../../inputFiles/actuator_test_disk.xyz')
    #     CFDSolver.addActuatorRegion(actuatorFile, np.array([0,0,0]),np.array([1,0,0]), 'actuator', thrust=600 )
    #     CFDSolver(self.ap, writeSolution=False)

    #     for DV in self.ap.DVs:
    #         xDvDot = {DV:1}

    #         resDot, funcsDot, fDot, hfDot = CFDSolver.computeJacobianVectorProductFwd(
    #         xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='CS')
    #         print(np.real(resDot[resDot != 0 ]))
    #         print(np.real(np.sum(resDot)))

    #         # np.testing.assert_allclose(resDot_FD,resDot, rtol=1e-8)
    #         # for func in funcsDot:
    #         #     np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], rtol=1e-8)

    #         # np.testing.assert_allclose(fDot_FD,fDot, rtol=1e-8)
    #         # np.testing.assert_allclose(hfDot_FD,hfDot, rtol=1e-8)



    #     del(CFDSolver)
    #     del(ADFLOW_C)

    #     # print(resDot)
    #     # print(funcsDot)
    #     # print(fDot)
    #     # print(hfDot)

if __name__ == '__main__':
    unittest.main()