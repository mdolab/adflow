# built-ins
import unittest
import numpy as np
import os
import copy

# MACH classes
from adflow import ADFLOW

# import the testing utilities that live a few directories up
import reg_test_utils as utils
from reg_default_options import adflowDefOpts
from reg_aeroproblems import ap_actuator_pipe
import reg_test_classes

baseDir = os.path.dirname(os.path.abspath(__file__))

class ActuatorBasicTests(unittest.TestCase):
    """
    Tests for the actuator zone.
    """

    N_PROCS = 2

    def setUp(self):

        self.options = {'gridfile': os.path.join(baseDir, '../../inputFiles/actuator_test_pipe.cgns'),
                        # "outputdirectory": os.path.join(baseDir, "../output_files"),
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
                        'l2convergence': 1e-13}

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

        # add thrust as an AP DV
        self.ap.setBCVar("Thrust", 0.0, "actuator_region")
        self.ap.addDV("Thrust", family="actuator_region", units="N", name="thrust")

        # also add flowpower as an AZ function
        CFDSolver.addFunction("flowpower", "actuator_region", name="flowpower_az")

    def test_actuator_momentum_flux(self):
        "Tests if the correct amount of momentum is added to the flow by the actuator"

        # set the az force
        az_force = 600.
        self.ap.setDesignVars({"thrust": az_force})

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


    def test_fwd(self):
        self.CFDSolver(self.ap)
        "note the torque dv does not effect anything so the deriv is and should be zero"
        for DV in self.ap.DVs:
            xDvDot = {DV:1}

            resDot_FD, funcsDot_FD, fDot_FD, hfDot_FD = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='FD', h=1e-3)

            resDot, funcsDot, fDot, hfDot = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True, mode='AD')


            np.testing.assert_allclose(resDot_FD,resDot, atol=1e-10)
            for func in funcsDot:
                np.testing.assert_allclose(funcsDot_FD[func],funcsDot[func], atol=1e-7)

            np.testing.assert_allclose(fDot_FD,fDot, atol=1e-7)
            np.testing.assert_allclose(hfDot_FD,hfDot, atol=1e-10)


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


    def test_bwd(self):
        """ tests the bwd AD for actuator regions """
        self.CFDSolver(self.ap)

        for DV in self.ap.DVs:
            xDvDot = {DV:1}

            resDot, _, _, _ = self.CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True, hfDeriv=True)

            dwBar = self.CFDSolver.getStatePerturbation(123)
            _, _, xDvBar = self.CFDSolver.computeJacobianVectorProductBwd(
                resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

            np.testing.assert_array_almost_equal(np.dot(xDvBar[DV],xDvDot[DV]),\
                                                    np.dot(resDot, dwBar), decimal=14)


if __name__ == '__main__':
    unittest.main()