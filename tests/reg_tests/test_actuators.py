from __future__ import print_function
import unittest
import numpy as np
from baseclasses import AeroProblem
import sys
import os
baseDir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(baseDir,'../../'))
from adflow.pyADflow import ADFLOW
from pprint import pprint as pp
import copy


class ActuatorBasicTests(unittest.TestCase):

    def setUp(self):
        gridFile = os.path.join(baseDir, '../../inputFiles/actuator_test_pipe.cgns')
        self.options = {'gridfile': gridFile,
                        "outputdirectory": os.path.join(baseDir, "../output_files"),
                        'mgcycle':'sg',
                        'ncycles':1000,
                        'useanksolver':True,
                        'usenksolver':True,
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


        self.force = 600
        actuatorFile = os.path.join(baseDir, '../../inputFiles/actuator_test_disk.xyz')
        self.CFDSolver.addActuatorRegion(actuatorFile, np.array([0,0,0]),np.array([1,0,0]), 'actuator', thrust=self.force )

    def test_mom_flux(self):
        "Tests if the correct amount of momentum is added to the flow by the actuator"

        self.CFDSolver(self.ap)
        funcs = {}
        self.CFDSolver.evalFunctions(self.ap, funcs)

        force_cfd =  -funcs[self.ap.name + '_mdot_out']*funcs[self.ap.name + '_vx_out'] + funcs[self.ap.name + '_ps_out']*funcs[self.ap.name +'_area_out'] \
                 - ( funcs[self.ap.name + '_mdot_in']*funcs[self.ap.name + '_vx_in'] + funcs[self.ap.name + '_ps_in']*funcs[self.ap.name +'_area_in'])

        # The low accuracy is because the intgrated quantities don't have a lot of precision
        np.testing.assert_allclose(force_cfd, self.force, rtol=1e-3)

    def test_set_thrust(self):
        "Tests if thrust is set accurately"

        force = 100

        ap = copy.deepcopy(self.ap)
        ap.setActuatorVar('Thrust',  force, 'actuator')
        ap.addDV('Thrust', familyGroup='actuator', name='actuator_thrust')
        self.CFDSolver(ap)
        funcs = {}
        self.CFDSolver.evalFunctions(ap, funcs)


        force_cfd =  -funcs[self.ap.name + '_mdot_out']*funcs[self.ap.name + '_vx_out'] + funcs[self.ap.name + '_ps_out']*funcs[self.ap.name +'_area_out'] \
                 - ( funcs[self.ap.name + '_mdot_in']*funcs[self.ap.name + '_vx_in'] + funcs[self.ap.name + '_ps_in']*funcs[self.ap.name +'_area_in'])

        # The low accuracy is because the intgrated quantities don't have a lot of precision
        np.testing.assert_allclose(force_cfd, force, rtol=1e-3)


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