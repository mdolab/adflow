from __future__ import print_function
############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI
from baseclasses import AeroProblem
from mdo_regression_helper import *
from commonUtils import *

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from adflow import ADFLOW' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pyADflow import ADFLOW
# ###################################################################

# ======================================================================
#         Input Information -- Modify accordingly!
# ======================================================================
outputDirectory = './'

gridFile = '../inputFiles/conic_conv_nozzle_mb.cgns'
restartFile = None #'restart_mb.cgns'
mgcycle = '3w'


options = copy.copy(adflowDefOpts)

options = {
    # Common Parameters
    'gridFile':gridFile,
    # Physics Parameters
    'equationType':'euler',
    'smoother':'dadi',
    'nsubiter':3,
    'CFL':4.0,
    'CFLCoarse':1.25,
    'MGCycle':'3w',
    'MGStartLevel':-1,
    'nCyclesCoarse':250,
    'nCycles':1000,
    'nkcfl0':1e10,
    'monitorvariables':['cpu', 'resrho','cl','cd'],
    'volumevariables':['blank'],
    'surfacevariables':['mach', 'cp', 'vx', 'vy','vz', 'blank'],
    'useNKSolver':True,
    'nkswitchtol':.01,
    'nkadpc':True,
    'nkjacobianlag':5,
    'nkouterpreconits':3,
    'nkinnerpreconits':2,
    # Convergence Parameters
    'L2Convergence':1e-10,
    'L2ConvergenceCoarse':1e-4,
    'adjointl2convergence':1e-6,
    'forcesAsTractions':True,
    'debugzipper':True,
    'nearwalldist':.001,
    'nkls':'none',
    'solutionprecision':'double',
    'adjointsubspacesize':200,
    'outerpreconits':3,
    'zipperSurfaceFamily':'output_fam',
    'flowtype':'internal',
    }

solve = True
if 'solve' not in sys.argv:
    options['restartfile'] = gridFile
    solve = False

alpha = 90.0
mach = 0.5
areaRef = 1.0
chordRef = 1.0
altitude = 0
name = 'nozzle'

# Aerodynamic problem description
ap = AeroProblem(name=name, alpha=alpha, mach=mach, altitude=altitude,
                 areaRef=areaRef, chordRef=chordRef, R=287.87,
                 evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane',
                            'mavgptot_up', 'mavgptot_down', 'mavgptot_plane',
                            'aavgptot_up', 'aavgptot_down', 'aavgptot_plane',
                            'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                            'mavgps_up', 'mavgps_down', 'mavgps_plane',
                            'aavgps_up', 'aavgps_down', 'aavgps_plane',
                            ])


ap.setBCVar('Pressure',  79326.7, 'downstream')
ap.addDV('Pressure', family='downstream')

ap.setBCVar('PressureStagnation',  100000.0, 'upstream')
ap.addDV('PressureStagnation', family='upstream')

ap.setBCVar('TemperatureStagnation',  500.0, 'upstream')
ap.addDV('TemperatureStagnation', family='upstream')


def setup_cb(comm):

    solver = ADFLOW(options=options, comm=comm, debug=True)


    solver.addIntegrationSurface('../inputFiles/integration_plane_viscous.fmt', 'viscous_plane')
    solver.finalizeUserIntegrationSurfaces()

    solver.addFamilyGroup('upstream',['inlet'])
    solver.addFamilyGroup('downstream',['outlet'])
    solver.addFamilyGroup('all_flow',['inlet', 'outlet'])
    solver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

    solver.addFunction('mdot', 'upstream', name="mdot_up")
    solver.addFunction('mdot', 'downstream', name="mdot_down")
    solver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

    solver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
    solver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
    solver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

    solver.addFunction('aavgptot', 'downstream', name="aavgptot_down")
    solver.addFunction('aavgptot', 'upstream', name="aavgptot_up")
    solver.addFunction('aavgptot', 'viscous_plane', name="aavgptot_plane")

    solver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
    solver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
    solver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

    solver.addFunction('mavgps', 'downstream', name="mavgps_down")
    solver.addFunction('mavgps', 'upstream', name="mavgps_up")
    solver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

    solver.addFunction('aavgps', 'downstream', name="aavgps_down")
    solver.addFunction('aavgps', 'upstream', name="aavgps_up")
    solver.addFunction('aavgps', 'viscous_plane', name="aavgps_plane")

    solver.setOption('ncycles',1000)

    return solver, None, None, None


if __name__ == "__main__":

    CFDSolver, _, _, _ = setup_cb(MPI.COMM_WORLD)

    CFDSolver(ap)

    # Check the residual
    res = CFDSolver.getResidual(ap)
    totalR0, totalRStart, totalRFinal = CFDSolver.getResNorms()
    res /= totalR0

    parPrint('Norm of residual')
    reg_par_write_norm(res, 1e-10, 1e-10)

    # Get and check the states
    parPrint('Norm of state vector')
    reg_par_write_norm(CFDSolver.getStates(), 1e-10, 1e-10)


    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print('Eval Functions:')
        reg_write_dict(funcs, 1e-10, 1e-10)


