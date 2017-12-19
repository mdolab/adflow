from __future__ import print_function
############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI
from baseclasses import AeroProblem
sys.path.append(os.path.abspath('./'))
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

gridFile = '../inputFiles/conic_conv_nozzle.cgns'

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
    'MGCycle':'sg',
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
                 areaRef=areaRef, chordRef=chordRef,
                 evalFuncs=['mdot_up', 'mdot_down', #'mdot_plane',
                            'mavgptot_up', 'mavgptot_down',# 'mavgptot_plane',
                            'mavgttot_up', 'mavgttot_down',# 'mavgttot_plane',
                            'mavgps_up', 'mavgps_down', #'mavgps_plane'
                            'sigmamn_up',  #'sigmamn_plane'
                            'sigmaptot_up', #'sigmaptot_plane'
                            ])


ap.setBCVar('Pressure',  79326.7, 'downstream')
ap.addDV('Pressure', family='downstream')

ap.setBCVar('PressureStagnation',  100000.0, 'upstream')
ap.addDV('PressureStagnation', family='upstream')

ap.setBCVar('TemperatureStagnation',  500.0, 'upstream')
ap.addDV('TemperatureStagnation', family='upstream')


CFDSolver = ADFLOW(options=options, debug=True)

#CFDSolver.addIntegrationSurface('integration_plane.fmt', 'coarse_plane')
#CFDSolver.addIntegrationSurface('integration_plane_fine.fmt', 'fine_plane')
#CFDSolver.addIntegrationSurface('integration_plane_viscous.fmt', 'viscous_plane')

CFDSolver.addFamilyGroup('upstream',['inlet'])
CFDSolver.addFamilyGroup('downstream',['outlet'])
CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
#CFDSolver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
#CFDSolver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
#CFDSolver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")
#CFDSolver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

CFDSolver.addFunction('sigmamn', 'upstream', name="sigmamn_up")
CFDSolver.addFunction('sigmaptot', 'upstream', name="sigmaptot_up")


CFDSolver.setOption('ncycles',1000)

# Check the residual
res = CFDSolver.getResidual(ap)
#TODO: getResNorms() doesn't work for overset?
# totalR0, totalRStart, totalRFinal = CFDSolver.getResNorms()
# print res, totalR0, totalRStart, totalRFinal
# res /= totalR0

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
