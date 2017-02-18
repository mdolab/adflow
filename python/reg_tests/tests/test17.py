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
                 areaRef=areaRef, chordRef=chordRef,
                 evalFuncs=['mdot_up', 'mdot_down', 'mdot_plane', 
                            'mavgptot_up', 'mavgptot_down', 'mavgptot_plane', 
                            'mavgttot_up', 'mavgttot_down', 'mavgttot_plane',
                            'mavgps_up', 'mavgps_down', 'mavgps_plane',     
                            'pk_up', 'pk_down', 'pk_plane',
                            'edot_up', 'edot_down', 'edot_plane',
                            'sigmamn_up',  'sigmamn_plane',
                            'sigmaptot_up', 'sigmaptot_plane',
                            ])


ap.setBCVar('Pressure',  79326.7, 'downstream')
ap.addDV('Pressure', family='downstream')

ap.setBCVar('PressureStagnation',  100000.0, 'upstream')
ap.addDV('PressureStagnation', family='upstream')

ap.setBCVar('TemperatureStagnation',  500.0, 'upstream')
ap.addDV('TemperatureStagnation', family='upstream')

 
CFDSolver = ADFLOW(options=options, debug=True)

CFDSolver.addIntegrationSurface('../inputFiles/integration_plane_viscous.fmt', 'viscous_plane')

CFDSolver.addFamilyGroup('upstream',['inlet'])
CFDSolver.addFamilyGroup('downstream',['outlet'])
CFDSolver.addFamilyGroup('all_flow',['inlet', 'outlet'])
CFDSolver.addFamilyGroup('output_fam',['all_flow', 'allWalls'])

CFDSolver.addFunction('mdot', 'upstream', name="mdot_up")
CFDSolver.addFunction('mdot', 'downstream', name="mdot_down")
CFDSolver.addFunction('mdot', 'viscous_plane', name="mdot_plane")

CFDSolver.addFunction('mavgptot', 'downstream', name="mavgptot_down")
CFDSolver.addFunction('mavgptot', 'upstream', name="mavgptot_up")
CFDSolver.addFunction('mavgptot', 'viscous_plane', name="mavgptot_plane")

CFDSolver.addFunction('mavgttot', 'downstream', name="mavgttot_down")
CFDSolver.addFunction('mavgttot', 'upstream', name="mavgttot_up")
CFDSolver.addFunction('mavgttot', 'viscous_plane', name="mavgttot_plane")

CFDSolver.addFunction('mavgps', 'downstream', name="mavgps_down")
CFDSolver.addFunction('mavgps', 'upstream', name="mavgps_up")
CFDSolver.addFunction('mavgps', 'viscous_plane', name="mavgps_plane")

CFDSolver.addFunction('pk', 'downstream', name="pk_down")
CFDSolver.addFunction('pk', 'upstream', name="pk_up")
CFDSolver.addFunction('pk', 'viscous_plane', name="pk_plane")

CFDSolver.addFunction('edot', 'downstream', name="edot_down")
CFDSolver.addFunction('edot', 'upstream', name="edot_up")
CFDSolver.addFunction('edot', 'viscous_plane', name="edot_plane")

CFDSolver.addFunction('sigmamn', 'upstream', name="sigmamn_up")
CFDSolver.addFunction('sigmamn', 'viscous_plane', name="sigmamn_plane")

CFDSolver.addFunction('sigmaptot', 'upstream', name="sigmaptot_up")
CFDSolver.addFunction('sigmaptot', 'viscous_plane', name="sigmaptot_plane")


CFDSolver.setOption('ncycles',1000)

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
    print 'Eval Functions:'
    reg_write_dict(funcs, 1e-10, 1e-10)


