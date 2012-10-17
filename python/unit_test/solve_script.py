# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys

# =============================================================================
# External Python modules
# =============================================================================
import numpy
# =============================================================================
# Extension modules
# =============================================================================

from mdo_import_helper import import_modules, MPI

# Import PETSc so it is initialized on ALL Procs
import petsc4py 
petsc4py.init()

# MDO_Lab Imports
exec(import_modules('pyAero_problem','pyAero_reference','pyAero_geometry'))
exec(import_modules('pyAero_flow'))
exec(import_modules('pySUMB'))

# Pull out the options. These are given by run_unit_test.py so we
# don't have to particular about checking them
grid_file = sys.argv[1]
sol_typ = sys.argv[2]
solver = sys.argv[3]
eq_mode = sys.argv[4]
adj_type = sys.argv[5]

if grid_file == 'bump':
    grid_file = '../examples/Bump_CGNS'
    Mach    = 0.8395
    Area_ref = 1.0
    Span_ref = 1.0
    Chord_ref =  1.0
    mgcycle = 'sg'
else:
    grid_file = '../examples/oneram6_2k'
    Mach    = 0.8395
    Area_ref = .772893541
    Span_ref = 14.2
    Chord_ref =  .64607
    mgcycle = 'sg'
# end if

# Do some processing on options
if solver == 'NK':
    use_NK = True
else:
    use_NK = False
# end if

if eq_mode == 'TS':
    eq_mode = 'Time Spectral'
else:
    eq_mode = 'Steady'
# end if

if adj_type == 'fwd':
    useReverse = False
else:
    useReverse = True

# ================================================================
#                   INPUT INFORMATION  

output_directory = './'
gcomm = MPI.COMM_WORLD

# ================================================================
#               Set Options for each solver

aeroOptions = {
    # Common Paramters
    'gridFile':grid_file+'.cgns',
    'outputDir':'./',
    'writeSolution':False,

    # Physics Paramters
    'vis4':.0156,
    'vis2':0.50,
    'vis2Coarse':0.50, 
    'restrictionRelaxation':1.0,
    'equationMode': eq_mode,

    # Common Paramters
    'CFL':1.0, # 1.0 for TS
    'CFLCoarse':1.0,
    'MGCycle':mgcycle,
    'MGStartLevel':-1, # Start on coarest
    
    # Convergence Paramters
    'L2Convergence':1e-8,
    'L2ConvergenceCoarse':.01,

    # Newton-Krylov Paramters
    'useNKSolver':use_NK,
    'NKLinearSolver':'gmres',
    'NKSwitchTol':.1,
    'NKSubspaceSize':50,
    'NKPC':'Additive Schwartz',
    'NKASMOverlap':1,
    'NKPCILUFill':1,
    'NKLocalPCOrdering':'RCM',
    'NKJacobianLag':10,

    # Load Balance Paramters
    'blockSplitting':True,
    'loadImbalance':0.1,
       
    # Misc Paramters
    'printIterations':True,
    'printTiming':False,
    'monitorVariables':['resrho','cl','cd','totalR'],
    'surfaceVariables':['vx','vy','vz','rho','P','mach','cp'],

    # Time Spectral Paramters
    'timeIntervals': 3,
    'alphaMode':False,
    'betaMode':False,
    'machMode':False,
    'pMode':False,
    'qMode':True,
    'rMode':False,
    'altitudeMode':False,
    'windAxis':False,
    'rotCenter':[0.0,0.0,0.0],
    'TSStability': False,

    # Adjoint Paramters
    'adjointL2Convergence':1e-8,
    'approxPC': True,
    'restartAdjoint':False,
    'adjointSolver': 'GMRES',
    'adjointMaxIter': 500,
    'adjointSubspaceSize' : 60,
    'adjointMonitorStep': 10,
    'preconditionerSide': 'RIGHT',
    'matrixOrdering': 'RCM',
    'globalPreconditioner': 'Additive Schwartz',
    'localPreconditioner' : 'ILU',
    'ILUFill':1,
    'ASMOverlap':1,
    'finiteDifferencePC':True,
    'useReverseModeAD':useReverse,
    }

#  Setup AeroProblem
flow = Flow(name='Base Case',mach=Mach,alpha=3.06,beta=0.0,liftIndex=2,
            degreeFourier=1, omegaFourier=6.28, cosCoefFourier=[0.0,0.0],
            sinCoefFourier=[0.01])#34]) 

ref = Reference('Baseline Reference',Area_ref,Span_ref,Chord_ref) #area,span,chord 
aeroProblem = AeroProblem(name='SUMB Test',flow_set=flow,ref_set=ref)

#Setup the Solver
CFDsolver = SUMB(comm=gcomm, init_petsc=False, options=aeroOptions, mesh=None)
CFDsolver.initialize(aeroProblem)

# Add aeroDVs
CFDsolver.addAeroDV('aofa')
CFDsolver.addAeroDV('mach')

# Solve Problem
CFDsolver(aeroProblem,500)

# Get the solution and print:
sol = CFDsolver.getSolution()
if gcomm.rank == 0:
    for key in sorted(sol.keys()):
        print 'sol[%s]:'%(key),sol[key]
    # end for
# end if

# We are going to run as many things as possible if for no other
# reason to make sure they run

# Get the forces...these are the sumb forces:
forces = CFDsolver.getForces()

# Each processor computes sum of their forces
F_sum = numpy.sum(forces.flatten())

# Gather to root: (since order when printing in paralell is undefined)
F = gcomm.gather(F_sum, root=0)
if gcomm.rank == 0:
    print 'proc, sum(F)'
    for i in xrange(len(F)):
        print i, F[i]


# Solve Adjoint(s)
CFDsolver.solveAdjoint('cl')

# Get res norms:
r, rstart, rfinal = CFDsolver.getResNorms()
if gcomm.rank == 0:
    print 'Res norm, Start Res Norm, Final Res Norm:',r, rstart, rfinal

# Get total aero derivative
dIda = CFDsolver.totalAeroDerivative('cl')
if gcomm.rank == 0:
    print 'Aero Derivatives:'
    print dIda

# Get dRdxvPsi since we can't do the total spatial derivative:
CFDsolver.setupSpatialMatrices()
drdxvpsi = CFDsolver.getdRdXvPsi(None, 'cl')
nrm = numpy.linalg.norm(drdxvpsi)
nrm = gcomm.gather(nrm, root=0)
if gcomm.rank == 0:
    print 'proc, norm of component of dRdxv^T * psi'
    for i in xrange(len(nrm)):
        print i, nrm[i]
    # end for
# end if

# get dIdx
didx = CFDsolver.getdIdx('cl')
nrm = numpy.linalg.norm(didx)
nrm = gcomm.gather(nrm, root=0)
if gcomm.rank == 0:
    print 'proc, norm of component of dIdx'
    for i in xrange(len(nrm)):
        print i, nrm[i]
    # end for
# end if

# Setup coupling matrices
CFDsolver.setupCouplingMatrices()

# Reset Adjoint
CFDsolver.resetAdjoint('cl')

# Rest flow
CFDsolver.resetFlow()

# Dump Adjoint memory
CFDsolver.releaseAdjointMemory()
