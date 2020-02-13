############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE SUMB
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI
from baseclasses import AeroProblem
from idwarp import USMesh
from mdo_regression_helper import *
from commonUtils import *

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from sumb import SUMB' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pySUmb import SUMB

printHeader('Test 15d: NACA 0012 2D Time-Accurate, Forced motion, Warping of Mesh, ALE external loop - DADI Smoother')
gridFile = '../inputFiles/naca0012_rans-L2.cgns'

options = copy.copy(sumbDefOpts)

k = 0.0808
M = 0.6
gamma = 1.4
R = 287.085
T = 280.0
c = 1.0
alpha_m = 2.77 # 2.89 #2.77 #Modified numbers
alpha_0 = 2.34 # 2.41 #2.34

omega = 2*M*numpy.sqrt(gamma*R*T)*k/c 
deltaAlpha = -alpha_0*numpy.pi/180.0 

# Set forcing frequency and other information
f = omega/(2*numpy.pi) # [Hz] Forcing frequency of the flow
period = 1.0/f # [sec]
nStepPerPeriod = 8
nPeriods = 1
nfineSteps = nStepPerPeriod*nPeriods
dt = period / nStepPerPeriod # [s] The actual timestep

name = '0012pitching'

# Mesh object
optMesh = {
    'gridFile':           gridFile,
    'warpType':           'algebraic',
    'solidSolutionType':  'linear', # or nonLinear or steppedLinear
    'fillType':           'linear', # or cubic
    'solidWarpType':      'n',
    'n':                  7,
    
    # Solution Parameters
    'nSolutionSteps':     5,      # only for nonLinear or steppedLinear
    'ksp_its':            25,
    'warp_tol':           1e-10,  # overall tolerance for nonLinear/steppedLinear 
    'reduction_factor_for_reform': 0.5, 
    'useSolutionMonitor': False,
    'useKSPMonitor':      False,
    
    # Physical parameters
    'nu':                 0.0,
    'e_exp':              1.0,
    'skew_exp':           0.0,
    'useRotationCorrection': False,
    }
mesh = USMesh(options = optMesh)

options.update(
    {'gridfile': '../inputFiles/naca0012_rans-L2.cgns',
     'writevolumesolution':False,
     'vis4':.025,
     'vis2':0.5,
     'restrictionrelaxation':.5,
     'smoother':'dadi',
     'equationtype':'RANS',
     'equationmode':'unsteady',
     'timeIntegrationscheme':'bdf',
     'ntimestepsfine':nfineSteps,
     'deltat':dt,
     'nsubiterturb':10,
     'nsubiter':5,         
     'useale':True,
     'usegridmotion':False,
     'coupledsolution':True,
     'cfl':2.5,
     'cflcoarse':1.2,
     'ncycles':2000,            
     'mgcycle':'3w',
     'mgstartlevel':1,
     'monitorvariables':['cpu','resrho','cl','cd','cmz'],                            
     'usenksolver':False,
     'l2convergence':1e-6,
     'l2convergencecoarse':1e-4,
     'qmode':False,
     'alphafollowing': False,
     'blockSplitting':True,
     }
)
ap = AeroProblem(name=name, alpha=alpha_m,  mach=M, machRef=M, reynolds=4800000.0,reynoldsLength=c, T=T, R=R,
                 areaRef=1.0, chordRef=c, evalFuncs=['cl','cd','cmz'],xRef=0.25,xRot=0.25)

def callback(refGrid, t, ts):
    newGrid = numpy.copy(refGrid)
    x = refGrid[:, 0]
    y = refGrid[:, 1]
    p = deltaAlpha*numpy.sin(omega*t)
    c = numpy.cos(p)
    s = numpy.sin(p)
    newGrid[:, 0] = c * (x-0.25) - s * y + 0.25
    newGrid[:, 1] = s * (x-0.25) + c * y
    return newGrid

CFDSolver = SUMB(options=options)
CFDSolver.setMesh(mesh)
CFDSolver.addSlices('z',[0.5])
CFDSolver(ap)

refCoor = CFDSolver.getSurfaceCoordinates('allWalls')
for tdx in xrange(1, nfineSteps+1):
    curTime, curTimeStep = CFDSolver.advanceTimeStepCounter()
    newCoor = callback(refCoor, curTime, curTimeStep)
    # Set displacements
    CFDSolver.setSurfaceCoordinates(newCoor, 'allwalls')
    CFDSolver.updateGeometryInfo()
    # Solve current time step
    CFDSolver.solveTimeStep()

funcs = {}
CFDSolver.evalFunctions(ap, funcs)
CFDSolver.checkSolutionFailure(ap, funcs)
if MPI.COMM_WORLD.rank == 0:
    print 'Eval Functions:'
    reg_write_dict(funcs, 1e-6, 1e-6)

os.system('rm  0012pitching*')
