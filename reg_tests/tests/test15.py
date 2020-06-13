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
from adflow import ADFLOW

printHeader('Test 15: NACA 0012 2D Time-Accurate, Forced motion, Rigid Rotation of Mesh - DADI Smoother')
gridFile = '../inputFiles/naca0012_rans-L2.cgns'

options = copy.copy(adflowDefOpts)

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
f = 10.0 # [Hz] Forcing frequency of the flow
period = 1.0/f # [sec]
nStepPerPeriod = 8
nPeriods = 1
nfineSteps = nStepPerPeriod*nPeriods
dt = period / nStepPerPeriod # [s] The actual timestep

name = '0012pitching'

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
     'useale':False,
     'usegridmotion':True,
     'cfl':2.5,
     'cflcoarse':1.2,
     'ncycles':2000,
     'mgcycle':'3w',
     'mgstartlevel':1,
     'monitorvariables':['cpu','resrho','cl','cd','cmz'],
     'usenksolver':False,
     'l2convergence':1e-6,
     'l2convergencecoarse':1e-4,
     'qmode':True,
     'alphafollowing': False,
     'blockSplitting':True,
     'useblockettes':False,
     }
)
ap = AeroProblem(name=name, alpha=alpha_m,  mach=M, machRef=M, reynolds=4800000.0,reynoldsLength=c, T=T, R=R,
                 areaRef=1.0, chordRef=c, evalFuncs=['cl','cd','cmz'],xRef=0.25,xRot=0.25,
                 degreePol=0,coefPol=[0.0],degreeFourier=1,omegaFourier=omega,
                 cosCoefFourier=[0.0,0.0],sinCoefFourier=[deltaAlpha])

def setup_cb(comm): 
    CFDSolver = ADFLOW(options=options, comm=comm)
    CFDSolver.addSlices('z',[0.5])
    CFDSolver(ap)

    return CFDSolver, None, None, None

if __name__ == "__main__":

    CFDSolver, _, _, _ = setup_cb(MPI.COMM_WORLD)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print('Eval Functions:')
        reg_write_dict(funcs, 1e-6, 1e-6)

    os.system('rm  0012pitching*')
