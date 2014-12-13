############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE SUMB
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI
from mdo_regression_helper import *
from baseclasses import AeroProblem
from pygeo import DVGeometry
from pywarp import MBMesh_C as MBMesh
import pyspline

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from sumb import SUMB' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pySUmb_C import SUMB_C as SUMB
# ###################################################################

# First thing we will do is define a complete set of default options
# that will be reused as we do differnt tests.  These are the default
# options as Dec 8, 2014.

defOpts = {
    # Common Paramters
    'gridfile': 'default.cgns',
    'restartfile': 'default_restart.cgns',
    'solrestart': False,

    # Output Parameters
    'storerindlayer': True,
    'outputdirectory': './',
    'writesymmetry': True,
    'writefarfield': False,
    'writesurfacesolution':True,
    'writevolumesolution':True,
    'solutionprecision':'single',
    'gridprecision':'double',
    'isosurface': {},
    'isovariables': [],
    'viscoussurfacevelocities': True,

    # Physics Paramters
    'discretization': 'central plus scalar dissipation',
    'coarsediscretization': 'central plus scalar dissipation',
    'limiter': 'vanalbeda',
    'smoother': 'runge kutta',
    'equationtype':  'euler',
    'equationmode':  'steady',
    'flowtype': 'external',
    'turbulencemodel': 'sa',
    'turbulenceorder': 'first order',
    'usewallfunctions': False,
    'useapproxwalldistance': True,
    'walltreatment': 'linear pressure extrapolation',
    'dissipationscalingexponent': 0.67,
    'vis4': 0.0156,
    'vis2': 0.25,
    'vis2coarse': 0.5,
    'restrictionrelaxation': .80,
    'liftindex': 2,
    'lowspeedpreconditioner': False,
    'turbresscale': 10000.0,

    # Common Paramters
    'ncycles': 500,
    'ncyclescoarse': 500,
    'nsubiterturb': 1,
    'nsubiter': 1,
    'cfl': 1.7,
    'cflcoarse': 1.0,
    'mgcycle': '3w',
    'mgstartlevel': -1,
    'resaveraging':'alternateresaveraging',
    'smoothparameter': 1.5,
    'cfllimit': 1.5,

    # Unsteady Paramters
    'timeintegrationscheme': 'bdf',
    'timeaccuracy': 2,
    'ntimestepscoarse': 48,
    'ntimestepsfine': 400,
    'deltat': .010,
    
    # Time Spectral Paramters
    'timeintervals':  1,
    'alphamode': False,
    'betamode': False,
    'machmode': False,
    'pmode': False,
    'qmode': False,
    'rmode': False,
    'altitudemode': False,
    'windaxis': False,
    'tsstability':  False,

    # Convergence Paramters
    'l2convergence': 1e-6,
    'l2convergencerel': 1e-16,
    'l2convergencecoarse': 1e-2,
    'maxl2deviationfactor': 1.0,
    'coeffconvcheck': False,
    'miniterationnum': 0,

    # Newton-Krylov Paramters
    'usenksolver': False,
    'nklinearsolver': 'gmres',
    'nkswitchtol': 2.5e-4,
    'nksubspacesize': 60,
    'nklinearsolvetol': 0.3,
    'nkuseew': True,
    'nkpc': 'additive schwartz',
    'nkadpc': False,
    'nkviscpc': False,
    'nkasmoverlap': 1,
    'nkpcilufill': 2,
    'nklocalpcordering': 'rcm',
    'nkjacobianlag': 20,
    'rkreset': False,
    'nrkreset': 5,
    'applypcsubspacesize': 10,
    'nkinnerpreconits': 1,
    'nkouterpreconits': 1,
    'nkls': 'cubic',
    
    # Load Balance/partitioning parameters
    'blocksplitting': True,
    'loadimbalance': 0.1,
    'loadbalanceiter': 10,
    'partitiononly': False,

    # Misc Paramters
    'metricconversion': 1.0,
    'autosolveretry': False,
    'autoadjointretry': False,
    'storehistory': False,
    'numbersolutions': True,
    'printiterations': True,
    'printtiming': True,
    'setmonitor': True,
    'printwarnings': True,
    'monitorvariables': ['cpu','resrho','cl', 'cd'],
    'surfacevariables': ['cp','vx', 'vy','vz', 'mach'],
    'volumevariables': ['resrho'],
    
    # Multidisciplinary Coupling Parameters:
    'forcesastractions': True,

    # Adjoint Paramters
    'adjointl2convergence': 1e-6,
    'adjointl2convergencerel': 1e-16,
    'adjointl2convergenceabs': 1e-16,
    'adjointdivtol': 1e5,
    'approxpc':  True,
    'adpc':  False,
    'viscpc':False,
    'usediagtspc': True,
    'restartadjoint': True,
    'adjointsolver':  'gmres',
    'adjointmaxiter':  500,
    'adjointsubspacesize' :  100,
    'adjointmonitorstep':  10,
    'dissipationlumpingparameter': 6.0,
    'preconditionerside':  'right',
    'matrixordering':  'rcm',
    'globalpreconditioner':  'additive schwartz',
    'localpreconditioner' :  'ilu',
    'ilufill':  2,
    'asmoverlap' :  1,
    'innerpreconits': 1,
    'outerpreconits': 3,
    'usereversemodead': False,
    'applyadjointpcsubspacesize': 20,
    'frozenturbulence': True,
    'usematrixfreedrdw': False,
    'usematrixfreedrdx': False,

    # ADjoint debugger
    'firstrun': True,
    'verifystate': True,
    'verifyspatial': True,
    'verifyextra': True,
}
evalFuncs=['lift','drag','cl','cd','fx','fy','fz','cfx','cfy',
           'cfz','mx','my','mz','cmx','cmy','cmz','sepsensor']
# First thing we will test is the euler mesh of the MDO tutorial. This
# is a very small mesh that can be run very quickly. We therefore do a
# lot of testing with it.

def printHeader(testName):
    if MPI.COMM_WORLD.rank == 0:
        print '+' + '-'*78 + '+'
        print '| Test Name: ' + '%-66s'%testName + '|'
        print '+' + '-'*78 + '+'
h = 1e-40

def test1():
    # ****************************************************************************
    printHeader('MDO tutorial Euler Mesh Alpha Perturbation')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridFile': '../inputFiles/mdo_tutorial_euler.cgns',
         'MGCycle':'2w',
         'CFL':1.5,
         'CFLCoarse':1.25,
         'nCyclesCoarse':250,
         'nCycles':1000,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'useNKSolver':True,
         'L2Convergence':1e-12,
         'L2ConvergenceCoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-12,
         'nkls':'none',
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8+h*1j, mach=0.80,
                     altitude=10000.0, areaRef=45.5, chordRef=3.25,
                     evalFuncs=evalFuncs)
    CFDSolver = SUMB(options=aeroOptions)
    CFDSolver(ap,writeSolution=False)
    sol = {}
    CFDSolver.evalFunctions(ap,sol)
    if MPI.COMM_WORLD.rank == 0:
        for key in sorted(sol.keys()):
            print 'sol[%s]:'%key
            reg_write(numpy.real(sol[key]),1e-10,1e-10)
            deriv = numpy.imag(sol[key])/h
            reg_write(deriv,1e-6,1e-6)
  
    # Clean up:
    del CFDSolver

def test2():
    # ****************************************************************************
    printHeader('MDO tutorial Euler Mesh TS-Mode Alpha Perturbation')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridFile': '../inputFiles/mdo_tutorial_euler_l2.cgns',
         'MGCycle':'sg',
         'CFL':1.5,
         'CFLCoarse':1.25,
         'nCyclesCoarse':250,
         'nCycles':10000,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'L2Convergence':1e-12,
         'L2ConvergenceCoarse':1e-2,
         'adjointl2convergence': 1e-12,
         'timeIntervals':3,
         'equationMode':'Time Spectral',
         'tsstability':True,
         'alphamode':True,
         'useNKSolver':True,
         'NKSwitchTol':1e-2,
         'usediagtspc':False,
         'asmoverlap':2,
         'ilufill':3,
         'adjointmaxiter':1000,
     }
    )
    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8+h*1j, mach=0.80, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25, 
                     degreeFourier=1, omegaFourier=6.28, degreePol=0,
                     cosCoefFourier=[0.0,0.0], sinCoefFourier=[0.01], 
                     coefPol=[0],
                     evalFuncs=['cl', 'cl0','clalpha', 'clalphadot'])
    CFDSolver = SUMB(options=aeroOptions)
    CFDSolver(ap,writeSolution=False)
    sol = {}
    CFDSolver.evalFunctions(ap,sol)
    if MPI.COMM_WORLD.rank == 0:
        for key in sorted(sol.keys()):
            print 'sol[%s]:'%key
            reg_write(numpy.real(sol[key]),1e-10,1e-10)
            deriv = numpy.imag(sol[key])/h
            reg_write(deriv,1e-6,1e-6)
            
    # Clean up:
    del CFDSolver
  
def test3():
    # ****************************************************************************
    printHeader('MDO tutorial RANS Mesh Twist Perturbation')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridFile': '../inputFiles/mdo_tutorial_rans.cgns',
         'MGCycle':'2w',
         'equationType':'RANS',
         'smoother':'dadi',
         'CFL':1.5,
         'CFLCoarse':1.25,
         'resaveraging':'noresaveraging',
         'nsubiter':3,
         'nsubiterturb':3,
         'nCyclesCoarse':100,
         'nCycles':1000,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'useNKSolver':True,
         'L2Convergence':1e-15,
         'L2ConvergenceCoarse':1e-4,
         'nkswitchtol':1e-3,
         'adjointl2convergence': 1e-13,
         'frozenTurbulence':False,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25,evalFuncs=evalFuncs)
    CFDSolver = SUMB(options=aeroOptions)
    DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt',complex=True)
    nTwist = 2
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_rans.cgns'})
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    # Aeroproblem must be set before we can call DVGeo.setDesignVars
    CFDSolver.setAeroProblem(ap)
    DVGeo.setDesignVars({'twist':[0,h*1j]})
    CFDSolver(ap,writeSolution=False)
    sol = {}
    CFDSolver.evalFunctions(ap,sol)
    if MPI.COMM_WORLD.rank == 0:
        for key in sorted(sol.keys()):
            print 'sol[%s]:'%key
            reg_write(numpy.real(sol[key]),1e-10,1e-10)
            deriv = numpy.imag(sol[key])/h
            reg_write(deriv,1e-6,1e-6)
    del CFDSolver
    del mesh
    del DVGeo

def test4():
    # ****************************************************************************
    printHeader('MDO tutorial RANS Mesh Alpha Perturbation')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridFile': '../inputFiles/mdo_tutorial_rans.cgns',
         'MGCycle':'2w',
         'equationType':'RANS',
         'smoother':'dadi',
         'CFL':1.5,
         'CFLCoarse':1.25,
         'resaveraging':'noresaveraging',
         'nsubiter':3,
         'nsubiterturb':3,
         'nCyclesCoarse':100,
         'nCycles':1000,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'useNKSolver':True,
         'L2Convergence':1e-15,
         'L2ConvergenceCoarse':1e-4,
         'nkswitchtol':1e-3,
         'adjointl2convergence': 1e-13,
         'frozenTurbulence':False,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8+h*1j, mach=0.80, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25,evalFuncs=evalFuncs)
    CFDSolver = SUMB(options=aeroOptions)

    CFDSolver(ap,writeSolution=False)
    sol = {}
    CFDSolver.evalFunctions(ap,sol)
    if MPI.COMM_WORLD.rank == 0:
        for key in sorted(sol.keys()):
            print 'sol[%s]:'%key
            reg_write(numpy.real(sol[key]),1e-10,1e-10)
            deriv = numpy.imag(sol[key])/h
            reg_write(deriv,1e-6,1e-6)
    del CFDSolver
 
if __name__ == '__main__':
    if len(sys.argv) == 1:
        test1()
        test2()
        test3()
        test4()
    else:
        # Run individual ones
        if 'test1' in sys.argv:
            test1()
        if 'test2' in sys.argv:
            test2()
        if 'test3' in sys.argv:
            test3()
        if 'test4' in sys.argv:
            test4()

