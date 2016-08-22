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
from pywarp import MBMesh
import pyspline

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from sumb import SUMB' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pySUmb import SUMB
# ###################################################################

# First thing we will do is define a complete set of default options
# that will be reused as we do differnt tests.  These are the default
# options as Oct 14, 2015.

defOpts = {
    # Common Paramters
    'gridFile':'default.cgns',
    'restartfile':None,

    # Output Parameters
    'storerindlayer':True,
    'outputdirectory':'./',
    'writesymmetry':True,
    'writefarfield':False,
    'writesurfacesolution':True,
    'writevolumesolution':True,
    'nsavevolume':1,
    'nsavesurface':1,
    'solutionprecision':'single',
    'gridprecision':'double',
    'isosurface':{},
    'isovariables':[],
    'viscoussurfacevelocities':True,

    # Physics Paramters
    'discretization':'central plus scalar dissipation',
    'coarsediscretization':'central plus scalar dissipation',
    'limiter':'vanalbeda',
    'smoother':'runge kutta',
    'equationtype':'euler',
    'equationmode':'steady',
    'flowtype':'external',
    'turbulencemodel':'sa',
    'turbulenceorder':'first order',
    'turbresscale':10000.0,
    'usewallfunctions':False,
    'useapproxwalldistance':True,
    'eulerwalltreatment':'linear pressure extrapolation',
    'dissipationscalingexponent':0.67,
    'vis4':0.0156,
    'vis2':0.25,
    'vis2coarse':0.5,
    'restrictionrelaxation':.80,
    'liftindex':2,
    'lowspeedpreconditioner':False,

    # Common Paramters
    'ncycles':500,
    'ncyclescoarse':500,
    'nsubiterturb':1,
    'nsubiter':1,
    'cfl':1.7,
    'cflcoarse':1.0,
    'mgcycle':'3w',
    'mgstartlevel':-1,
    'resaveraging':'alternateresaveraging',
    'smoothparameter':1.5,
    'cfllimit':1.5,

    # Unsteady Paramters
    'timeintegrationscheme':'bdf',
    'timeaccuracy':2,
    'ntimestepscoarse':48,
    'ntimestepsfine':400,
    'deltat':.010,

    # Grid motion parameters
    'useale':True,
    'usegridmotion':False,

    # Time Spectral Paramters
    'timeintervals':1,
    'alphamode':False,
    'betamode':False,
    'machmode':False,
    'pmode':False,
    'qmode':False,
    'rmode':False,
    'altitudemode':False,
    'windaxis':False,
    'alphafollowing':True,
    'tsstability':False,

    # Convergence Paramters
    'l2convergence':1e-6,
    'l2convergencerel':1e-16,
    'l2convergencecoarse':1e-2,
    'maxl2deviationfactor':1.0,

    # Newton-Krylov Paramters
    'usenksolver':False,
    'nkswitchtol':2.5e-4,
    'nksubspacesize':60,
    'nklinearsolvetol':0.3,
    'nkuseew':True,
    'nkadpc':False,
    'nkviscpc':False,
    'nkasmoverlap':1,
    'nkpcilufill':2,
    'nkjacobianlag':20,
    'rkreset':False,
    'nrkreset':5,
    'applypcsubspacesize':10,
    'nkinnerpreconits':1,
    'nkouterpreconits':1,
    'nkls':'cubic',
    'nkcfl0':1000000000000.0,

    # Approximate Newton-Krylov Parameters
    'useanksolver':False,
    'ankswitchtol':1e-2,
    'anksubspacesize':5,
    'anklinearsolvetol':0.5,
    'ankasmoverlap':1,
    'ankpcilufill':1,
    'ankjacobianlag':20,
    'ankinnerpreconits':1,

    # LoadBalance/partitioning Parameters
    'blocksplitting':True,
    'loadimbalance':0.1,
    'loadbalanceiter':10,
    'partitiononly':False,

    # Misc Paramters
    'autosolveretry':False,
    'autoadjointretry':False,
    'numbersolutions':True,
    'printiterations':True,
    'printtiming':True,
    'setmonitor':True,
    'printwarnings':True,
    'monitorvariables':['cpu','resrho','cl','cd'],
    'surfacevariables':['cp','vx','vy','vz','mach'],
    'volumevariables':['resrho'],

    # Multidisciplinary Coupling Parameters:
    'forcesastractions':True,

    # Adjoint Paramters
    'adjointl2convergence':1e-6,
    'adjointl2convergencerel':1e-16,
    'adjointl2convergenceabs':1e-16,
    'adjointdivtol':1e5,
    'approxpc':True,
    'adpc':False,
    'viscpc':False,
    'usediagtspc':True,
    'restartadjoint':True,
    'adjointsolver':'gmres',
    'adjointmaxiter':500,
    'adjointsubspacesize':100,
    'adjointmonitorstep':10,
    'dissipationlumpingparameter':6.0,
    'preconditionerside':'right',
    'matrixordering':'rcm',
    'globalpreconditioner':  'additive schwartz',
    'localpreconditioner':'ilu',
    'ilufill':2,
    'asmoverlap':1,
    'innerpreconits':1,
    'outerpreconits':3,
    'applyadjointpcsubspacesize':20,
    'frozenturbulence':True,
    'usematrixfreedrdw':True,

    # ADjoint debugger
    'firstrun':True,
    'verifystate':True,
    'verifyspatial':True,
    'verifyextra':True,

    # Function Parmeters
    'sepsensoroffset':0.0,
    'sepsensorsharpness':10.0,
}

# First thing we will test is the euler mesh of the MDO tutorial. This
# is a very small mesh that can be run very quickly. We therefore do a
# lot of testing with it.

def printHeader(testName):
    if MPI.COMM_WORLD.rank == 0:
        print '+' + '-'*78 + '+'
        print '| Test Name: ' + '%-66s'%testName + '|'
        print '+' + '-'*78 + '+'

def test1():
    # ****************************************************************************
    printHeader('MDO tutorial Euler Mesh - Python functionality testing')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_euler.cgns',
         'mgcycle':'2w',
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':10000,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'usenksolver':True,
         'l2convergence':1e-14,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-14,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cl','cd','cfx','cfy','cfz',
                                                             'cmx','cmy','cmz'])
    ap.addDV('alpha')
    ap.addDV('mach')
    ap.addDV('altitude')
    CFDSolver = SUMB(options=aeroOptions, debug=True)
    DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
    nTwist = 6
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_euler.cgns'}, debug=False)
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    CFDSolver.addLiftDistribution(10, 'z')
    CFDSolver.addSlices('z', [0.1, 1, 10.0], sliceType='relative')
    CFDSolver.addSlices('z', [5.0], sliceType='absolute')
    surf = CFDSolver.getTriangulatedMeshSurface()
    if MPI.COMM_WORLD.rank == 0:
        print 'Sum of Triangulated Surface:'
        reg_write(numpy.sum(surf[0]))
        reg_write(numpy.sum(surf[1]))
        reg_write(numpy.sum(surf[2]))
  
    CFDSolver(ap)
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-10, 1e-10)

    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions Sens:'
        reg_write_dict(funcsSens, 1e-10, 1e-10)

    # Get the forces...these are the sumb forces:
    forces = CFDSolver.getForces()
    forces = MPI.COMM_WORLD.allreduce(numpy.sum(forces), MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        reg_write(forces,1e-10, 1e-10)
    CFDSolver.setOption('forcesAsTractions', False)

    forces = CFDSolver.getForces()
    forces = MPI.COMM_WORLD.allreduce(numpy.sum(forces), MPI.SUM)
    if MPI.COMM_WORLD.rank == 0:
        reg_write(forces,1e-10, 1e-10)
    CFDSolver.writeForceFile('forces.txt')

    # Now test the different discretization options:
    printHeader('MDO tutorial Euler Mesh - Matrix dissipation')
    CFDSolver.setOption('discretization','central plus matrix dissipation')
    CFDSolver.setOption('coarsediscretization','central plus matrix dissipation')
    CFDSolver.setOption('vis4',0.1)
    CFDSolver.setOption('CFLCoarse',0.75)
    CFDSolver.resetFlow(ap)
    CFDSolver(ap)
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        reg_write_dict(funcs, 1e-10, 1e-10)

    printHeader('MDO tutorial Euler Mesh - Upwind dissipation')
    CFDSolver.setOption('discretization','upwind')
    CFDSolver.setOption('coarseDiscretization','upwind')
    CFDSolver.resetFlow(ap)
    CFDSolver(ap)
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        reg_write_dict(funcs, 1e-10, 1e-10)

    # Test the solve CL Routine
    printHeader('MDO tutorial Euler Mesh - SolveCL Check')
    CFDSolver.resetFlow(ap)
    CFDSolver.solveCL(ap, 0.475, alpha0=0, delta=0.1, tol=1e-4, autoReset=True)
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs, evalFuncs=['cl'])
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'CL-CL*'
        reg_write(funcs['mdo_tutorial_cl'] - 0.475, 1e-4, 1e-4)

    # Clean up:
    del CFDSolver
    del mesh
    del DVGeo
    os.system('rm -fr *.cgns *.dat')

def test2():
    # ****************************************************************************
    printHeader('MDO tutorial Random Euler Mesh')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_euler_random.cgns',
         'mgcycle':'2w',
         'smoother':'dadi',
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':10000,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'usenksolver':True,
         'l2convergence':1e-14,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-14,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cl','cd','cfx','cfy','cfz',
                                                             'cmx','cmy','cmz'])
    ap.addDV('alpha')
    ap.addDV('mach')
    CFDSolver = SUMB(options=aeroOptions)
    DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
    nTwist = 6
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_euler_random.cgns'})
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    CFDSolver(ap)
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-10, 1e-10)

    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions Sens:'
        reg_write_dict(funcsSens, 1e-10, 1e-10)

    # Clean up:
    del CFDSolver
    del mesh
    del DVGeo
    os.system('rm -fr *.cgns *.dat')


def test3():
    # ****************************************************************************
    printHeader('MDO tutorial DADI Euler Mesh')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_euler.cgns',
         'mgcycle':'2w',
         'smoother':'dadi',
         'cfl':10.0,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':10000,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'l2convergence':1e-12,
         'l2convergencecoarse':1e-2,
         'adjointl2convergence': 1e-14,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cl','cd','cfx','cfy','cfz',
                                                             'cmx','cmy','cmz'])
    ap.addDV('alpha')
    ap.addDV('mach')
    CFDSolver = SUMB(options=aeroOptions)
    CFDSolver(ap)
    # Just check the functions
    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-10, 1e-10)
   
    # Clean up:
    del CFDSolver
    os.system('rm -fr *.cgns *.dat')

def test4():
    # ****************************************************************************
    printHeader('MDO tutorial 1-Processor Test')
    # ****************************************************************************
    from petsc4py import PETSc
    if MPI.COMM_WORLD.rank == 0:
        aeroOptions = copy.deepcopy(defOpts)

        # Now set the options that need to be overwritten for this example:
        aeroOptions.update(
            {'gridfile': '../inputFiles/mdo_tutorial_euler.cgns',
             'mgcycle':'2w',
             'cfl':1.5,
             'cflcoarse':1.25,
             'ncyclescoarse':250,
             'ncycles':10000,
             'monitorvariables':['resrho','cl','cd','cmz','totalr'],
             'l2convergence':1e-14,
             'l2convergenceCoarse':1e-2,
             'adjointl2convergence': 1e-14,
             'usenksolver':True,
             'nkswitchtol':1e-2,
         }
        )

        # Setup aeroproblem, cfdsolver, mesh and geometry.
        ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, altitude=10000.0,
                         areaRef=45.5, chordRef=3.25, evalFuncs=['cl', 'drag'])
        ap.addDV('alpha')
        ap.addDV('mach')
        CFDSolver = SUMB(options=aeroOptions, comm=MPI.COMM_SELF)
        DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
        nTwist = 6
        DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                                y=numpy.zeros(nTwist),
                                                z=numpy.linspace(0,14, nTwist), k=2))
        def twist(val, geo):
            for i in xrange(nTwist):
                geo.rot_z['wing'].coef[i] = val[i]

        DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
        DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
        mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_euler.cgns'}, comm=MPI.COMM_SELF)
        CFDSolver.setMesh(mesh)
        CFDSolver.setDVGeo(DVGeo)
        CFDSolver(ap)
        print 'Eval Functions:'
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        CFDSolver.checkSolutionFailure(ap, funcs)
        reg_write_dict(funcs, 1e-10, 1e-10)
        print 'Eval Functions Sens:'
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)
        reg_write_dict(funcsSens, 1e-10, 1e-10)

        # Clean up:
        del CFDSolver
        del mesh
        del DVGeo
        os.system('rm -fr *.cgns *.dat')

def test5():
    # THIS TEST NEEDS TO BE VERIFIED WITH CS AND IT IS NEEDS TO BE READDED TO REGRESSIONS
    # ****************************************************************************
    printHeader('MDO tutorial Euler Time Spectral Test')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_euler_l2.cgns',
         'mgcycle':'sg',
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':10000,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'l2convergence':1e-14,
         'l2convergencecoarse':1e-2,
         'adjointl2convergence': 1e-14,
         'timetntervals':3,
         'equationmode':'Time Spectral',
         'tsstability':True,
         'alphamode':True,
         'usenksolver':True,
         'nkswitchtol':1e-2,
         'usediagtspc':False,
         'asmoverlap':2,
         'ilufill':3,
         'adjointmaxiter':1000,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25, 
                     degreeFourier=1, omegaFourier=6.28, degreePol=0,
                     cosCoefFourier=[0.0,0.0], sinCoefFourier=[0.01], 
                     coefPol=[0],
                     evalFuncs=['cl', 'cl0','clalpha', 'clalphadot'])
    ap.addDV('alpha')
    # Note that the mach number derivative for the timspectral is
    # broken and needs to be fixed. 
    CFDSolver = SUMB(options=aeroOptions)
    DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
    nTwist = 6
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_euler_l2.cgns'})
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    CFDSolver(ap)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-10, 1e-10)

    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions Sens:'
        reg_write_dict(funcsSens, 1e-10, 1e-10)

    # Clean up:
    del CFDSolver
    del mesh
    del DVGeo
    os.system('rm -fr *.cgns* *.dat')
    
def test6():
    # ****************************************************************************
    printHeader('MDO tutorial Viscous Mesh')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_rans.cgns',
         'mgcycle':'2w',
         'equationtype':'Laminar NS',
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':10000,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-15,
         'l2Convergencecoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-14,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.50, 
                     reynolds=50000.0, reynoldsLength=3.25, T=293.15,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cd','lift'])
    ap.addDV('alpha')
    ap.addDV('mach')
    CFDSolver = SUMB(options=aeroOptions)
    DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
    nTwist = 6
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
    CFDSolver(ap)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-10, 1e-10)

    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions Sens:'
        reg_write_dict(funcsSens, 1e-10, 1e-10)

    # Clean up:
    del CFDSolver
    del mesh
    del DVGeo
    os.system('rm -fr *.cgns *.dat')

def test7():
    # ****************************************************************************
    printHeader('MDO tutorial RANS Mesh')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_rans.cgns',
         'mgcycle':'2w',
         'equationtype':'RANS',
         'smoother':'dadi',
         'cfl':1.5,
         'cflcoarse':1.25,
         'resaveraging':'noresaveraging',
         'nsubiter':3,
         'nsubiterturb':3,
         'ncyclescoarse':100,
         'ncycles':1000,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-14,
         'l2convergencecoarse':1e-4,
         'nkswitchtol':1e-3,
         'adjointl2convergence': 1e-14,
         'frozenturbulence':False,
     }
    )
    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, 
                     reynolds=50000.0, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cd','lift','cmz'])
    ap.addDV('alpha')
    ap.addDV('mach')
    CFDSolver = SUMB(options=aeroOptions)
    DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
    nTwist = 6
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
    CFDSolver(ap)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-9, 1e-9)
    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions Sens:'
        reg_write_dict(funcsSens, 1e-9, 1e-9)

    # Clean up:
    del CFDSolver
    del mesh
    del DVGeo
    os.system('rm -fr *.cgns *.dat')

def test8():
    # ****************************************************************************
    printHeader('MDO tutorial Random RANS mesh')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_rans_random.cgns',
         'mgcycle':'2w',
         'equationtype':'RANS',
         'smoother':'dadi',
         'cfl':1.5,
         'cflcoarse':1.25,
         'resaveraging':'noresaveraging',
         'nsubiter':3,
         'nsubiterturb':3,
         'ncyclescoarse':100,
         'ncycles':1000,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-14,
         'l2convergencecoarse':1e-4,
         'nkswitchtol':1e-3,
         'adjointl2convergence': 1e-14,
         'frozenturbulence':False,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, 
                     reynolds=50000.0, altitude=10000.0,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cd','lift','cmz'])
    ap.addDV('alpha')
    ap.addDV('mach')
    CFDSolver = SUMB(options=aeroOptions)
    DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
    nTwist = 6
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_rans_random.cgns'})
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    CFDSolver(ap)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-9, 1e-9)
    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions Sens:'
        reg_write_dict(funcsSens, 1e-9, 1e-9)

    # Clean up:
    del CFDSolver
    del mesh
    del DVGeo
    os.system('rm -fr *.cgns *.dat')

def test9():
    # ****************************************************************************
    printHeader('CRM WBT Euler Mesh')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/dpw4_38k.cgns',
         'mgcycle':'sg',
         'cfl':1.5,
         'cflcoarse':1.25,
         'resaveraging':'noresaveraging',
         'ncycles':1000,
         'monitorvariables':['resrho','cl','cd','cmy','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-14,
         'l2convergencecoarse':1e-4,
         'nkswitchtol':1e-1,
         'adjointl2convergence': 1e-14,
         'liftindex':3,
     }
    )

    # Setup aeroproblem, cfdsolver, mesh and geometry.
    ap = AeroProblem(name='crm', alpha=1.0, mach=0.85, 
                     reynolds=50000.0, altitude=8000.0,
                     areaRef=594720*.0254**2/2.0, chordRef=275.8*.0254, 
                     evalFuncs=['cd','lift','cmy'])
    ap.addDV('alpha')
    ap.addDV('mach')
    CFDSolver = SUMB(options=aeroOptions)
    DVGeo = DVGeometry('../inputFiles/CRM_ffd.fmt')

    # Setup curves for ref_axis
    leRoot = numpy.array([25.22, 3.08, 4.46])
    leTip = numpy.array([45.1735, 29.4681, 4.91902])
    rootChord = 11.83165
    breakChord = 7.25894
    tipChord = 2.7275295

    # We have to be careful with the reference axis...We need to ensure
    # that there is a point on the ref axis *right* at the symmetry plane.
    # X1 information is taken directly from the CRM section data table in
    # Vassberg.
    
    X1 = numpy.array([22.97+.25*13.62, 0, 4.42])
    X2 = leRoot + numpy.array([1.0, 0.0, 0.0])*rootChord*.25
    X3 = leTip + numpy.array([1.0, 0.0, 0.0])*tipChord*.25
    X = []
    X.append(X1)
    nTwist = 8
    for i in range(nTwist):
        fact = float(i)/(nTwist-1)
        X.append((1.0-fact)*X2 + fact*X3)

    c1 = pyspline.Curve(X=X, k=2)
    DVGeo.addRefAxis('wing', c1, volumes=[0,5])

    x = numpy.array([2365.0 , 2365.0])*.0254
    y = numpy.array([0, 840/2.0])*.0254
    z = numpy.array([255.0, 255.0])*.0254
    c2 = pyspline.Curve(x=x, y=y, z=z, k=2)
    DVGeo.addRefAxis('tail', c2, volumes=[25])

    def twist(val, geo):
        # Set all the twist values
        for i in xrange(nTwist-1):
            geo.rot_y['wing'].coef[i+2] = val[i]

    def tailTwist(val, geo):
        # Set one twist angle for the tail
        geo.rot_y['tail'].coef[:] = val[0]

    def ssd(val, geo):
        # Span-sweep-dihedreal --- move the tip in any direction
        C = geo.extractCoef('wing')
        s = geo.extractS('wing')
        # Get a second s that only scaled from the 2nd pt out:
        ss = (s[1:] - s[1])/(s[-1] - s[1])
        for i in xrange(len(C)-1):
            C[i+1, 0] = C[i+1, 0] + ss[i]*val[0]
            C[i+1, 1] = C[i+1, 1] + ss[i]*val[1]
            C[i+1, 2] = C[i+1, 2] + ss[i]*val[2]
        geo.restoreCoef(C, 'wing')

    def chord(val, geo):
        geo.scale['wing'].coef[:] = val[0]

    DVGeo.addGeoDVGlobal('twist', numpy.zeros(nTwist-1), twist, lower=-10, upper=10)
    DVGeo.addGeoDVGlobal('ssd', [0,0,0], ssd, lower=-20, upper=20)
    DVGeo.addGeoDVGlobal('chord', [1.0], chord, lower=0.75, upper=1.25)
    DVGeo.addGeoDVLocal('shape', lower=-1.0, upper=1.0, axis='z', scale=1.0, 
                        volList=[0])
    meshOptions ={
        'gridFile':'../inputFiles/dpw4_38k.cgns',
        'warpType':'solid',
        'solidWarpType':'n',
        'n':3, 
        'fillType':'linear'}

    mesh = MBMesh(options=meshOptions)
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    CFDSolver(ap)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-10, 1e-10)
    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions Sens:'
        reg_write_dict(funcsSens, 1e-10, 1e-10)

    # Clean up:
    del CFDSolver
    del mesh
    del DVGeo
    os.system('rm -fr *.cgns *.dat')

def test10():
    # ****************************************************************************
    printHeader('NACA 0012 2D Time-Accurate, Forced motion, Rigid Rotation of Mesh - DADI Smoother')
    # ****************************************************************************
    
    # THIS TEST NEEDS TO BE REFINED. NEED TO
    # * check values at each timestep
    # * drive the residual for each timestep further down using the newton solver 
    #   and compare more digits
    # * add better mesh
    
    aeroOptions = copy.deepcopy(defOpts)

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

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
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
     }
    )

    ap = AeroProblem(name=name, alpha=alpha_m,  mach=M, machRef=M, reynolds=4800000.0,reynoldsLength=c, T=T, R=R,
                     areaRef=1.0, chordRef=c, evalFuncs=['cl','cd','cmz'],xRef=0.25,xRot=0.25,
                     degreePol=0,coefPol=[0.0],degreeFourier=1,omegaFourier=omega,
                     cosCoefFourier=[0.0,0.0],sinCoefFourier=[deltaAlpha])

    CFDSolver = SUMB(options=aeroOptions)
    CFDSolver.addSlices('z',[0.5])
    CFDSolver(ap)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-6, 1e-6)

    # Clean up:
    del CFDSolver
    os.system('rm -fr *.cgns *.dat *.cgns*')
        

def test11():
    # ****************************************************************************
    printHeader('BSCW 3D Time-Accurate Rigid Rotation - Case 1 of AePW-2')
    # ****************************************************************************

    # DISABLED FOR NOW SINCE THIS TAKES TO LONG TO RUN.

    # Set forcing frequency and other information
    f = 10.0 # [Hz] Forcing frequency of the flow
    period = 1.0/f # [sec]
    nStepPerPeriod = 64
    nPeriods = 4
    nfineSteps = nStepPerPeriod*nPeriods
    dt = period / nStepPerPeriod # [s] The actual timestep

    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/bscw_rans_finer_L3.cgns',
         'liftindex':3,
         'writevolumesolution':False,
         'usenksolver':False,
         'equationtype':'RANS',
         'equationmode':'unsteady',
         'timeintegrationscheme':'bdf',
         'ntimestepsfine':nfineSteps,
         'deltat':dt,
         'smoother':'dadi',
         'nsubiterturb':3, 
         'nsubiter':5, 
         'useale':False,
         'usegridmotion':True,
         'cfl':1.5,
         'cflcoarse':1.25,
         'mgcycle':'2w',
         'mgstartlevel':1,
         'ncycles': 30, 
         'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cdp','cdv','cmy','yplus','totalr'],
         'l2convergence':1e-8,
         'rmode':True,  
         'alphaFollowing': False,
     }
    )

    meshOptions = {
        'gridFile':'../inputFiles/bscw_rans_finer_L3.cgns',
        }

    ap = AeroProblem(name='bscw', alpha=3.0, mach=0.7,
                     areaRef=512*0.0254**2, chordRef=16*0.0254, 
                     reynoldsLength=16*0.0254, reynolds=4.56e6,
                     evalFuncs=['cl','cd','cmy'],
                     gamma=1.113, Pr=0.683, T=302.9788888889,
                     R=84.3517245542, SSuthDim=243.3722222222, muSuthDim=1.12E-005, TSuthDim=273.0,
                     xRef=0.3*16*0.0254, yRef=0.0, zRef=0.0,
                     xRot=0.3*16*0.0254, yRot=0.0, zRot=0.0,
                     degreePol=0, coefPol=[0.0], degreeFourier=1, omegaFourier=2*numpy.pi*f,
                     cosCoefFourier=[0.0,0.0], sinCoefFourier=[-1.0*numpy.pi/180.0])
                     
    # Generate a mesh object
    Mesh = MBMesh(options=meshOptions)
                     
    # Create solver
    CFDSolver = SUMB(options=aeroOptions)
    CFDSolver.setMesh(Mesh)
    CFDSolver.addFamilyGroup("upperSurface",["upperSurface"])
    CFDSolver.addFamilyGroup("lowerSurface",["lowerSurface"])
    CFDSolver.addSlices('y', [.48768, .77216], groupName="upperSurface")
    CFDSolver.addSlices('y', [.48768, .77216], groupName="lowerSurface")

    # Solve
    CFDSolver(ap)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs)
    CFDSolver.checkSolutionFailure(ap, funcs)
    if MPI.COMM_WORLD.rank == 0:
        print 'Eval Functions:'
        reg_write_dict(funcs, 1e-10, 1e-10)    

    # Clean up:
    del CFDSolver
    del mesh
    os.system('rm -fr *.cgns *.dat')


if __name__ == '__main__':
    if len(sys.argv) == 1:
        test1()
        test2()
        test3()
        test4()
        test6()
        test7()
        test8()
        test9()
        test10()
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
        if 'test5' in sys.argv:
            test5()
        if 'test6' in sys.argv:
            test6()
        if 'test7' in sys.argv:
            test7()
        if 'test8' in sys.argv:
            test8()
        if 'test9' in sys.argv:
            test9()
        if 'test10' in sys.argv:
            test10()            
