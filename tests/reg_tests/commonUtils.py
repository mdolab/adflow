# This file defines a few common variables that can be used for all
# regression test scripts. This includes default option lists, default
# function lists, and some priting routines.
import numpy

defaultFuncList = ['lift', 'drag', 'cl', 'cd', 'fx', 'fy', 'fz', 'cfx', 'cfy', 'cfz',
                   'mx', 'my', 'mz', 'cmx', 'cmy', 'cmz', 'sepsensor', 'sepsensoravgx',
                   'sepsensoravgy', 'sepsensoravgz']

defaultAeroDVs = ['alpha', 'beta', 'mach', 'P', 'T', 'xRef', 'yRef', 'zRef']

# Note that the option keys here are all consistently LOWERCASE.
adflowDefOpts = {
    # Common Paramters
    'gridfile':'default.cgns',
    'restartfile':None,

    # Output Parameters
    'storerindlayer':True,
    'outputdirectory':'./',
    'writesurfacesolution':True,
    'writevolumesolution':True,
    'nsavevolume':1,
    'nsavesurface':1,
    'solutionprecision':'double',
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
    'blocksplitting':False,
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
    'monitorvariables':['cpu', 'resrho', 'resturb', 'cl', 'cd'],
    'surfacevariables':['cp', 'vx', 'vy', 'vz', 'mach'],
    'volumevariables':[],

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


pyWarpDefOpts = {}

IDWarpDefOpts = {
    'gridFile':None,
    'fileType':'cgns',
    'specifiedSurfaces':None,
    'symmetrySurfaces':None,
    'symmetryPlanes':None,
    'aExp': 3.0,
    'bExp': 5.0,
    'LdefFact':1.0,
    'alpha':0.25,
    'errTol':0.0005,
    'evalMode':'fast',
    'symmTol':1e-6,
    'useRotations':True,
    'zeroCornerRotations':True,
    'cornerAngle':30.0,
    'restartFile':None,
    'bucketSize':8,
}

def solution_test(handler, CFDSolver, ap):

    # Standard test for solving the problem.
    for dv in defaultAeroDVs:
        ap.addDV(dv)

    CFDSolver(ap)

    # Check the residual
    res = CFDSolver.getResidual(ap)
    totalR0, totalRStart, totalRFinal = CFDSolver.getResNorms()
    res /= totalR0

    handler.root_print('Norm of residual')
    handler.par_add_norm(res, 1e-10, 1e-10)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs, defaultFuncList)
    handler.root_print('Eval Functions:')
    handler.root_add_dict(funcs, 1e-10, 1e-10)

    # Get and check the states
    handler.root_print('Norm of state vector')
    handler.par_add_norm(CFDSolver.getStates(), 1e-10, 1e-10)

def adjoint_test(handler, CFDSolver, ap):

    # Standard test for solving multiple adjoints and going right back
    # to the DVs. This solves for whatever functions are in the
    # aeroProblem.
    for dv in defaultAeroDVs:
        ap.addDV(dv)

    res = CFDSolver.getResidual(ap)
    handler.root_print('Norm of residual')
    totalR0 = CFDSolver.getFreeStreamResidual(ap)
    res /= totalR0
    handler.par_add_norm(res, 1e-10, 1e-10)

    funcsSens = {}
    CFDSolver.evalFunctionsSens(ap, funcsSens)
    handler.root_print('Eval Functions Sens:')
    handler.root_add_dict(funcsSens, 1e-10, 1e-10)

def standard_test(handler, CFDSolver, ap, solve):
    # Run a standard set of tests which can be run an any steady grid
    # simulation. This will load in the solution (or solve if
    # solve==True), evaluate the default list of functions, check the
    # forces, and then do a full suite of forward mode, reverse mode
    # and dot-product tests.

    # Now a few simple checks
    handler.root_print('Total number of state DOF')
    handler.par_add_sum(CFDSolver.getStateSize())

    handler.root_print('Total number of adjoint state DOF')
    handler.par_add_sum(CFDSolver.getAdjointStateSize())

    handler.root_print('Total number of spatial DOF')
    handler.par_add_sum(CFDSolver.getSpatialSize())

    for dv in defaultAeroDVs:
        ap.addDV(dv)

    if solve:
        # We are told that we must first solve the problem, most likely
        # for a training run.
        CFDSolver(ap)

    # Check the residual
    res = CFDSolver.getResidual(ap)
    totalR0 = CFDSolver.getFreeStreamResidual(ap)
    res /= totalR0
    handler.root_print('Norm of residual')
    handler.par_add_norm(res, 1e-10, 1e-10)

    funcs = {}
    CFDSolver.evalFunctions(ap, funcs, defaultFuncList)
    handler.root_print('Eval Functions:')
    handler.root_add_dict(funcs, 1e-10, 1e-10)

    # Check the tractions/forces
    forces = CFDSolver.getForces()
    handler.root_print('Sum of Tractions x')
    handler.par_add_sum(forces[:, 0], 1e-10, 1e-10)
    handler.root_print('Sum of Tractions y')
    handler.par_add_sum(forces[:, 1], 1e-10, 1e-10)
    handler.root_print('Sum of Tractions z')
    handler.par_add_sum(forces[:, 2], 1e-10, 1e-10)

    CFDSolver.setOption('forcesAsTractions', False)

    forces = CFDSolver.getForces()
    handler.root_print('Sum of Forces x')
    handler.par_add_sum(forces[:, 0], 1e-10, 1e-10)
    handler.root_print('Sum of Forces y')
    handler.par_add_sum(forces[:, 1], 1e-10, 1e-10)
    handler.root_print('Sum of Forces z')
    handler.par_add_sum(forces[:, 2], 1e-10, 1e-10)

    # Reset the option
    CFDSolver.setOption('forcesAsTractions', True)

    # Make sure we can write the force file.
    CFDSolver.writeForceFile('forces.txt')

    # Get and check the states
    handler.root_print('Norm of state vector')
    states = CFDSolver.getStates()
    handler.par_add_norm(states, 1e-10)

    # Now for the most fun part. Checking the derivatives. These are
    # generally the most important things to check. However, since the
    # checking requires random seeds, it is quite tricky to ensure that we
    # are actually doing the same thing in parallel as in serial. In order
    # to account for this we have to set the random set of numbers to
    # correspond to the full CGNS mesh volume ordering and then scatter
    # back to the correct locations. We have a special routine built into
    # adflow specifically for this purpose.

    # Our "full" jacobian looks like the following:
    #                       residuals  objectives   forces
    #                     +----------+------------+--------+
    #     state Variables |          |            |        |
    #                     +----------+----------- +--------+
    #        volume Nodes |          |            |        |
    #                     +----------+----------- +--------+
    #   "extra" Variables |          |            |        |
    #                     +----------+------------+--------+
    #
    # This defines everything that goes into adflow and everything we care
    # about coming back out. We will check all the derivatives using
    # forward mode AD, reverse mode AD as well as with the dot product
    # test.
    #
    handler.root_print('# ---------------------------------------------------#')
    handler.root_print('#             Forward mode testing                   #')
    handler.root_print('# ---------------------------------------------------#')

    handler.root_print('-> Derivatives with respect to states. wDot, ')
    wDot = CFDSolver.getStatePerturbation(314)

    resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
        wDot=wDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

    handler.root_print('||dR/dw * wDot||')
    handler.par_add_norm(resDot, 1e-10, 1e-10)

    handler.root_print('dFuncs/dw * wDot')
    handler.root_add_dict(funcsDot, 1e-10, 1e-10)

    handler.root_print('||dF/dw * wDot||')
    handler.par_add_norm(fDot, 1e-10, 1e-10)

    handler.root_print('-> Derivatives with respect to nodes')
    xVDot = CFDSolver.getSpatialPerturbation(314)

    resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
        xVDot=xVDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

    handler.root_print('||dR/dXv * xVDot||')
    handler.par_add_norm(resDot, 1e-10, 1e-10)

    # These can be finiky sometimes so a bigger tolerance.
    handler.root_print('dFuncs/dXv * xVDot')
    handler.root_add_dict(funcsDot, 1e-9, 1e-9)

    handler.root_print('||dF/dXv * xVDot||')
    handler.par_add_norm(fDot, 1e-10, 1e-10)

    handler.root_print('-> Derivatives with respect to extra variables')
    for key in defaultAeroDVs:
        handler.root_print('  -> %s'%key)
        xDvDot = {key:1.0}

        resDot, funcsDot, fDot = CFDSolver.computeJacobianVectorProductFwd(
            xDvDot=xDvDot, residualDeriv=True, funcDeriv=True, fDeriv=True)

        handler.root_print('||dR/d%s||'%key)
        handler.par_add_norm(resDot, 1e-10, 1e-10)

        handler.root_print('dFuncs/d%s'%key)
        handler.root_add_dict(funcsDot, 1e-10, 1e-10)

        handler.root_print('||dF/d%s||'%key)
        handler.par_add_norm(fDot, 1e-10, 1e-10)

    handler.root_print('# ---------------------------------------------------#')
    handler.root_print('#             Reverse mode testing                   #')
    handler.root_print('# ---------------------------------------------------#')

    handler.root_print('-> Res bar Seed')
    dwBar = CFDSolver.getStatePerturbation(314)

    wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
        resBar=dwBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

    handler.root_print('||dwBar^T * dR/dw||')
    handler.par_add_norm(wBar, 1e-10, 1e-10)

    handler.root_print('||dwBar^T * dR/dXv||')
    norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
    handler.root_add_val(norm, 1e-10, 1e-10)

    handler.root_print('||dwBar^T * dR/xDv||')
    handler.root_add_dict(xDvBar, 1e-10, 1e-10)

    handler.root_print('-> F Bar Seed')
    fBar = CFDSolver.getSurfacePerturbation(314)

    wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
        fBar=fBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

    handler.root_print('||FBar^T * dF/dw||')
    handler.par_add_norm(wBar, 1e-10, 1e-10)

    handler.root_print('||FBar^T * dF/dXv||')
    norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
    handler.root_add_val(norm, 1e-10, 1e-10)

    handler.root_print('||FBar^T * dF/xDv||')
    handler.root_add_dict(xDvBar, 1e-10, 1e-10)

    handler.root_print(' -> Objective Seeds')

    for key in defaultFuncList:
        handler.root_print('  -> %s'%key)
        funcsBar = {key:1.0}

        wBar, xVBar, xDvBar = CFDSolver.computeJacobianVectorProductBwd(
            funcsBar=funcsBar, wDeriv=True, xVDeriv=True, xDvDerivAero=True)

        handler.root_print('||d%s/dw||'%key)
        handler.par_add_norm(wBar, 1e-10, 1e-10)

        handler.root_print('||d%s/dXv||'%key)
        norm = CFDSolver.getUniqueSpatialPerturbationNorm(xVBar)
        handler.root_add_val(norm, 1e-10, 1e-10)

        handler.root_print('||d%s/dXdv||'%key)
        handler.root_add_dict(xDvBar, 1e-10, 1e-10)

    handler.root_print('# ---------------------------------------------------#')
    handler.root_print('#                 Dot product Tests                  #')
    handler.root_print('# ---------------------------------------------------#')

    # Create a common set of seeds
    wDot = CFDSolver.getStatePerturbation(314)
    xVDot= CFDSolver.getSpatialPerturbation(314)
    dwBar = CFDSolver.getStatePerturbation(314)
    fBar = CFDSolver.getSurfacePerturbation(314)

    handler.root_print ('Dot product test for w -> R')

    dwDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, residualDeriv=True)
    wBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar,  wDeriv=True)

    dotLocal1 = numpy.sum(dwDot*dwBar)
    dotLocal2= numpy.sum(wDot*wBar)

    handler.par_add_sum(dotLocal1, 1e-10, 1e-10)
    handler.par_add_sum(dotLocal2, 1e-10, 1e-10)

    handler.root_print ('Dot product test for Xv -> R')
    dwDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, residualDeriv=True)
    xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar,  xVDeriv=True)

    dotLocal1 = numpy.sum(dwDot*dwBar)
    dotLocal2 = numpy.sum(xVDot*xVBar)

    # For laminar/Rans the DP test for nodes is generally appears a
    # little less accurate. This is ok. It has to do with the very
    # small offwall spacing on laminar/RANS meshes.
    handler.par_add_sum(dotLocal1, 1e-9, 1e-9)
    handler.par_add_sum(dotLocal2, 1e-9, 1e-9)

    handler.root_print ('Dot product test for w -> F')
    wBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar,  wDeriv=True)
    fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, fDeriv=True)

    dotLocal1 = numpy.sum(fDot.flatten()*fBar.flatten())
    dotLocal2 = numpy.sum(wDot*wBar)

    handler.par_add_sum(dotLocal1, 1e-10, 1e-10)
    handler.par_add_sum(dotLocal2, 1e-10, 1e-10)

    handler.root_print ('Dot product test for xV -> F')

    fDot = CFDSolver.computeJacobianVectorProductFwd(xVDot=xVDot, fDeriv=True)
    xVBar = CFDSolver.computeJacobianVectorProductBwd(fBar=fBar,  xVDeriv=True)

    dotLocal1 = numpy.sum(fDot.flatten()*fBar.flatten())
    dotLocal2 = numpy.sum(xVDot*xVBar)

    handler.par_add_sum(dotLocal1, 1e-10, 1e-10)
    handler.par_add_sum(dotLocal2, 1e-10, 1e-10)


    handler.root_print ('Dot product test for (w, xV) -> (dw, F)')

    dwDot, fDot = CFDSolver.computeJacobianVectorProductFwd(wDot=wDot, xVDot=xVDot,
                                                            residualDeriv=True, fDeriv=True)
    wBar, xVBar = CFDSolver.computeJacobianVectorProductBwd(resBar=dwBar, fBar=fBar,
                                                            wDeriv=True, xVDeriv=True)

    dotLocal1 = numpy.sum(dwDot*dwBar) + numpy.sum(fDot.flatten()*fBar.flatten())
    dotLocal2= numpy.sum(wDot*wBar) + numpy.sum(xVDot*xVBar)

    handler.par_add_sum(dotLocal1, 1e-10, 1e-10)
    handler.par_add_sum(dotLocal2, 1e-10, 1e-10)

