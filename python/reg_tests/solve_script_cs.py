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
import pyspline

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from sumb import SUMB' for regular scripts.
sys.path.append(os.path.abspath('../../'))
if 'complex' in sys.argv:
    from python.pySUmb_C import SUMB_C as SUMB
    from pywarp import MBMesh_C as MBMesh
else:
    from pywarp import MBMesh
    from python.pySUmb import SUMB


# ###################################################################

# First thing we will do is define a complete set of default options
# that will be reused as we do differnt tests.  These are the default
# options as Dec 8, 2014.

defOpts = {
    # Common Paramters
    'gridfile': 'default.cgns',
    'restartfile':None,

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
    'eulerwalltreatment': 'linear pressure extrapolation',
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

    # Newton-Krylov Paramters
    'usenksolver': False,
    'nkswitchtol': 2.5e-4,
    'nksubspacesize': 60,
    'nklinearsolvetol': 0.3,
    'nkuseew': True,
    'nkadpc': False,
    'nkviscpc': False,
    'nkasmoverlap': 1,
    'nkpcilufill': 2,
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
    'autosolveretry': False,
    'autoadjointretry': False,
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
    'applyadjointpcsubspacesize': 20,
    'frozenturbulence': True,
    'usematrixfreedrdw': False,

    # ADjoint debugger
    'firstrun': True,
    'verifystate': True,
    'verifyspatial': True,
    'verifyextra': True,
}

def printHeader(testName):
    if MPI.COMM_WORLD.rank == 0:
        print '+' + '-'*78 + '+'
        print '| Test Name: ' + '%-66s'%testName + '|'
        print '+' + '-'*78 + '+'
h = 1e-40



def test1():
    # ****************************************************************************
    printHeader('MDO tutorial Euler Aerodynamic Variables')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_euler.cgns',
         'mgcycle':'2w',
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':400,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'usenksolver':True,
         'l2convergence':1e-15,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-15,
         'nkls':'non monotone',
     }
    )

    # Setup aeroproblem, cfdsolver
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8+h*1j, mach=0.80,
                     altitude=10000.0, areaRef=45.5, chordRef=3.25,
                     evalFuncs=['cd','lift','cmz'])
    ap.addDV('alpha')
    ap.addDV('mach')
    ap.addDV('altitude')
    CFDSolver = SUMB(options=aeroOptions, debug=True)

    if not 'complex' in sys.argv:
        # Solve system
        CFDSolver(ap, writeSolution=False)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        # Solve sensitivities
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        # Write values and derivatives out:
        if MPI.COMM_WORLD.rank == 0:
            for key in ['cd','cmz','lift']:
                print 'funcs[%s]:'%key
                reg_write(funcs['mdo_tutorial_%s'%key],1e-10,1e-10)
            # Now write the derivatives in the same order the CS will do them:
            print ('Alpha Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['alpha_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['alpha_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_lift']['alpha_mdo_tutorial'], 1e-10,1e-10)

            print ('Mach Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_lift']['mach_mdo_tutorial'], 1e-10,1e-10)

            print ('Altitude Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['altitude_mdo_tutorial'], 1e-8,1e-8)
            reg_write(funcsSens['mdo_tutorial_cmz']['altitude_mdo_tutorial'], 1e-8,1e-8)
            reg_write(funcsSens['mdo_tutorial_lift']['altitude_mdo_tutorial'], 1e-8,1e-8)

    else:
        # For the complex....we just do successive perturbation
        for ii in range(3):
            ap.alpha = 1.8
            ap.mach = 0.80
            ap.altitude = 10000.0
            if ii == 0:
                ap.alpha += h*1j
            elif ii == 1:
                ap.mach += h*1j
            else:
                ap.altitude += h*1j

            CFDSolver.resetFlow(ap)
            CFDSolver(ap, writeSolution=False)
            funcs = {}
            CFDSolver.evalFunctions(ap, funcs)

            if MPI.COMM_WORLD.rank == 0:
                if ii == 0:
                    for key in ['cd','cmz','lift']:
                        print 'funcs[%s]:'%key
                        reg_write(numpy.real(funcs['mdo_tutorial_%s'%key]),1e-10,1e-10)

                if ii == 0:
                    print ('Alpha Derivatives:')
                    for key in ['cd','cmz','lift']:
                        deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                        reg_write(deriv,1e-10,1e-10)

                elif ii == 1:
                    print ('Mach Derivatives:')
                    for key in ['cd','cmz','lift']:
                        deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                        reg_write(deriv,1e-10,1e-10)

                elif ii == 2:
                    print ('Altitude Derivatives:')
                    for key in ['cd','cmz','lift']:
                        deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                        reg_write(deriv,1e-8,1e-8)



    del CFDSolver

def test2():
    # ****************************************************************************
    printHeader('MDO tutorial Euler Geometric Variables')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_euler.cgns',
         'mgcycle':'2w',
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':400,
         'monitorvariables':['resrho','cl','cd','cmz','totalr'],
         'usenksolver':True,
         'l2convergence':1e-15,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-15,
         'nkls':'non monotone',
     }
    )

    # Setup aeroproblem, cfdsolver
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80,
                     altitude=10000.0, areaRef=45.5, chordRef=3.25,
                     evalFuncs=['cl','cmz','drag'])

    CFDSolver = SUMB(options=aeroOptions)
    if 'complex' in sys.argv:
        DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt', complex=True)
    else:
        DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt', complex=False)

    nTwist = 2
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    def span(val, geo):
        # Span
        C = geo.extractCoef('wing')
        s = geo.extractS('wing')
        for i in xrange(len(C)-1):
            C[-1, 2] = C[-1, 2] + val[0]
        geo.restoreCoef(C, 'wing')

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVGlobal('span', [0], span, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_euler.cgns'})
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    #Aeroproblem must be set before we can call DVGeo.setDesignVars
    CFDSolver.setAeroProblem(ap)
    if not 'complex' in sys.argv:
        # Solve system
        CFDSolver(ap, writeSolution=False)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        # Solve sensitivities
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        # Write values and derivatives out:
        if MPI.COMM_WORLD.rank == 0:
            for key in ['cl','cmz','drag']:
                print 'funcs[%s]:'%key
                reg_write(funcs['mdo_tutorial_%s'%key],1e-10,1e-10)
            # Now write the derivatives in the same order the CS will do them:
            print ('Twist[0] Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['twist'][0][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['twist'][0][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['twist'][0][0], 1e-10,1e-10)

            print ('Span Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['span'][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['span'][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['span'][0], 1e-10,1e-10)

            print ('shape[13] Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['shape'][0][13], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['shape'][0][13], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['shape'][0][13], 1e-10,1e-10)
    else:
        # For the complex....we just do successive perturbation
        for ii in range(3):
            xRef = {'twist':[0.0, 0.0], 'span':[0.0], 'shape':numpy.zeros(72, dtype='D')}
            if ii == 0:
                xRef['twist'][0] += h*1j
            elif ii == 1:      
                xRef['span'][0] += h*1j
            else:
                xRef['shape'][13] += h*1j

            CFDSolver.resetFlow(ap)
            DVGeo.setDesignVars(xRef)
            CFDSolver(ap, writeSolution=False)
            funcs = {}
            CFDSolver.evalFunctions(ap, funcs)

            if MPI.COMM_WORLD.rank == 0:
                if ii == 0:
                    for key in ['cl','cmz','drag']:
                        print 'funcs[%s]:'%key
                        reg_write(numpy.real(funcs['mdo_tutorial_%s'%key]),1e-10,1e-10)

                if ii == 0:
                    print ('Twist[0] Derivatives:')
                elif ii == 1:
                    print ('Span Derivatives:')
                elif ii == 2:
                    print ('shape[13] Derivatives:')

                for key in ['cl','cmz','drag']:
                    deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                    reg_write(deriv,1e-10,1e-10)
    del CFDSolver
    del mesh

# Notes for the Viscous Derivatives:

def test3():
    # ****************************************************************************
    printHeader('MDO tutorial Viscous Aerodynamic Variables')
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
         'ncycles':500,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-16,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-15,
         'nkls':'non monotone',
     }
    )

    # Setup aeroproblem, cfdsolver
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.50, 
                     reynolds=50000.0, reynoldsLength=3.25, T=293.15,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cd','cmz','lift'])
    ap.addDV('alpha')
    ap.addDV('mach')

    CFDSolver = SUMB(options=aeroOptions)

    if not 'complex' in sys.argv:
        # Solve system
        CFDSolver(ap, writeSolution=False)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        # Solve sensitivities
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        # Write values and derivatives out:
        if MPI.COMM_WORLD.rank == 0:
            for key in ['cd','cmz','lift']:
                print 'funcs[%s]:'%key
                reg_write(funcs['mdo_tutorial_%s'%key],1e-10,1e-10)
            # Now write the derivatives in the same order the CS will do them:
            print ('Alpha Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['alpha_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['alpha_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_lift']['alpha_mdo_tutorial'], 1e-10,1e-10)

            print ('Mach Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_lift']['mach_mdo_tutorial'], 1e-10,1e-10)

    else:
        # For the complex....we just do successive perturbation
        for ii in range(2):
            ap.alpha = 1.8
            ap.mach = 0.50

            if ii == 0:
                ap.alpha += h*1j
            elif ii == 1:
                ap.mach += h*1j

            CFDSolver.resetFlow(ap)
            CFDSolver(ap, writeSolution=False)
            funcs = {}
            CFDSolver.evalFunctions(ap, funcs)

            if MPI.COMM_WORLD.rank == 0:
                if ii == 0:
                    for key in ['cd','cmz','lift']:
                        print 'funcs[%s]:'%key
                        reg_write(numpy.real(funcs['mdo_tutorial_%s'%key]),1e-10,1e-10)

                if ii == 0:
                    print ('Alpha Derivatives:')
                elif ii == 1:
                    print ('Mach Derivatives:')

                for key in ['cd','cmz','lift']:
                    deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                    reg_write(deriv,1e-10,1e-10)

    del CFDSolver

def test4():
    # ****************************************************************************
    printHeader('MDO tutorial Viscous Geometric Variables')
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
         'ncycles':500,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-16,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-2,
         'adjointl2convergence': 1e-15,
         'nkls':'non monotone',
     }
    )

    # Setup aeroproblem, cfdsolver
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.50, 
                     reynolds=50000.0, reynoldsLength=3.25, T=293.15,
                     areaRef=45.5, chordRef=3.25, evalFuncs=['cl','cmz','drag'])

    CFDSolver = SUMB(options=aeroOptions)
    if 'complex' in sys.argv:
        DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt', complex=True)
    else:
        DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt', complex=False)

    nTwist = 2
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    def span(val, geo):
        # Span
        C = geo.extractCoef('wing')
        s = geo.extractS('wing')
        for i in xrange(len(C)-1):
            C[-1, 2] = C[-1, 2] + val[0]
        geo.restoreCoef(C, 'wing')

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVGlobal('span', [0], span, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_rans.cgns'})
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    #Aeroproblem must be set before we can call DVGeo.setDesignVars
    CFDSolver.setAeroProblem(ap)
    if not 'complex' in sys.argv:
        # Solve system
        CFDSolver(ap, writeSolution=False)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        # Solve sensitivities
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        # Write values and derivatives out:
        if MPI.COMM_WORLD.rank == 0:
            for key in ['cl','cmz','drag']:
                print 'funcs[%s]:'%key
                reg_write(funcs['mdo_tutorial_%s'%key],1e-10,1e-10)
            # Now write the derivatives in the same order the CS will do them:
            print ('Twist[0] Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['twist'][0][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['twist'][0][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['twist'][0][0], 1e-10,1e-10)

            print ('Span Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['span'][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['span'][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['span'][0], 1e-10,1e-10)

            print ('shape[13] Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['shape'][0][13], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['shape'][0][13], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['shape'][0][13], 1e-10,1e-10)
    else:
        # For the complex....we just do successive perturbation
        for ii in range(3):
            xRef = {'twist':[0.0, 0.0], 'span':[0.0], 'shape':numpy.zeros(72, dtype='D')}
            if ii == 0:
                xRef['twist'][0] += h*1j
            elif ii == 1:      
                xRef['span'][0] += h*1j
            else:
                xRef['shape'][13] += h*1j

            CFDSolver.resetFlow(ap)
            DVGeo.setDesignVars(xRef)
            CFDSolver(ap, writeSolution=False)
            funcs = {}
            CFDSolver.evalFunctions(ap, funcs)

            if MPI.COMM_WORLD.rank == 0:
                if ii == 0:
                    for key in ['cl','cmz','drag']:
                        print 'funcs[%s]:'%key
                        reg_write(numpy.real(funcs['mdo_tutorial_%s'%key]),1e-10,1e-10)

                if ii == 0:
                    print ('Twist[0] Derivatives:')
                elif ii == 1:
                    print ('Span Derivatives:')
                elif ii == 2:
                    print ('shape[13] Derivatives:')

                for key in ['cl','cmz','drag']:
                    deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                    reg_write(deriv,1e-10,1e-10)

    del CFDSolver
    del mesh

# ----------- Most of this has been supersceeded. Derivatives are now
# ----------- accurate to 1e-10 across multiple compilers (gfortran/intel)

# Notes for the RANS Derivatives: It appears that the CS derviative
# itself flucuates in the 10 to 11 digit or so even when the residual
# is bouncing around machine precision. This si the limit of what can
# be verified wtih CS. This is the reason for putting the check
# tolerance at 1e-8. However, the altidue derivatives seem to be quite
# a bit less accurate, but this thought to be caused by the fact that
# they so small in magniude. These are only requried to match to 1e-4.
# Note that the NKJacobian lag is set to 5 which is much lower than is
# norally used. The reason for this is with the RANS analysis,
# sometimes the CS derivative will "blow up -- increase by 6 or 7
# orders of magnitude during the solution. It eventually comes back to
# the right order of magnitude, but it now quite inaccurate. It can
# also happen if the the nkswtich tol is too large. It is now 1e4
# which seems to work ok. Futhermore, the intel compiler does a
# *TERRIBLE* job with accuracy with complex variables. It really is
# awful. When you are playing around with verifying derivatives it
# really is important to use -fp-model precise. For the verification,
# we can only reliably get about 8 digits since the actual values
# won't match much better than that. You should be able to get 10
# digits+ if you are using the precise fp model. Also, the reason for
# doing the aerodynamic and geometric derivatives separately is that
# the when the nodes are embedded in the FFD there may be slight
# differences which will result in differences on the order of
# ~1e-12. This way, the aerodynamic derivatives can be tested
# independent of the goemetric ones.


def test5():
    # ****************************************************************************
    printHeader('MDO tutorial RANS Aerodynamic Variables')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_rans.cgns',
         'mgcycle':'2w',
         'equationtype':'RANS',
         'smoother':'dadi',
         'nsubiterturb':3,
         'nsubiter':3,
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':750,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-17,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-4,
         'adjointl2convergence': 1e-16,
         'nkls': 'non monotone',
         'frozenturbulence':False,
         'nkjacobianlag':2,
     }
    )

    # Setup aeroproblem, cfdsolver
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, 
                     altitude=10000.0, areaRef=45.5, chordRef=3.25, evalFuncs=['cd','cmz','lift'])
    ap.addDV('alpha')
    ap.addDV('mach')
    ap.addDV('altitude')
    CFDSolver = SUMB(options=aeroOptions,debug=True)

    if not 'complex' in sys.argv:
        # Solve system
        CFDSolver(ap, writeSolution=False)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        # Solve sensitivities
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        # Write values and derivatives out:
        if MPI.COMM_WORLD.rank == 0:
            for key in ['cd','cmz','lift']:
                print 'funcs[%s]:'%key
                reg_write(funcs['mdo_tutorial_%s'%key],1e-10,1e-10)
            # Now write the derivatives in the same order the CS will do them:
            print ('Alpha Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['alpha_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['alpha_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_lift']['alpha_mdo_tutorial'], 1e-10,1e-10)

            print ('Mach Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_lift']['mach_mdo_tutorial'], 1e-10,1e-10)

            print ('Altitude Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cd']['altitude_mdo_tutorial'], 1e-8,1e-8)
            reg_write(funcsSens['mdo_tutorial_cmz']['altitude_mdo_tutorial'], 1e-8,1e-8)
            reg_write(funcsSens['mdo_tutorial_lift']['altitude_mdo_tutorial'], 1e-8,1e-8)

    else:
        # For the complex....we just do successive perturbation
        for ii in range(3):
            ap.alpha = 1.8
            ap.mach = 0.80
            ap.altitude = 10000.0
            if ii == 0:
                ap.alpha += h*1j
            elif ii == 1:
                ap.mach += h*1j
            else:
                ap.altitude += h*1j

            CFDSolver.resetFlow(ap)
            CFDSolver(ap, writeSolution=False)
            funcs = {}
            CFDSolver.evalFunctions(ap, funcs)

            if MPI.COMM_WORLD.rank == 0:
                if ii == 0:
                    for key in ['cd','cmz','lift']:
                        print 'funcs[%s]:'%key
                        reg_write(numpy.real(funcs['mdo_tutorial_%s'%key]),1e-10,1e-10)

                if ii == 0:
                    print ('Alpha Derivatives:')
                    for key in ['cd','cmz','lift']:
                        deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                        reg_write(deriv,1e-10,1e-10)

                elif ii == 1:
                    print ('Mach Derivatives:')
                    for key in ['cd','cmz','lift']:
                        deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                        reg_write(deriv,1e-10,1e-10)

                else:
                    print ('AltitudeDerivatives:')
                    for key in ['cd','cmz','lift']:
                        deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                        reg_write(deriv,1e-10,1e-10)


    del CFDSolver

def test6():
    # ****************************************************************************
    printHeader('MDO tutorial RANS Geometric Variables')
    # ****************************************************************************
    aeroOptions = copy.deepcopy(defOpts)

    # Now set the options that need to be overwritten for this example:
    aeroOptions.update(
        {'gridfile': '../inputFiles/mdo_tutorial_rans.cgns',
         'mgcycle':'2w',
         'equationtype':'RANS',
         'smoother':'dadi',
         'nsubiterturb':3,
         'nsubiter':3,
         'cfl':1.5,
         'cflcoarse':1.25,
         'ncyclescoarse':250,
         'ncycles':750,
         'monitorvariables':['resrho','resturb','cl','cd','cmz','yplus','totalr'],
         'usenksolver':True,
         'l2convergence':1e-17,
         'l2convergencecoarse':1e-2,
         'nkswitchtol':1e-4,
         'adjointl2convergence': 1e-16,
         'nkls': 'non monotone',
         'frozenturbulence':False,
         'nkjacobianlag':2,
     }
    )

    # Setup aeroproblem, cfdsolver
    ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, 
                     altitude=10000.0, areaRef=45.5, chordRef=3.25, evalFuncs=['cl','cmz','drag'])
    
    ap.addDV('alpha')
    ap.addDV('mach')
    CFDSolver = SUMB(options=aeroOptions)
    if 'complex' in sys.argv:
        DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt', complex=True)
    else:
        DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt', complex=False)

    nTwist = 2
    DVGeo.addRefAxis('wing', pyspline.Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                                            y=numpy.zeros(nTwist),
                                            z=numpy.linspace(0,14, nTwist), k=2))
    def twist(val, geo):
        for i in xrange(nTwist):
            geo.rot_z['wing'].coef[i] = val[i]

    def span(val, geo):
        # Span
        C = geo.extractCoef('wing')
        s = geo.extractS('wing')
        for i in xrange(len(C)-1):
            C[-1, 2] = C[-1, 2] + val[0]
        geo.restoreCoef(C, 'wing')

    DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVGlobal('span', [0], span, lower=-10, upper=10, scale=1.0)
    DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
    mesh = MBMesh(options={'gridFile':'../inputFiles/mdo_tutorial_rans.cgns'})
    CFDSolver.setMesh(mesh)
    CFDSolver.setDVGeo(DVGeo)
    #Aeroproblem must be set before we can call DVGeo.setDesignVars
    CFDSolver.setAeroProblem(ap)
    if not 'complex' in sys.argv:
        # Solve system
        CFDSolver(ap, writeSolution=False)
        funcs = {}
        CFDSolver.evalFunctions(ap, funcs)
        # Solve sensitivities
        funcsSens = {}
        CFDSolver.evalFunctionsSens(ap, funcsSens)

        # Write values and derivatives out:
        if MPI.COMM_WORLD.rank == 0:
            for key in ['cl','cmz','drag']:
                print 'funcs[%s]:'%key
                reg_write(funcs['mdo_tutorial_%s'%key],1e-10,1e-10)
            # Now write the derivatives in the same order the CS will do them:
            print ('Twist[0] Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['twist'][0][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['twist'][0][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['twist'][0][0], 1e-10,1e-10)

            print ('Span Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['span'][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['span'][0], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['span'][0], 1e-10,1e-10)

            print ('shape[13] Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['shape'][0][13], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['shape'][0][13], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['shape'][0][13], 1e-10,1e-10)

            print ('mach Derivatives:')
            reg_write(funcsSens['mdo_tutorial_cl']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_cmz']['mach_mdo_tutorial'], 1e-10,1e-10)
            reg_write(funcsSens['mdo_tutorial_drag']['mach_mdo_tutorial'], 1e-10,1e-10)

    else:
        # For the complex....we just do successive perturbation
        for ii in range(4):
            xRef = {'twist':[0.0, 0.0], 'span':[0.0], 'shape':numpy.zeros(72, dtype='D'), 'mach_mdo_tutorial':0.8}
            if ii == 0:
                xRef['twist'][0] += h*1j
            elif ii == 1:      
                xRef['span'][0] += h*1j
            elif ii == 2:
                xRef['shape'][13] += h*1j
            else:
                xRef['mach_mdo_tutorial']+=h*1j

            ap.setDesignVars(xRef)
            CFDSolver.resetFlow(ap)
            DVGeo.setDesignVars(xRef)
            CFDSolver(ap, writeSolution=False)
            funcs = {}
            CFDSolver.evalFunctions(ap, funcs)

            if MPI.COMM_WORLD.rank == 0:
                if ii == 0:
                    for key in ['cl','cmz','drag']:
                        print 'funcs[%s]:'%key
                        reg_write(numpy.real(funcs['mdo_tutorial_%s'%key]),1e-10,1e-10)

                if ii == 0:
                    print ('Twist[0] Derivatives:')
                elif ii == 1:
                    print ('Span Derivatives:')
                elif ii == 2:
                    print ('shape[13] Derivatives:')
                elif ii == 3:
                    print ('mach Derivatives:')

                for key in ['cl','cmz','drag']:
                    deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                    reg_write(deriv,1e-10,1e-10)

    del CFDSolver
    del mesh

if __name__ == '__main__':
    if len(sys.argv) == 1 or (len(sys.argv) == 2 and 'complex' in sys.argv):
        test1()
        test2()
        test3()
        test4()
        test5()
        test6()
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

            
