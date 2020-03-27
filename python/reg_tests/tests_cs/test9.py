 ############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE ADFLOW
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI

from baseclasses import AeroProblem
from pygeo import DVGeometry
import pyspline

from mdo_regression_helper import *
from commonUtils import *

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from adflow import ADFLOW' for regular scripts.
sys.path.append(os.path.abspath('../../'))

if 'complex' in sys.argv:
    from python.pyADflow_C import ADFLOW_C as ADFLOW
    from idwarp import USMesh_C as USMesh
else:
    from python.pyADflow import ADFLOW
    from idwarp import  USMesh
# ###################################################################
# ****************************************************************************
printHeader('MDO tutorial RANS Geometric Variables with IDWarp')
# ****************************************************************************
aeroOptions = copy.deepcopy(adflowDefOpts)

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
     'blocksplitting': True
 }
)
h = 1e-40
# Setup aeroproblem, cfdsolver
ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, R=287.87,
                 altitude=10000.0, areaRef=45.5, chordRef=3.25, evalFuncs=['cl','cmz','drag'])

ap.addDV('alpha')
ap.addDV('mach')
CFDSolver = ADFLOW(options=aeroOptions)
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
mesh = USMesh(options={'gridFile':'../inputFiles/mdo_tutorial_rans.cgns'})
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