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
if 'complex' in sys.argv:
    from adflow import ADFLOW_C as ADFLOW
else:
    from adflow import ADFLOW

# ###################################################################

printHeader('MDO tutorial Viscous Aerodynamic Variables')
# ****************************************************************************
aeroOptions = copy.deepcopy(adflowDefOpts)
h = 1e-40

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
     'blocksplitting': True
 }
)


# Setup aeroproblem, cfdsolver
ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.50, R=287.87,
                 reynolds=50000.0, reynoldsLength=3.25, T=293.15,
                 areaRef=45.5, chordRef=3.25, evalFuncs=['cd','cmz','lift'])
ap.addDV('alpha')
ap.addDV('mach')

CFDSolver = ADFLOW(options=aeroOptions)

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
            print('funcs[%s]:'%key)
            reg_write(funcs['mdo_tutorial_%s'%key],1e-10,1e-10)
        # Now write the derivatives in the same order the CS will do them:
        print('Alpha Derivatives:')
        reg_write(funcsSens['mdo_tutorial_cd']['alpha_mdo_tutorial'], 1e-10,1e-10)
        reg_write(funcsSens['mdo_tutorial_cmz']['alpha_mdo_tutorial'], 1e-10,1e-10)
        reg_write(funcsSens['mdo_tutorial_lift']['alpha_mdo_tutorial'], 1e-10,1e-10)

        print('Mach Derivatives:')
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
                    print('funcs[%s]:'%key)
                    reg_write(numpy.real(funcs['mdo_tutorial_%s'%key]),1e-10,1e-10)

            if ii == 0:
                print('Alpha Derivatives:')
            elif ii == 1:
                print('Mach Derivatives:')

            for key in ['cd','cmz','lift']:
                deriv = numpy.imag(funcs['mdo_tutorial_%s'%key])/h
                reg_write(deriv,1e-10,1e-10)
