############################################################
# DO NOT USE THIS SCRIPT AS A REFERENCE FOR HOW TO USE SUMB
# THIS SCRIPT USES PRIVATE INTERNAL FUNCTIONALITY THAT IS
# SUBJECT TO CHANGE!!
############################################################
import sys, os, copy
from mpi4py import MPI
from baseclasses import AeroProblem
from pyspline import Curve
from pygeo import DVGeometry
from pywarpustruct import USMesh
from mdo_regression_helper import *
from commonUtils import *

# ###################################################################
# DO NOT USE THIS IMPORT STRATEGY! THIS IS ONLY USED FOR REGRESSION
# SCRIPTS ONLY. Use 'from sumb import SUMB' for regular scripts.
sys.path.append(os.path.abspath('../../'))
from python.pySUmb import SUMB
# ###################################################################

# ****************************************************************************
printHeader('Test14: MDO tutorial -- Rans -- Adjoint Test')

# ****************************************************************************
gridFile = '../inputFiles/mdo_tutorial_rans_scalar_jst.cgns'

options = copy.copy(sumbDefOpts)

options.update(
    {'gridfile': gridFile,
     'restartFile':gridFile,
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
     'monitorvariables':['cpu', 'resrho','resturb','cl','cd','cmz','yplus','totalr'],
     'usenksolver':True,
     'l2convergence':1e-14,
     'l2convergencecoarse':1e-4,
     'nkswitchtol':1e-3,
     'adjointl2convergence': 1e-14,
     'frozenturbulence':False,
     'blockSplitting':False,
 }
)
ap = AeroProblem(name='mdo_tutorial', alpha=1.8, mach=0.80, P=20000.0, T=220.0,
                 areaRef=45.5, chordRef=3.25, beta=0.0, 
                 xRef=0.0, yRef=0.0, zRef=0.0, evalFuncs=['fx', 'mz'])

ap.addDV('alpha')
ap.addDV('mach')
CFDSolver = SUMB(options=options)
DVGeo = DVGeometry('../inputFiles/mdo_tutorial_ffd.fmt')
nTwist = 6
DVGeo.addRefAxis('wing', Curve(x=numpy.linspace(5.0/4.0, 1.5/4.0+7.5, nTwist), 
                               y=numpy.zeros(nTwist),
                               z=numpy.linspace(0,14, nTwist), k=2))
def twist(val, geo):
    for i in xrange(nTwist):
        geo.rot_z['wing'].coef[i] = val[i]

DVGeo.addGeoDVGlobal('twist', [0]*nTwist, twist, lower=-10, upper=10, scale=1.0)
DVGeo.addGeoDVLocal('shape', lower=-0.5, upper=0.5, axis='y', scale=10.0)
mesh = USMesh(options={'gridFile':gridFile, 
                   })

CFDSolver.setMesh(mesh)
CFDSolver.setDVGeo(DVGeo)
adjointTest(CFDSolver, ap)

