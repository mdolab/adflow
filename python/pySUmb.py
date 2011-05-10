#!/usr/local/bin/python
'''
pySUmb - A Python interface to SUmb.

Copyright (c) 2008 by Mr.C.A (Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 03/07/2008 11:00$


Developers:
-----------
- Mr. C.A.(Sandy) Mader (SM)
- Dr. Ruben E. Perez (RP)

History
-------
v. 1.0  - Original pyAero Framework Implementation (RP,SM 2008)
'''

__version__ = '$Revision: $'

'''
To Do:
- 
'''


# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import pdb
import time
import copy

# =============================================================================
# External Python modules
# =============================================================================
import numpy
from numpy import real,pi,sqrt

# =============================================================================
# Extension modules
# =============================================================================
from mdo_import_helper import *
exec(import_modules('pyAero_solver'))

# =============================================================================
# SUMB Class
# =============================================================================
class SUMB(AeroSolver):
    
    '''
    SUmb Aerodynamic Analysis Class - Inherited from Solver Abstract Class
    '''
    
    def __init__(self, *args, **kwargs):

        '''
        SUMB Class Initialization
        
        Documentation last updated:  July. 03, 2008 - C.A.(Sandy) Mader
        '''
        
        name = 'SUMB'
        category = 'Three Dimensional CFD'
        def_opts = {
            # Common Paramters
            'gridFile':[str,'default.cgns'],
            'restartFile':[str,'default_restart.cgns'],
            'probName':[str,''],
            'outputDir':[str,'./'],
            'solRestart':[bool,False],
            'writeSolution':[bool,True],

            # Physics Paramters
            'Discretization':[str,'Central plus scalar dissipation'],
            'Smoother':[str,'Runge-Kutta'],
            'equationType': [str,'Euler'],
            'equationMode': [str,'Steady'],
            'flowType':[str,'External'],
            'turblenceModel':[str,'SA'], 
            'useWallFunctions':[bool,False],
            'reynoldsNumber':[float,1e6], 
            'reynoldsLength':[float,1.0], 
            'wallTreatment':[str,'Linear Pressure Extrapolation'],
            'dissipationScalingExponent':[float,0.67],
            'vis4':[float,0.0156],
            'vis2':[float,0.5],
            'vis2Coarse':[float,0.5], 
            'restrictionRelaxation':[float,1.0],

            # Common Paramters
            'nCycles':[int,500],
            'nCyclesCoarse':[int,500],
            'CFL':[float,1.7],
            'CFLCoarse':[float,1.0],
            'Mach':[float,0.5],
            'machCoef':[float,0.5],
            'machGrid':[float,0.0],
            'MGCycle':[str,'3w'],
            'MGStartLevel':[int,-1],

            # Unsteady Paramters           
            'timeIntegrationScheme':[str,'BDF'],
            'timeAccuracy':[int,2],
            'nTimeStepsCoarse':[int,48],
            'nTimeStepsFine':[int,400],
            'deltaT':[float,.010],            

            # Time Spectral Paramters
            'timeIntervals': [int,1],
            'alphaMode':[bool,False],
            'betaMode':[bool,False],
            'machMode':[bool,False],
            'pMode':[bool,False],
            'qMode':[bool,False],
            'rMode':[bool,False],
            'altitudeMode':[bool,False],
            'windAxis':[bool,False],
            'familyRot':[str,''],
            'rotCenter':[list,[0.0,0.0,0.0]],
            'rotRate':[list,[0.0,0.0,0.0]],
            'TSStability': [bool,False],

            # Convergence Paramters
            'L2Convergence':[float,1e-6],
            'L2ConvergenceRel':[float,1e-16],
            'L2ConvergenceCoarse':[float,1e-2], 
            'maxL2DeviationFactor':[float,1.0],

            # Newton-Krylov Paramters
            'useNKSolver':[bool,False],
            'NKLinearSolver':[str,'gmres'],
            'NKSwitchTol':[float,1e-2],
            'NKSubspaceSize':[int,60],
            'NKLinearSolveTol':[float,1e-1],
            'NKPC':[str,'Additive Schwartz'],
            'NKASMOverlap':[int,3],
            'NKPCILUFill':[int,3],
            'NKLocalPCOrdering':[str,'RCM'],
            'NKJacobianLag':[int,10],
            'RKReset':[bool,False],
            'nRKReset':[int,5],

            # Load Balance Paramters
            'blockSplitting':[bool,False],
            'loadImbalance':[float,0.1],
       
            # Misc Paramters
            'metricConversion':[float,1.0],
            'autoSolveRetry':[bool,False],
            'autoAdjointRetry':[bool,False],
            'storeHistory':[bool,False],
            'numberSolutions':[bool,False],
            'printIterations':[bool,False],
            'printTiming':[bool,True],
            'setMonitor':[bool,True],        
            'monitorVariables':[list,['resrho','cl','cd']],
            'surfaceVariables':[list,['cp','vx','vy','vz','mach']],
            'volumeVariables':[list,['resrho']],

            # Reference Paramters
            'refTemp':[float,398.0],
            'refPressure':[float,17654.0],
            'refDensity':[float,.28837],
            'areaAxis':[list,[0.0,1.0,0.0]],

            # Adjoint Paramters
            'adjointL2Convergence':[float,1e-10],
            'adjointL2ConvergenceRel':[float,1e-16],
            'adjointL2ConvergenceAbs':[float,1e-16],
            'adjointDivTol':[float,1e5],
            'approxPC': [bool,False],
            'restartAdjoint':[bool,False],
            'adjointSolver': [str,'GMRES'],
            'adjointMaxIter': [int,500],
            'adjointSubspaceSize' : [int,80],
            'adjointMonitorStep': [int,10],
            'dissipationLumpingParameter':[float,6.0],
            'preconditionerSide': [str,'LEFT'],
            'matrixOrdering': [str,'Nested Dissection'],
            'globalPreconditioner': [str,'Additive Schwartz'],
            'localPreconditioner' : [str,'ILU'],
            'ILUFill': [int,2],
            'ASMOverlap' : [int,5],
            'subKSPSubspaceSize':[int,10],
            'finiteDifferencePC':[bool,True],
            }

        informs = {
            }

        # Load the compiled module using MExt, which allow multiple imports
        if 'sumb' in kwargs:
            self.sumb = kwargs['sumb']
        else:
            sumb_mod = MExt('sumb_parallel')
            self.sumb = sumb_mod._module
        # end if
        
        # Next set the MPI Communicators
        if 'comm' in kwargs:
            self.sumb.communication.sumb_comm_world = kwargs['comm'].py2f()
            self.sumb.communication.sumb_petsc_comm_world = \
                kwargs['comm'].py2f()
            self.sumb.communication.sumb_comm_self  = MPI.COMM_SELF.py2f()
            self.sumb_comm_world = kwargs['comm']
            self.comm = kwargs['comm']
        else: # Set communicators to comm_woprld by default
            self.sumb.communication.sumb_comm_world = MPI.COMM_WORLD.py2f()
            self.sumb.communication.sumb_comm_self  = MPI.COMM_SELF.py2f()
            self.sumb_comm_world = MPI.COMM_WORLD
            self.comm = MPI.COMM_WORLD
        # end if
	
        if 'init_petsc' in kwargs:
            if kwargs['init_petsc']:
                self.sumb.initializepetsc()
        # end if

        # Determine the rank sumb_comm_world size
        self.myid = self.sumb.communication.myid = self.sumb_comm_world.rank
        self.nproc = self.sumb.communication.nproc = self.sumb_comm_world.size

        try:
            self.sumb.communication.sendrequests = numpy.zeros(
                (self.sumb_comm_world.size))
            self.sumb.communication.recvrequests = numpy.zeros(
                (self.sumb_comm_world.size))
        except:
            print "Memory allocation failure for SENDREQUESTS " \
                "and RECVREQUESTS."
            sys.exit(1)
        # end try

        self.callCounter = -1
               
        self.sumb.iteration.standalonemode = False
        self.sumb.iteration.deforming_grid = False

        # Set the frompython flag to true
        self.sumb.killsignals.frompython = True

        # This is SUmb's internal mapping for cost functions
        self.SUmbCostfunctions = \
            {'Lift':self.sumb.costfunctions.costfunclift,
             'Drag':self.sumb.costfunctions.costfuncdrag,
             'Cl'  :self.sumb.costfunctions.costfuncliftcoef,
             'Cd'  :self.sumb.costfunctions.costfuncdragcoef,
             'Fx'  :self.sumb.costfunctions.costfuncforcex,
             'Fy'  :self.sumb.costfunctions.costfuncforcey,
             'Fz'  :self.sumb.costfunctions.costfuncforcez,
             'cFx' :self.sumb.costfunctions.costfuncforcexcoef,
             'cFy' :self.sumb.costfunctions.costfuncforceycoef,
             'cFz' :self.sumb.costfunctions.costfuncforcezcoef,
             'Mx'  :self.sumb.costfunctions.costfuncmomx,
             'My'  :self.sumb.costfunctions.costfuncmomy,
             'Mz'  :self.sumb.costfunctions.costfuncmomz,
             'cMx':self.sumb.costfunctions.costfuncmomxcoef,
             'cMy':self.sumb.costfunctions.costfuncmomycoef,
             'cMz':self.sumb.costfunctions.costfuncmomzcoef,
             'cM0':self.sumb.costfunctions.costfunccm0,
             'cMzAlpha':self.sumb.costfunctions.costfunccmzalpha,
             'cMzAlphaDot':self.sumb.costfunctions.costfunccmzalphadot,
             'cl0':self.sumb.costfunctions.costfunccl0,
             'clAlpha':self.sumb.costfunctions.costfuncclalpha,
             'clAlphaDot':self.sumb.costfunctions.costfuncclalphadot,
             'cfy0':self.sumb.costfunctions.costfunccfy0,
             'cfyAlpha':self.sumb.costfunctions.costfunccfyalpha,
             'cfyAlphaDot':self.sumb.costfunctions.costfunccfyalphadot,
             'cd0':self.sumb.costfunctions.costfunccd0,
             'cdAlpha':self.sumb.costfunctions.costfunccdalpha,
             'cdAlphaDot':self.sumb.costfunctions.costfunccdalphadot,
             'cmzq':self.sumb.costfunctions.costfunccmzq,
             'cmzqDot':self.sumb.costfunctions.costfunccmzqdot,
             'clq':self.sumb.costfunctions.costfuncclq,
             'clqDot':self.sumb.costfunctions.costfuncclqdot,
             }
        
        self.possibleObjectives = \
            { 'lift':'Lift',
              'drag':'Drag',
              'cl':'Cl','cd':'Cd',
              'fx':'Fx','fy':'Fy','fz':'Fz',
              'cfx':'cFx','cfy':'cFy','cfz':'cFz',
              'mx':'Mx','my':'My','mz':'Mz',
              'cmx':'cMx','cmy':'cMy','cmz':'cMz',
              'cm0':'cM0',
              'cmzalpha':'cMzAlpha',
              'cmzalphadot':'cMzAlphaDot',
              'cl0':'cl0',
              'clalpha':'clAlpha',
              'clalphadot':'clAlphaDot',
              'cd0':'cd0',
              'cdalpha':'cdAlpha',
              'cdalphadot':'cdAlphaDot',
              'cmzq':'cmzq',
              'cmzqdot':'cmzqdot',
              'clq':'clq',
              'clqdot':'clqDot',
              'area':'area',
              'volume':'volume'
              }

        self.possibleAeroDVs = {
            'aofa':'adjointvars.ndesignaoa',
            'ssa':'adjointvars.ndesignssa',
            'mach':'adjointvars.ndesignmach',
            'machgrid':'adjointvars.ndesignmachgrid',
            'rotx':'adjointvars.ndesignrotx',
            'roty':'adjointvars.ndesignroty',
            'rotz':'adjointvars.ndesignrotz',
            'rotcenx':'adjointvars.ndesignrotcenx',
            'rotceny':'adjointvars.ndesignrotceny',
            'rotcenz':'adjointvars.ndesignrotcenz',
            'pointrefx':'adjointvars.ndesignpointrefx',
            'pointrefy':'adjointvars.ndesignpointrefy',
            'pointrefz':'adjointvars.ndesignpointrefzy',
            'lengthref':'adjointvars.ndesignlengthref',
            'surfaceref':'adjointvars.ndesignsurfaceref'
            }

        self.aeroDVs = []

        if 'optionMap' in kwargs:
            self.optionMap = kwargs['optionMap']
        else:
            self.optionMap = {\
                # Common Paramters
                'gridFile':{'location':'inputio.gridfile',
                            'len':self.sumb.constants.maxstringlen},
                'restartFile':{'location':'inputio.restartfile',
                               'len':self.sumb.constants.maxstringlen},
                'solRestart':{'location':'inputio.restart'},
                # Physics Paramters
                'Discretization':{'Central plus scalar dissipation':
                                      self.sumb.inputdiscretization.dissscalar,
                                  'Central plus matrix dissipation':
                                      self.sumb.inputdiscretization.dissmatrix,
                                  'Central plus CUSP dissipation':
                                      self.sumb.inputdiscretization.disscusp,
                                  'Upwind':
                                      self.sumb.inputdiscretization.upwind,
                                  'location':
                                      'inputdiscretization.spacediscr'},
                'Smoother':{'Runge-Kutta':
                                self.sumb.inputiteration.rungekutta,
                            'LU-SGS':
                                self.sumb.inputiteration.nllusgs,
                            'LU-SGS Line':
                                self.sumb.inputiteration.nllusgsline,
                            'location':
                                'inputiteration.smoother'},

                'equationType':{'Euler':
                                    self.sumb.inputphysics.eulerequations,
                                'Laminar NS':
                                    self.sumb.inputphysics.nsequations,
                                'RANS':
                                     self.sumb.inputphysics.ransequations,
                                'location':
                                    'inputphysics.equations'},
                'equationMode':{'Steady':
                                   self.sumb.inputphysics.steady,
                               'Unsteady':
                                   self.sumb.inputphysics.unsteady,
                               'Time Spectral':
                                   self.sumb.inputphysics.timespectral,
                               'location':
                                   'inputphysics.equationmode'},
                'flowType':{'Internal':
                                self.sumb.inputphysics.internalflow,
                            'External':
                                self.sumb.inputphysics.externalflow,
                            'location':
                                'inputphysics.flowtype'},
                'turblenceModel':{'Baldwin Lomax':
                                      self.sumb.inputphysics.baldwinlomax,
                                  'SA':
                                      self.sumb.inputphysics.spalartallmaras,
                                  'SAE':
                                      self.sumb.inputphysics.spalartallmarasedwards,
                                  'K Omega Wilcox':
                                      self.sumb.inputphysics.komegawilcox,
                                  'K Omega Modified':
                                      self.sumb.inputphysics.komegamodified,
                                  'Ktau':
                                      self.sumb.inputphysics.ktau,
                                  'Menter SST':
                                      self.sumb.inputphysics.mentersst,
                                  'v2f':
                                      self.sumb.inputphysics.v2f,
                                  'location':
                                      'inputphysics.turbmodel'},
                'useWallFunctions':{'location':'inputphysics.wallfunctions'},
                'reynoldsNumber':{'location':'inputphysics.reynolds'},
                'reynoldsLength':{'location':'inputphysics.reynoldslength'},
                'wallTreatment':{'Linear Pressure Extrapolation':
                                     self.sumb.inputdiscretization.linextrapolpressure,
                                 'Constant Pressure Extrapolation':
                                     self.sumb.inputdiscretization.constantpressure,
                                 'Quadratic Pressure Extrapolation':
                                     self.sumb.inputdiscretization.quadextrapolpressure,
                                 'Normal Momentum':
                                     self.sumb.inputdiscretization.normalmomentum,
                                 'location':
                                     'inputdiscretization.wallbctreatment'},
          


                'dissipationScalingExponent':{'location':'inputdiscretization.adis'},
                'vis4':{'location':'inputdiscretization.vis4'},
                'vis2':{'location':'inputdiscretization.vis2'},
                'vis2Coarse':{'location':'inputdiscretization.vis2coarse'},
                'restrictionRelaxation':{'location':'inputiteration.fcoll'},
                
                # Common Paramters
                'nCycles':{'location':'inputiteration.ncycles'},
                'nCyclesCoarse':{'location':'inputiteration.ncyclescoarse'},
                'CFL':{'location':'inputiteration.cfl'},        
                'CFLCoarse':{'location':'inputiteration.cflcoarse'},        
                'Mach':{'location':'inputphysics.mach'},
                'machCoef':{'location':'inputphysics.machcoef'},
                'machGrid':{'location':'inputphysics.machgrid'},
                'MGCycle':{'location':'localmg.mgdescription',
                           'len':self.sumb.constants.maxstringlen},
                'MGStartLevel':{'location':'inputiteration.mgstartlevel'},

                # Unsteady Params
                'timeIntegrationScheme':{'BDF':
                                             self.sumb.inputunsteady.bdf,
                                         'explicitRK':
                                             self.sumb.inputunsteady.explicitrk,
                                         'inplicitRK':
                                             self.sumb.inputunsteady.implicitrk,
                                         'MD':
                                             self.sumb.inputunsteady.md,
                                         'location':
                                             'inputunsteady.timeintegrationscheme'},
                'timeAccuracy':{'location':'inputunsteady.timeaccuracy'},
                'nTimeStepsCoarse':{'location':'inputunsteady.ntimestepscoarse'},
                'nTimeStepsFine':{'location':'inputunsteady.ntimestepsfine'},
                'deltaT':{'location':'inputunsteady.deltat'},

                # Time Spectral Paramters
                'timeIntervals':{'location':'inputtimespectral.ntimeintervalsspectral'},
                'alphaMode':{'location':'inputtsstabderiv.tsalphamode'},
                'betaMode':{'location':'inputtsstabderiv.tsbetamode'},
                'machMode':{'location':'inputtsstabderiv.tsmachmode'},
                'pMode':{'location':'inputtsstabderiv.tspmode'},
                'qMode':{'location':'inputtsstabderiv.tsqmode'},
                'rMode':{'location':'inputtsstabderiv.tsrmode'},
                'altitudeMode':{'location':'inputtsstabderiv.tsaltitudemode'},
                'windAxis':{'location':'inputtsstabderiv.usewindaxis'},
                'rotCenter':{'location':'inputmotion.rotpoint'},
                'rotRate':{'location':'inputmotion.rotrate'},
                'TSStability':{'location':'inputtsstabderiv.tsstability'},

                # Convergence Paramters
                'L2Convergence':{'location':'inputiteration.l2conv'},
                'L2ConvergenceRel':{'location':'inputiteration.l2convrel'},
                'L2ConvergenceCoarse':{'location':'inputiteration.l2convcoarse'},
                'maxL2DeviationFactor':{'location':'inputiteration.maxl2deviationfactor'},

                # Newton-Krylov Paramters
                'useNKSolver':{'location':'nksolvervars.usenksolver'},
                'NKLinearSolver':{'gmres':'gmres',
                                  'fgmres':'fgmres',
                                  'location':
                                      'nksolvervars.ksp_solver_type',
                                  'len':self.sumb.constants.maxstringlen},

                'NKSwitchTol':{'location':'nksolvervars.nk_switch_tol'},
                'NKSubspaceSize':{'location':'nksolvervars.ksp_subspace'},
                'NKLinearSolveTol':{'location':'nksolvervars.ksp_rtol'},
                'NKPC':{'BlockJacobi':'bjacobi',
                         'Jacobi':'jacobi',
                        'Additive Schwartz':'asm',
                        'location':
                            'nksolvervars.global_pc_type',
                        'len':self.sumb.constants.maxstringlen},
                'NKASMOverlap':{'location':'nksolvervars.asm_overlap'},
                'NKPCILUFill':{'location':'nksolvervars.local_pc_ilu_level'},               
                'NKLocalPCOrdering':{'Natural':'natural',
                                     'RCM':'rcm',
                                     'Nested Dissection':'nd',
                                     'One Way Dissection':'1wd',
                                     'Quotient Minimum Degree':'qmd',
                                     'location':
                                         'nksolvervars.local_pc_ordering',
                                     'len':self.sumb.constants.maxstringlen},
                'NKMaxLinearKspIts':{'location':'nksolvervars.ksp_max_it'},
                'NKJacobianLag':{'location':'nksolvervars.jacobian_lag'},
                'RKReset':{'location':'nksolvervars.rkreset'},
                'nRKReset':{'location':'nksolvervars.nrkreset'},
                'NKFiniteDifferencePC':{'location':'nksolvervars.nkfinitedifferencepc'},
                
                # Load Balance Paramters
                'blockSplitting':{'location':'inputparallel.splitblocks'},
                'loadImbalance':{'location':'inputparallel.loadimbalance'},

                # Misc Paramters

                'printIterations':{'location':'inputiteration.printiterations'},
                'printTiming':{'location':'inputadjoint.printtiming'},
                'setMonitor':{'location':'inputadjoint.setmonitor'},

                # Reference Params
                'refTemp':{'location':'flowvarrefstate.tref'},
                'refPressure':{'location':'flowvarrefstate.pref'},
                'refDensity':{'location':'flowvarrefstate.rhoref'},
               
                # Adjoint Params
                'adjointL2Convergence':{'location':'inputadjoint.adjreltol'},
                'adjointL2ConvergenceRel':{'location':'inputadjoint.adjreltolrel'},
                'adjointL2ConvergenceAbs':{'location':'inputadjoint.adjabstol'},
                'adjointDivTol':{'location':'inputadjoint.adjdivtol'},
                'approxPC':{'location':'inputadjoint.approxpc'},
                'restartAdjoint':{'location':'inputadjoint.restartadjoint'},
                'adjointSolver':{'GMRES':
                                     self.sumb.inputadjoint.petscgmres,
                                 'FGMRES':
                                     self.sumb.inputadjoint.petscfgmres,
                                 'BiCGStab':
                                     self.sumb.inputadjoint.petscbicgstab,
                                 'CG':
                                     self.sumb.inputadjoint.petsccg,
                                 'location':
                                     'inputadjoint.adjointsolvertype'},
                'adjointMaxIter':{'location':'inputadjoint.adjmaxiter'},
                'adjointSubspaceSize':{'location':'inputadjoint.adjrestart'},
                'adjointMonitorStep':{'location':'inputadjoint.adjmonstep'},
                'dissipationLumpingParameter':{'location':'inputdiscretization.sigma'},
                'preconditionerSide':{'LEFT':
                                           self.sumb.inputadjoint.left,
                                       'RIGHT':
                                           self.sumb.inputadjoint.right,
                                       'location':
                                          'inputadjoint.pcside'},
                'matrixOrdering':{'Natural':
                                      self.sumb.inputadjoint.natural,
                                  'RCM':
                                      self.sumb.inputadjoint.reversecuthillmckee,
                                  'Nested Dissection':
                                      self.sumb.inputadjoint.nesteddissection,
                                  'One Way Dissection':
                                       self.sumb.inputadjoint.onewaydissection,
                                  'Quotient Minimum Degree':
                                      self.sumb.inputadjoint.quotientminimumdegree,
                                  'location':
                                      'inputadjoint.matrixordering'},
                'globalPreconditioner':{'BlockJacobi':
                                            self.sumb.inputadjoint.blockjacobi,
                                        'Jacobi':
                                            self.sumb.inputadjoint.jacobi,
                                        'Additive Schwartz':
                                            self.sumb.inputadjoint.additiveschwartz,
                                        'location':
                                            'inputadjoint.precondtype'},
                'localPreconditioner':{'ILU':
                                           self.sumb.inputadjoint.ilu,
                                       'ICC':
                                           self.sumb.inputadjoint.icc,
                                       'LU':
                                           self.sumb.inputadjoint.lu,
                                       'Cholesky':
                                           self.sumb.inputadjoint.cholesky,
                                       'location':
                                           'inputadjoint.localpctype'},
                'ILUFill':{'location':'inputadjoint.filllevel'},
                'ASMOverlap':{'location':'inputadjoint.overlap'},
                'finiteDifferencePC':{'location':'inputadjoint.finitedifferencepc'},
                'subKSPSubspaceSize':{'location':'inputadjoint.subkspsubspacesize'},
                
                }                
        # end if

        # These "ignore_options" are NOT actually, ignore, rather,
        # they DO NOT GET SET IN THE FORTRAN CODE. Rather, they are
        # used strictly in Python
        if 'ignore_options' in kwargs:
            self.ignore_options = kwargs['ignore_options']
        else:
            self.ignore_options = [
                'defaults',
                'storeHistory',
                'numberSolutions',
                'writeSolution',
                'familyRot',  # -> Not sure how to do
                'areaAxis',
                'autoSolveRetry',
                'autoAdjointRetry'
                ]
        # end if
        
        if 'special_options' in kwargs:
            self.special_options = kwargs['special_options']
        else:
            self.special_options = ['surfaceVariables',
                                    'volumeVariables',
                                    'monitorVariables',
                                    'metricConversion',
                                    'outputDir',
                                    'probName']
        # end if

        self.storedADjoints = {}

        # Set default values --- actual options will be set when
        # aero_solver is initialized
        self.sumb.setdefaultvalues()

        # Initialize the inherited aerosolver
        AeroSolver.__init__(\
            self, name, category, def_opts, informs, *args, **kwargs)
        self.sumb.inputio.autoparameterupdate = False

        # Setup the external Mesh Warping
        self._update_geom_info = False
        if 'mesh' in kwargs and kwargs['mesh']:
            self.mesh = kwargs['mesh']
        else:
            self.mesh = SUmbDummyMesh()
        # end if

        # Set Flags that are used to keep of track of what is "done"
        # in fortran
        self.allInitialized = False    # All flow solver initialization   
        self.adjointPreprocessed = False
        self.adjointMatrixSetup = False # Adjoint matrices assembled
        self.adjointRHS         = None # When this is setup, it has
                                       # the current objective
        
        self._update_geom_info = True
        self._update_period_info = True
        self._update_vel_info = True
        self.solve_failed = False
        self.adjoint_failed = False
        self.dtype = 'd'
        # Write the intro message
        self.sumb.writeintromessage()

        return

    def initialize(self, aero_problem, *args,**kwargs):
        '''
        Run High Level Initialization 
        
        Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
        '''
        
        if self.allInitialized == True:
            return
        
        # Set periodic paramters
        self.setPeriodicParams(aero_problem)
  
        # Make sure all the params are ok
        for option in self.options:
            if option != 'defaults':
                self.setOption(option, self.options[option][1])
            # end if
        # end for
      
        # Do the remainder of the operations that would have been done
        # had we read in a param file
        self.sumb.iteration.deforming_grid = True
        self.sumb.dummyreadparamfile()

        self.setMachNumber(aero_problem)
        self.setPeriodicParams(aero_problem)

        #This is just to flip the -1 to 1 possibly a memory issue?
        self.sumb.inputio.storeconvinneriter = \
            abs(self.sumb.inputio.storeconvinneriter)

        if(self.myid ==0):
            print ' -> Partitioning and Reading Grid'
        
        self.sumb.partitionandreadgrid()
       
        if(self.myid==0):
            print ' -> Preprocessing'
        self.sumb.preprocessing()

        if(self.myid==0):
            print ' -> Initializing flow'
        self.sumb.initflow()


        # Create dictionary of variables we are monitoring
        nmon = self.sumb.monitor.nmon
        self.monnames = {}
        for i in range(nmon):
            self.monnames[string.strip(
                    self.sumb.monitor.monnames[i].tostring())] = i
        # end for
      
        # Setup External Warping
        meshInd = self.getMeshIndices()
        self.mesh.setExternalMeshIndices(meshInd)
       
        forceInd = self.getForceIndices()
        self.mesh.setExternalForceIndices(forceInd)
        
        # Solver is initialize
        self.allInitialized = True
        self.sumb.preprocessingadjoint()
        return

    def setInflowAngle(self,aero_problem):
        '''
        Set the alpha and beta fromthe desiggn variables
        '''
        
        [velDir,liftDir,dragDir] = self.sumb.adjustinflowangleadj(\
            (aero_problem._flows.alpha*(pi/180.0)),
            (aero_problem._flows.beta*(pi/180.0)),
            aero_problem._flows.liftIndex)
        self.sumb.inputphysics.veldirfreestream = velDir
        self.sumb.inputphysics.liftdirection = liftDir
        self.sumb.inputphysics.dragdirection = dragDir

        if self.sumb.inputiteration.printiterations:
            if self.myid == 0:
                print '-> Alpha...',
                print aero_problem._flows.alpha*(pi/180.0),
                print aero_problem._flows.alpha

        #update the flow vars
        self.sumb.updateflow()
        self._update_vel_info = True
        return
    
    def setReferencePoint(self,aero_problem):
        '''
        Set the reference point for rotations and moment calculations
        '''
        
        self.sumb.inputphysics.pointref[0] = aero_problem._refs.xref\
            *self.metricConversion
        self.sumb.inputphysics.pointref[1] = aero_problem._refs.yref\
            *self.metricConversion
        self.sumb.inputphysics.pointref[2] = aero_problem._refs.zref\
            *self.metricConversion
#         self.sumb.inputmotion.rotpoint[0] = aero_problem._refs.xrot\
#             *self.metricConversion
#         self.sumb.inputmotion.rotpoint[1] = aero_problem._refs.yrot\
#             *self.metricConversion
#         self.sumb.inputmotion.rotpoint[2] = aero_problem._refs.zrot\
#             *self.metricConversion
        #update the flow vars
        self.sumb.updatereferencepoint()
        self._update_vel_info = True
        return

    def setRotationRate(self,aero_problem):
        '''
        Set the rotational rate for the grid
        '''
        a  = sqrt(self.sumb.flowvarrefstate.gammainf*\
                      self.sumb.flowvarrefstate.pinfdim/ \
                      self.sumb.flowvarrefstate.rhoinfdim)
        V = (self.sumb.inputphysics.machgrid+self.sumb.inputphysics.mach)*a
        
        p = aero_problem._flows.phat*V/aero_problem._refs.bref
        q = aero_problem._flows.qhat*2*V/aero_problem._refs.cref
        r = aero_problem._flows.rhat*V/aero_problem._refs.bref

        self.sumb.updaterotationrate(p,r,q)
        self._update_vel_info = True
        return
    
    def setRefArea(self,aero_problem):
        self.sumb.inputphysics.surfaceref = aero_problem._refs.sref*self.metricConversion**2
        self.sumb.inputphysics.lengthref = aero_problem._refs.cref*self.metricConversion
        
        return

    def setPeriodicParams(self,aero_problem):
        '''
        Set the frequecy and amplitude of the oscillations
        '''
        if  self.getOption('alphaMode'):
            self.sumb.inputmotion.degreepolalpha = int(aero_problem._flows.degreePol)
            self.sumb.inputmotion.coefpolalpha = aero_problem._flows.coefPol
            self.sumb.inputmotion.omegafouralpha   = aero_problem._flows.omegaFourier
            #if self.myid==0: print 'omega',aero_problem._flows.omegaFourier,self.sumb.inputmotion.gridmotionspecified
            self.sumb.inputmotion.degreefouralpha  = aero_problem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffouralpha = aero_problem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffouralpha = aero_problem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif  self.getOption('betaMode'):
            self.sumb.inputmotion.degreepolmach = int(aero_problem._flows.degreePol)
            self.sumb.inputmotion.coefpolmach = aero_problem._flows.coefPol
            self.sumb.inputmotion.omegafourbeta   = aero_problem._flows.omegaFourier
            self.sumb.inputmotion.degreefourbeta  = aero_problem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourbeta = aero_problem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourbeta = aero_problem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('machMode'):
            self.sumb.inputmotion.degreepolmach = int(aero_problem._flows.degreePol)
            self.sumb.inputmotion.coefpolmach = aero_problem._flows.coefPol
            self.sumb.inputmotion.omegafourmach   = aero_problem._flows.omegaFourier
            self.sumb.inputmotion.degreefourmach  = aero_problem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourmach = aero_problem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourmach = aero_problem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif  self.getOption('pMode'):
            ### add in lift axis dependence
            self.sumb.inputmotion.degreepolxrot = int(aero_problem._flows.degreePol)
            self.sumb.inputmotion.coefpolxrot = aero_problem._flows.coefPol
            self.sumb.inputmotion.omegafourxrot = aero_problem._flows.omegaFourier
            self.sumb.inputmotion.degreefourxrot  = aero_problem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourxrot = aero_problem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourxrot = aero_problem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('qMode'):
            self.sumb.inputmotion.degreepolzrot = int(aero_problem._flows.degreePol)
            self.sumb.inputmotion.coefpolzrot = aero_problem._flows.coefPol
            self.sumb.inputmotion.omegafourzrot = aero_problem._flows.omegaFourier
            self.sumb.inputmotion.degreefourzrot  = aero_problem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourzrot = aero_problem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourzrot = aero_problem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('rMode'):
            self.sumb.inputmotion.degreepolyrot = int(aero_problem._flows.degreePol)
            self.sumb.inputmotion.coefpolyrot = aero_problem._flows.coefPol
            self.sumb.inputmotion.omegafouryrot = aero_problem._flows.omegaFourier
            self.sumb.inputmotion.degreefouryrot  = aero_problem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffouryrot = aero_problem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffouryrot = aero_problem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        #endif

        self._update_period_info = True
        self._update_geom_info = True
        self._update_vel_info = True
 
        return

    def setMachNumber(self,aero_problem):
        '''
        Set the mach number for the problem...
        '''
        if self.getOption('familyRot')!='':
            Rotating = True
        else:
            Rotating = False
        #endif
        #print 'Rotating',Rotating,self.getOption('familyRot'),self.getOption('equationMode').lower()
        if Rotating or self.getOption('equationMode').lower()=='time spectral':
            self.sumb.inputphysics.mach = 0.0
            self.sumb.inputphysics.machcoef = aero_problem._flows.mach
            self.sumb.inputphysics.machgrid = aero_problem._flows.mach
            self.sumb.inputmotion.gridmotionspecified = True
        else:
            self.sumb.inputphysics.mach = aero_problem._flows.mach 
            self.sumb.inputphysics.machcoef = aero_problem._flows.mach
            self.sumb.inputphysics.machgrid = 0.0
        #end

        return

    def resetFlow(self):
        '''
        Reset the flow for the complex derivative calculation
        '''
        mgLvlSave =  self.sumb.inputiteration.mgstartlevel
        self.sumb.inputiteration.mgstartlevel = 1
        self.sumb.setuniformflow()
        self.sumb.inputiteration.mgstartlevel = mgLvlSave
        
        return

    def getDensity(self):
        '''
        Get the density for this flow solution
        '''
        rho = self.sumb.flowvarrefstate.rhoinfdim

        return rho

    def getVelocity(self):
        '''
        Get the velocity for this flow solution
        '''
        a  = sqrt(self.sumb.flowvarrefstate.gammainf*\
                      self.sumb.flowvarrefstate.pinfdim/ \
                      self.sumb.flowvarrefstate.rhoinfdim)
        V = (self.sumb.inputphysics.machgrid+self.sumb.inputphysics.mach)*a

        return V

    def __solve__(self, aero_problem, nIterations=500, *args, **kwargs):
        
        '''
        Run Analyzer (Analyzer Specific Routine)
        
        Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
        '''
        # As soon as we run more iterations, adjoint matrices and
        # objective partials (dIdw,dIdx,dIda) are not valid so set
        # their flag to False
        self.adjointMatrixSetup = False 
        self.adjointRHS         = None
        self.callCounter += 1

        # Run Initialize, if already run it just returns.
        self.initialize(aero_problem,*args,**kwargs)

        #set inflow angle,refpoint etc.
        self.setMachNumber(aero_problem)
        self.setPeriodicParams(aero_problem)
        self.setInflowAngle(aero_problem)
        self.setReferencePoint(aero_problem)
        self.setRotationRate(aero_problem)
        self.setRefArea(aero_problem)

        # Run Solver
        t0 = time.time()

        storeHistory = self.getOption('storeHistory')

        # set the number of cycles for this call
        self.sumb.inputiteration.ncycles = nIterations

        # Cold Start -- First Run -- No Iterations Done
        if (self.sumb.monitor.niterold == 0 and 
            self.sumb.monitor.nitercur == 0 and 
            self.sumb.iteration.itertot == 0):

            if self.myid == 0:
                nn = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles
                self.sumb.deallocconvarrays()
                self.sumb.allocconvarrays(nn)
            # end if
        else:
            # More Time Steps / Iterations OR a restart
            # Reallocate convergence history array and time array
            # with new size, storing old values from previous runs
            if storeHistory:
                self._extendConvArray(
                    self.sumb.inputiteration.ncycles)
                self.sumb.monitor.niterold  = self.sumb.monitor.nitercur+1
            else:
                self.sumb.monitor.nitercur  = 0
                self.sumb.monitor.niterold  = 1
                self._clearConvArray()
            #endif

            self.sumb.iteration.itertot = 0
            self.sumb.inputiteration.mgstartlevel = 1
        # end if

        self.sumb.killsignals.routinefailed = False
        self._updatePeriodInfo()
        self._updateGeometryInfo()
        self._updateVelocityInfo()

        # Check to see if the above update routines failed.
        self.sumb.killsignals.routinefailed = \
            self.comm.allreduce(
            bool(self.sumb.killsignals.routinefailed), op=MPI.LOR)

        if self.sumb.killsignals.routinefailed:
            self.solve_failed = True
            return

        # Call the Solver
        if ('MDCallBack' in kwargs):
            self.sumb.solverunsteadymd(kwargs['MDCallBack'])
        else:
            self.sumb.solver()
        # end if


        # If the solve failed, reset the flow for the next time through
        if self.sumb.killsignals.routinefailed:
            mpiPrint('Resetting flow due to failed flow solve...')
            self.resetFlow() # Always reset flow if it failed
            if self.getOption('autoSolveRetry'): # Try the solver again
                self.sumb.solver()
                if self.sumb.killsignals.routinefailed:
                    self.solve_failed = True
                else:
                    self.solve_failed = False
                # end if
            # end if
        # end if 

        if self.sumb.killsignals.routinefailed:
            self.solve_failed = True
        else:
            self.solve_failed = False
        # end if

        sol_time = time.time() - t0

        if self.getOption('printTiming') and self.myid == 0:
            print 'Solution Time',sol_time
        # end if

        # Post-Processing -- Write Solutions
        if self.getOption('writeSolution'):
            base = self.getOption('outputDir') + '/' + self.getOption('probName')
            volname = base + '_vol.cgns'
            surfname = base + '_surf.cgns'

            if self.getOption('numberSolutions'):
                volname = base + '_vol%d.cgns'%(self.callCounter)
                surfname = base + '_surf%d.cgns'%(self.callCounter)
            #endif
            self.writeVolumeSolutionFile(volname)
            self.writeSurfaceSolutionFile(surfname)
        # end if
            
        if self.getOption('TSStability'):
            self.computeStabilityParameters()
        #endif
        
        return

    def _extendConvArray(self,nExtend):
        '''Generic function to extend the current convInfo array by nExtend'''
        if self.myid == 0:
            # Copy Current Values
            temp = copy.deepcopy(self.sumb.monitor.convarray)

            # Deallocate Array
            self.sumb.deallocconvarrays()

            # Allocate New Size
            self.sumb.allocconvarrays(temp.shape[0]+nExtend)

            # Copy temp data back in
            self.sumb.monitor.convarray[:temp.shape[0],:] = copy.deepcopy(temp)
        # end if

        return

    def _clearConvArray(self):
        '''Generic Function to clear the convInfo array. THIS KEEPS
        THE FIRST value for future reference!!!'''
        if self.myid == 0:
            temp=copy.deepcopy(self.sumb.monitor.convarray[0,:,:])
            self.sumb.deallocconvarrays()
            self.sumb.allocconvarrays(self.sumb.inputiteration.ncycles+1)
            self.sumb.monitor.convarray[0,:,:] = temp
        # end if
            
        return

    def solveCL(self,aeroProblem,CL_star,nIterations=500,alpha0=0,
                delta=0.5,tol=1e-3):
        '''This is a simple secant method search for solving for a
        fixed CL. This really should only be used to determine the
        starting alpha for a lift constraint in an optimization.

        Input:  aeroProblem -> aerodynamic problem definition
                CL_star     -> Target CL
                nIterations -> Number of CFD iterations to run (same as 
                               input to __solve__
                alpha0      -> Initial guess for secant search (deg)
                delta       -> Initial step direction (deg)
                tol         -> Absolute tolerance for CL convergence
        Output: aeroProblem._flows.alpha is updated with correct alpha
        '''

        anm2 = alpha0
        anm1 = alpha0 + delta

        # Solve for the n-2 value:
        aeroProblem._flows.alpha = anm2
        self.__solve__(aeroProblem,nIterations=nIterations)
        sol = self.getSolution()
        fnm2 =  sol['cl'] - CL_star

        for iIter in xrange(20):
            # We need to reset the flow since changing the alpha leads
            # to problems with the NK solver
            self.resetFlow()

            # Set current alpha
            aeroProblem._flows.alpha = anm1

            # Solve for n-1 value (anm1)
            self.__solve__(aeroProblem,nIterations=nIterations)
            sol = self.getSolution()
            fnm1 =  sol['cl'] - CL_star
            
            # Secant Update
            anew = anm1 - fnm1*(anm1-anm2)/(fnm1-fnm2)

            # Shift n-1 values to n-2 values
            fnm2 = fnm1
            anm2 = anm1
     
            # Se the n-1 alpha value from update
            anm1 = anew
            
            # Check for convergence
            if abs(fnm1) < tol:
                break
            # end if
        # end for

        # Set the new alpha in the aeroProblem
        aeroProblem._flows.alpha = anew
      
        return

    def getSurfaceCoordinates(self,group_name):
        ''' 
        See MultiBlockMesh.py for more info
        '''
        return self.mesh.getSurfaceCoordinates(group_name)

    def setSurfaceCoordinates(self,group_name,coordinates,reinitialize=True):
        ''' 
        See MultiBlockMesh.py for more info
        '''
        self._update_geom_info = True
        self.mesh.setSurfaceCoordinates(group_name,coordinates,reinitialize)

        return 

    def writeMeshFile(self,filename=None):
        self.sumb.monitor.writegrid = True
        self.sumb.monitor.writevolume = False
        self.sumb.monitor.writesurface = False

        if (filename):
            self.sumb.inputio.solfile[:] = ''
            self.sumb.inputio.solfile[0:len(filename)] = filename

            self.sumb.inputio.newgridfile[:] = ''
            self.sumb.inputio.newgridfile[0:len(filename)] = filename
        # end if

        self.sumb.writesol()

        return

    def writeVolumeSolutionFile(self,filename=None,writeGrid=True):
        """Write the current state of the volume flow solution to a CGNS file.
                Keyword arguments:
        
        filename -- the name of the file (optional)

        """

        self.sumb.monitor.writegrid = writeGrid
        self.sumb.monitor.writevolume = True
        self.sumb.monitor.writesurface = False

        if (filename):
            self.sumb.inputio.solfile[:] = ''
            self.sumb.inputio.solfile[0:len(filename)] = filename

            self.sumb.inputio.newgridfile[:] = ''
            self.sumb.inputio.newgridfile[0:len(filename)] = filename
        # end if

        self.sumb.writesol()

        return

    def writeSurfaceSolutionFile(self,*filename):
        """Write the current state of the surface flow solution to a CGNS file.
        
        Keyword arguments:
        
        filename -- the name of the file (optional)

        """
        if (filename):
            self.sumb.inputio.surfacesolfile[:] = ''
            self.sumb.inputio.surfacesolfile[0:len(filename[0])] = filename[0]
        self.sumb.monitor.writegrid=False
        self.sumb.monitor.writevolume=False
        self.sumb.monitor.writesurface=True
        self.sumb.writesol()

        return

    def writeForceFile(self,file_name,TSInstance=0):
        '''This function collects all the forces and locations and
        writes them to a file with each line having: X Y Z Fx Fy Fz.
        This can then be used to set a set of structural loads in TACS
        for structural only optimization'''
      
        pts = self.getForcePoints()
        nPts = pts.shape[1]
        if nPts > 0:
            forces =  self.sumb.getforces(pts.T).T
        else:
            forces = numpy.empty(pts.shape,dtype=self.dtype)
        # end if

        # Now take the desired time instance
        pts = pts[TSInstance,:,:]
        forces = forces[TSInstance,:,:]
        
        # Now we need to gather the data:
        pts = self.comm.gather(pts, root=0)
        forces = self.comm.gather(forces, root=0)

        # Write out Data only on root proc:
        if self.myid == 0:
            f = open(file_name,'w')
            for iproc in xrange(len(pts)):
                for ipt in xrange(len(pts[iproc])):
                    f.write('%20.15g %20.15g %20.15g %20.15g %20.15g %20.15g\n'%(
                            pts[iproc][ipt,0],pts[iproc][ipt,1],
                            pts[iproc][ipt,2],forces[iproc][ipt,0],
                            forces[iproc][ipt,1],forces[iproc][ipt,2]))
                # end for
            # end for
            f.close()
        # end if

        return 

    def getForces(self,group_name,cfd_force_pts=None):
        ''' Return the forces on this processor. Use
        cfd_force_pts to compute the forces if given
        
        '''

        if cfd_force_pts is None:
            cfd_force_pts = self.getForcePoints()
        # end if
        if len(cfd_force_pts) > 0:
            forces = self.sumb.getforces(cfd_force_pts.T).T
        else:
            forces = numpy.empty((0),dtype=self.dtype)
        # end if
            
        return self.mesh.solver_to_warp_force(group_name,forces)

    def getForcePoints(self):
        [npts,nTS] = self.sumb.getforcesize()
        if npts > 0:
            return self.sumb.getforcepoints(npts,nTS).T
        else:
            return numpy.empty((nTS,0,3),dtype=self.dtype)
        # end if

    def verifyForces(self,cfd_force_pts=None):

        # Adjoint must be initialized for force verification

        self.initAdjoint()
        if cfd_force_pts is None:
            cfd_force_pts = self.getForcePoints()
        # end if

        self.sumb.verifyforces(cfd_force_pts.T)

        return


    def globalNKPreCon(self,in_vec):
        '''This function is ONLY used as a preconditioner to the
        global Aero-Structural system'''

        out_vec = self.sumb.applypc(in_vec)
        
        return out_vec


    def verifydCdx(self,objective,**kwargs):
        '''
        call the reouttine to compare the partial dIda
        against FD
        '''
        # Check to see if the adjoint Matrix is setup:
        if self.myid==0: print 'setting up matrix'
        if not self.adjointMatrixSetup:
            self.sumb.createpetscvars()
            #self.setupAdjoint(forcePoints)
        # end if
        # Short form of objective--easier code reading
        if self.myid==0: print 'possible objectives',objective
        obj = self.possibleObjectives[objective.lower()]
        if self.myid==0: print 'obj',obj
        costFunc =  self.SUmbCostfunctions[obj]
        if self.myid==0: print 'costfunc',costFunc
        self.sumb.verifydcfdx(1,costFunc)
        
        return

    def verifydCdw(self,objective,**kwargs):
	
	self.sumb.verifydcdwfile(1)

	return
    def verifydIdw(self,objective,**kwargs):
        '''
        run compute obj partials, then print to a file...
        '''
        if self.myid==0: print 'setting up vector'
        if not self.adjointMatrixSetup:
            self.sumb.createpetscvars()
            #self.setupAdjoint(forcePoints)
        # end if
        if self.myid==0:print 'computing partials'
        self.computeObjPartials(objective)
        obj = self.possibleObjectives[objective.lower()]
        filename= self.getOption('outputDir') + '/' +'ADw%s'%(obj)
        if self.myid==0:print 'filename',filename
        costFunc =  self.SUmbCostfunctions[obj]
        level = 1
        if self.myid==0:print 'calling verify',level,costFunc,filename
        self.sumb.verifydidwfile(level,costFunc,filename)
        
        return
    def verifydIdx(self,objective,**kwargs):
        '''
        run compute obj partials, then print to a file...
        '''
        if self.myid==0: print 'setting up vector'
        if not self.adjointMatrixSetup:
            self.sumb.createpetscvars()
            #self.setupAdjoint(forcePoints)
        # end if
        if self.myid==0:print 'computing partials'
        self.computeObjPartials(objective)
        obj = self.possibleObjectives[objective.lower()]
        filename= self.getOption('outputDir') + '/' +'ADx%s'%(obj)
        if self.myid==0:print 'filename',filename
        costFunc =  self.SUmbCostfunctions[obj]
        level = 1
        if self.myid==0:print 'calling verify',level,costFunc,filename
        self.sumb.verifydidxfile(level,costFunc,filename)

        return

    def initAdjoint(self, *args, **kwargs):
        '''
        Initialize the Ajoint problem for this test case
        in SUMB
        '''
        
        # Set the index value of each nDesign to -1 -> Don't use by
        # default
        self.sumb.adjointvars.ndesignaoa = -1
        self.sumb.adjointvars.ndesignssa  = -1
        self.sumb.adjointvars.ndesignmach = -1
        self.sumb.adjointvars.ndesignmachgrid = -1
        self.sumb.adjointvars.ndesignrotx = -1
        self.sumb.adjointvars.ndesignroty = -1
        self.sumb.adjointvars.ndesignrotz = -1
        self.sumb.adjointvars.ndesignrotcenx = -1
        self.sumb.adjointvars.ndesignrotceny = -1
        self.sumb.adjointvars.ndesignrotcenz = -1
        self.sumb.adjointvars.ndesignpointrefx = -1
        self.sumb.adjointvars.ndesignpointrefy = -1
        self.sumb.adjointvars.ndesignpointrefz = -1
        self.sumb.adjointvars.ndesignlengthref = -1
        self.sumb.adjointvars.ndesignsurfaceref = -1
        
        # Set the required paramters for the aero-Only design vars:
        self.nDVAero = len(self.aeroDVs)#for debuggin with check all...

        self.sumb.adjointvars.ndesignextra = self.nDVAero
        
        if self.nDVAero >0:           
            self.sumb.adjointvars.dida = numpy.zeros(self.nDVAero)
            for i in xrange(self.nDVAero):
                exec_str = 'self.sumb.' + self.possibleAeroDVs[self.aeroDVs[i]] + \
                           '= %d'%(i)
                # Leave this zero-based since we only need to use it in petsc
                exec(exec_str)
            # end for
        #endif

        #Set the mesh level and timespectral instance for this
        #computation
        
        self.sumb.iteration.currentlevel=1
        self.sumb.iteration.groundlevel=1
        if not self.adjointPreprocessed:
            self.sumb.preprocessingadjoint()
        # end if

        #self.sumb.initializepetsc() -> I don't think we will need this

        return

    def addAeroDV(self,*dvs):
        '''Take in a list of DVs that the flow solver will use in
        addition to shape-type design variables'''
        for i in xrange(len(dvs)):
            if dvs[i] in self.possibleAeroDVs and not dvs[i] in self.aeroDVs:
                self.aeroDVs.append(dvs[i])
            else:
                print 'Warning: %s was not one of the possible AeroDVs'%(dvs[i])
                print 'Full AeroDV list is:'
                print self.possibleAeroDVs
            # end if
        # end for
        return

    def setupAdjoint(self, forcePoints=None, **kwargs):
        '''
        Setup the adjoint matrix for the current solution
        '''
        
        # Destroy the NKsolver to free memory -- Call this even if the
        # solver is not used...a safeguard check is done in Fortran
        self.sumb.destroynksolver()

        if not self.adjointMatrixSetup:
            self.sumb.createpetscvars()
 
            self.sumb.setupallresidualmatrices()

            if forcePoints is None:
                forcePoints = self.getForcePoints()
            # end if
            
            self.sumb.setupcouplingmatrixstruct(forcePoints.T)
            self.sumb.setuppetscksp()
            self.mesh.setupWarpDeriv()
            self.adjointMatrixSetup = True

        # end if

        return

    def printMatrixInfo(self, dRdwT=True, dRdwPre=True, dRdx=True,
                        dRda=True, dSdw=True, dSdx=True,
                        printLocal=False,printSum=True,printMax=False):
        
        # Call sumb matrixinfo function
        self.sumb.matrixinfo(dRdwT,dRdwPre,dRdx,dRda,dSdw,dSdx,
                             printLocal,printSum,printMax)

        return
    
    def releaseAdjointMemory(self):
        '''
        release the PETSc Memory...
        '''

        self.sumb.destroypetscvars()

        return

    def _on_adjoint(self,objective,forcePoints=None,*args,**kwargs):

        # Make sure adjoint is initialize
        self.initAdjoint()

        # Short form of objective--easier code reading
        obj = self.possibleObjectives[objective.lower()]

        # Check to see if the adjoint Matrix is setup:
        if not self.adjointMatrixSetup:
            self.setupAdjoint(forcePoints)
        # end if
       
        # Check to see if the RHS Partials have been computed
        if not self.adjointRHS == obj:
            self.computeObjPartials(obj,forcePoints)
        # end if
        # Check to see if we need to agument the RHS with a structural
        # adjoint:
        if 'structAdjoint' in kwargs and 'group_name' in kwargs:
            group_name = kwargs['group_name']
            phi = kwargs['structAdjoint']
            solver_phi = self.mesh.warp_to_solver_force(group_name,phi)
            self.sumb.agumentrhs(solver_phi)
        # end if

        # If we have saved adjoints, 
        if self.getOption('restartAdjoint'):
            # Objective is already stored, so just set it
            if obj in self.storedADjoints.keys():
                self.sumb.setadjoint(self.storedADjoints[obj])
            else:
                # Objective is not yet run, allocated zeros and set
                self.storedADjoints[obj]= numpy.zeros(self.getStateSize(),float)
                self.sumb.setadjoint(self.storedADjoints[obj])
            # end if
        # end if

        # Actually Solve the adjoint system
        self.sumb.solveadjointtransposepetsc()

        # Possibly try another solve
        if self.sumb.killsignals.adjointfailed and self.getOption('restartAdjoint'):
            # Only retry if the following conditions are met:
            # 1. restartAdjoint is true -> that is we were starting from a non-zero starting point
            # 2. The stored adjoint must have been already set at least once; that is we've already tried one solve
            
            if  obj in self.storedADjoints.keys():
                self.storedADjoints[obj][:] = 0.0 # Always reset a stored adjoint 

                if self.getOption('autoAdjointRetry'):
                    self.sumb.solveadjointtransposepetsc()
                # end if
            # end if
        # end if

        # Now set the flags and possibly reset adjoint
        if self.sumb.killsignals.adjointfailed == False:
            self.adjoint_failed = False
            # Copy out the adjoint to store
            if self.getOption('restartAdjoint'):
                self.storedADjoints[obj] =  \
                    self.sumb.getadjoint(self.getStateSize())
            # end if
        else:
            self.adjoint_failed = True
            # Reset stored adjoint
            if self.getOption('restartAdjoint'):
                self.storedADjoints[obj][:] = 0.0
            # end if
        # end if
       
        return

    def totalSurfaceDerivative(self,objective):
        # The adjoint vector is now calculated so perform the
        # following operation to produce dI/dX_surf:
        # (p represents partial, d total)
        # dI/dX_s = pI/pX_s - (dXv/dXs)^T * ( dRdX_v^T * psi)
        # 
        # The derivative wrt the surface captures the effect of ALL
        # GLOBAL Multidisciplinary variables -- any DV that changes
        # the surface. 

        obj = self.possibleObjectives[objective.lower()]
        
        if obj in ['area','volume']: # Possibly add more Direct objectives here...
            if obj == 'area':
                self.mesh.warp.computeareasensitivity(self.getOption('areaAxis'))
                dIdXs = self.mesh.getdXs('all')
            # end if

            if obj == 'volume':
                self.mesh.warp.computevolumesensitivity()
                dIdXs = self.mesh.getdXs('all')
            # end if

        else:
            if self.getOption('restartAdjoint'): # Selected stored adjoint
                self.sumb.setadjoint(self.storedADjoints[obj])
            # end if

            # Direct partial derivative contibution 
            dIdxs_1 = self.getdIdx(objective,'all')

            # dIdx contribution for drdx^T * psi
            dIdxs_2 = self.getdRdXvPsi('all')
        
            # Total derivative of the obective with surface coordinates
            dIdXs = dIdxs_1 - dIdxs_2 
        # end if

        return dIdXs

    def totalAeroDerivative(self,objective):
        # The adjoint vector is now calculated. This function as above
        # computes dI/dX_aero = pI/pX_aero - dR/dX_aero^T * psi. The
        # "aero" variables are intrinsic ONLY to the aero
        # discipline. Nothing in the structural process should depend
        # on these functions directly. 

        obj = self.possibleObjectives[objective.lower()]

        if obj in ['area','volume']: # Possibly add more Direct
                                     # objectives here...  These by
                                     # definition have zero dependance
            dIda = numpy.zeros(self.nDVAero)
        else:
            if self.getOption('restartAdjoint'):
                self.sumb.setadjoint(self.storedADjoints[obj])
            # end if

            # Direct partial derivative contibution 
            dIda_1 = self.getdIda(objective)

            # dIda contribution for drda^T * psi
            dIda_2 = self.getdRdaPsi()

            # Total derivative of the obective wrt aero-only DVs
            dIda = dIda_1 - dIda_2
        # end if

        return dIda
        
    def verifyPartials(self):
        ''' Run verifyResiduals to verify that dRdw,dRdx and dRda are
        computed correctly
        '''
        self.sumb.verifypartials()

        return
    
    def verifyResidual(self):
        ''' Run verifyRAdjoint to make sure the node-based-stencil
        routine gives the same results as the original routine
        '''
        self.sumb.verifyradj(1)

        return 

    def computeStabilityParameters(self):
        '''
        run the stability derivative driver to compute the stability parameters
        from the time spectral solution
        '''
        self.sumb.stabilityderivativedriver()
        return

    def _updateGeometryInfo(self):
        """Update the SUmb internal geometry info, if necessary."""
        if (self._update_geom_info):
            self.mesh.warpMesh()
            newGrid = self.mesh.getSolverGrid()
            if newGrid is not None:
                self.sumb.setgrid(self.mesh.getSolverGrid())
                
            self.sumb.updatecoordinatesalllevels()
            self.sumb.updatewalldistancealllevels()
            self.sumb.updateslidingalllevels()
            self.sumb.updatemetricsalllevels()
            self.sumb.updategridvelocitiesalllevels()
            self._update_geom_info = False
        # end if

        return 

    def _updatePeriodInfo(self):
        """Update the SUmb TS period info"""
        if (self._update_period_info):
            self.sumb.updateperiodicinfoalllevels()
            self._update_period_info = False
        # end if

        return 

    def _updateVelocityInfo(self):
        
        if (self._update_vel_info):
            self.sumb.updategridvelocitiesalllevels()
            self._update_vel_info = False
        # end if
        
        return 
            
    def getMonitoringVariables(self):
        """Return a list of the text strings describing the variables being
        monitored.

        """
        return self.monnames.keys()
    

    def getConvergenceHistory(self,name):
        """Return an array of the convergence history for a particular quantity.

        Keyword arguments:

        name -- the text string for a particular quantity

        """
        try:
            index = self.monnames[name]
        except KeyError:
            print "Error: No such quantity '%s'" % name
            return None
        if (self.myid == 0):
            if (self.sumb.monitor.niterold == 0 and
                self.sumb.monitor.nitercur == 0 and
                self.sumb.iteration.itertot == 0):
                history = None
            elif (self.sumb.monitor.nitercur == 0 and
                  self.sumb.iteration.itertot == 0):
                niterold = self.sumb.monitor.niterold[0]	    
                history = self.sumb.monitor.convarray[:niterold+1,index]
            else:
                history = self.sumb.monitor.convarray[:,index]
        else:
            history = None
        history = self.sumb_comm_world.bcast(history)
	
        return history

    def getAdjointResiduals(self):
        '''
        Return the following adjoint residual norms:
        initCFD Norm: Norm the adjoint starts with (zero adjoint)
        startCFD Norm: Norm at the start of adjoint call
        finalCFD Norm: Norm at the end of adjoint call
        '''
        
        startRes = self.sumb.adjointpetsc.adjreshist[0]
        finalIt  = self.sumb.adjointpetsc.adjconvits
        finalRes = self.sumb.adjointpetsc.adjreshist[finalIt]
        fail = self.sumb.killsignals.adjointfailed

        return startRes,finalRes,fail

    def getResNorms(self):
        '''Return the initial, starting and final Res Norms'''
        
        return \
            self.sumb.nksolvervars.totalr0, \
            self.sumb.nksolvervars.totalrstart,\
            self.sumb.nksolvervars.totalrfinal

    def getMeshIndices(self):
        ndof = self.sumb.getnumberlocalnodes()
        indices = self.sumb.getcgnsmeshindices(ndof)

        return indices
    
    def getForceIndices(self):
        ndof = self.sumb.getnumberlocalforcenodes()
        if ndof > 0:
            indices = self.sumb.getcgnsforceindices(ndof)
        else:
            indices = numpy.zeros(0,'intc')
        # end if

        return indices

    def getdRdXvPsi(self,group_name):
        ndof = self.sumb.adjointvars.nnodeslocal*3
        dxv_solver = self.sumb.getdrdxvpsi(ndof)
        self.mesh.WarpDeriv(dxv_solver)
        pforces = self.mesh.getdXs(group_name)
        
        return pforces

    def getdRdaPsi(self):
        if self.nDVAero > 0:
            dIda = self.sumb.getdrdapsi(self.nDVAero)
        else:
            dIda = numpy.zeros((0))
        # end if

        return dIda

    def getdRdwPsi(self):
        dRdwPsi = self.sumb.getdrdwtpsi(self.getStateSize())
        
        return dRdwPsi

    def getdFdxVec(self,group_name,vec):
        # Calculate dFdx * force_pts and return the result
        solver_vec = self.mesh.warp_to_solver_force(group_name,vec)

        dFdxVec = self.sumb.getdfdxvec(solver_vec)

        dFdxVec = self.mesh.solver_to_warp_force(group_name,solver_vec)

        return dFdxVec

    def computeObjPartials(self,objective,forcePoints=None):
        obj = self.possibleObjectives[objective.lower()]

        if forcePoints is None:
            forcePoints = self.getForcePoints()

        obj_num = self.SUmbCostfunctions[obj]
        if not obj == self.adjointRHS:
            self.sumb.computeobjpartials(obj_num,forcePoints.T)
            self.adjointRHS = obj
        # end if
        return 

    def getdIdx(self,objective,group_name,forcePoints=None):

        if forcePoints is None:
            forcePoints = self.getForcePoints()
        # end for
        self.computeObjPartials(objective,forcePoints)

        temp = forcePoints.shape
        
        sizeForcePoints = temp[1]*temp[2]#temp[0]*temp[1]#*temp[2]
        
        #if len(forcePoints > 0):
        if sizeForcePoints > 0:
            dIdpts = self.sumb.getdidx(sizeForcePoints)
        
            dIdpts.reshape(forcePoints[0,:,:].shape)
        
        else:
            dIdpts = numpy.zeros((0),self.dtype)
        # end if

        dIdpts = self.mesh.solver_to_warp_force(group_name,dIdpts)
        return dIdpts

    def getdIda(self,objective,forcePoints=None):
        if forcePoints is None:
            forcePoints = self.getForcePoints()
        # end if
        if self.nDVAero > 0:
            self.computeObjPartials(objective,forcePoints)
        
            dIda_local = self.sumb.adjointvars.dida
    
            # We must MPI all reuduce
            dIda = self.comm.allreduce(dIda_local, op=MPI.SUM)
        else:
            dIda = numpy.zeros((0))
    
        return dIda
        
    def finalizeAdjoint(self):
        '''
        destroy the PESTcKSP context
        '''
        
        self.releaseAdjointMemeory()
        
        return

    def getStateSize(self):
        '''Return the number of degrees of freedom (states) that are
        on this processor'''

        nw     = self.sumb.flowvarrefstate.nw
        ncells = self.sumb.adjointvars.ncellslocal
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return nw*ncells*ntime

    def getStates(self):
        '''Return the states on this processor. Used in aerostructural
        analysis'''
        states = self.sumb.getstates(self.getStateSize())

        return states

    def setStates(self,states):
        ''' Set the states on this processor. Used in aerostructural
        analysis'''

        self.sumb.setstates(states)

        return 

    def getResidual(self):

        '''Return the residual on this processor. Used in aerostructural
        analysis'''
        res    = self.sumb.getres(self.getStateSize())
        
        return res

    def getSolution(self,sps=1):
        '''
        retrieve the solution variables from the solver.
        '''

        # We should return the list of results that is the same as the
        # possibleObjectives list
        self.sumb.getsolution(sps)

        funcVals = self.sumb.costfunctions.functionvalue
        SUmbsolution =  \
            {'lift':funcVals[self.sumb.costfunctions.costfunclift-1],
             'drag':funcVals[self.sumb.costfunctions.costfuncdrag-1],
             'cl'  :funcVals[self.sumb.costfunctions.costfuncliftcoef-1],
             'cd'  :funcVals[self.sumb.costfunctions.costfuncdragcoef-1],
             'fx'  :funcVals[self.sumb.costfunctions.costfuncforcex-1],
             'fy'  :funcVals[self.sumb.costfunctions.costfuncforcey-1],
             'fz'  :funcVals[self.sumb.costfunctions.costfuncforcez-1],
             'cfx' :funcVals[self.sumb.costfunctions.costfuncforcexcoef-1],
             'cfy' :funcVals[self.sumb.costfunctions.costfuncforceycoef-1],
             'cfz' :funcVals[self.sumb.costfunctions.costfuncforcezcoef-1],
             'mx'  :funcVals[self.sumb.costfunctions.costfuncmomx-1],
             'my'  :funcVals[self.sumb.costfunctions.costfuncmomy-1],
             'mz'  :funcVals[self.sumb.costfunctions.costfuncmomz-1],
             'cmx' :funcVals[self.sumb.costfunctions.costfuncmomxcoef-1],
             'cmy' :funcVals[self.sumb.costfunctions.costfuncmomycoef-1],
             'cmz' :funcVals[self.sumb.costfunctions.costfuncmomzcoef-1],
             'cmzalphadot':funcVals[self.sumb.costfunctions.costfunccmzalphadot-1],
             'cmzalpha'   :funcVals[self.sumb.costfunctions.costfunccmzalpha-1],
             'cm0'        :funcVals[self.sumb.costfunctions.costfunccm0-1],
             'clalphadot' :funcVals[self.sumb.costfunctions.costfuncclalphadot-1],
             'clalpha'    :funcVals[self.sumb.costfunctions.costfuncclalpha-1],
             'cl0'        :funcVals[self.sumb.costfunctions.costfunccl0-1],
             'cfyalphadot':funcVals[self.sumb.costfunctions.costfunccfyalphadot-1],
             'cfyalpha'   :funcVals[self.sumb.costfunctions.costfunccfyalpha-1],
             'cfy0'       :funcVals[self.sumb.costfunctions.costfunccfy0-1],
             'cdalphadot' :funcVals[self.sumb.costfunctions.costfunccdalphadot-1],
             'cdalpha'    :funcVals[self.sumb.costfunctions.costfunccdalpha-1],
             'cd0'        :funcVals[self.sumb.costfunctions.costfunccd0-1],
             'cmzqdot'    :funcVals[self.sumb.costfunctions.costfunccmzqdot-1],
             'cmzq'       :funcVals[self.sumb.costfunctions.costfunccmzq-1],
             'clqdot'     :funcVals[self.sumb.costfunctions.costfuncclqdot-1],
             'clq'        :funcVals[self.sumb.costfunctions.costfuncclq-1]
             }
                                                 
        # Also add in 'direct' solutions. Area etc
        A_local = self.mesh.computeArea(self.getOption('areaAxis'))
        SUmbsolution['area'] = self.comm.allreduce(A_local,op=MPI.SUM)

        V_local = self.mesh.computeVolume()
        SUmbsolution['volume'] = self.comm.allreduce(V_local,op=MPI.SUM)
        
        return SUmbsolution
        
    def _on_setOption(self, name, value):
        
        '''
        Set Solver Option Value 
        '''

        # Ignored options do NOT get set in solver

        if name in self.ignore_options:
            return

        # Do special Options individually
        if name in self.special_options:
            if name in ['monitorVariables','surfaceVariables','volumeVariables']:
                varStr = ''
                for i in xrange(len(value)):
                    varStr = varStr + value[i] + '_'
                # end if
                varStr = varStr[0:-1] # Get rid of last '_'
                if name == 'monitorVariables':
                    self.sumb.monitorvariables(varStr)
                if name == 'surfaceVariables':
                    self.sumb.surfacevariables(varStr)
                if name == 'volumeVariables':
                    self.sumb.volumevariables(varStr)
            # end if
            if name == 'metricConversion':
                self.sumb.flowvarrefstate.lref = value
                self.sumb.flowvarrefstate.lrefspecified = True
                self.metricConversion = value
            # end if

            return
        # end if

        # All other options do genericaly by setting value in module:
        # Check if there is an additional mapping to what actually
        # has to be set in the solver

        temp = copy.copy(self.optionMap[name]) # This is the dictionary
        temp.pop('location')
        try:
            temp.pop('len')
        except:
            pass
        # end if

        # If temp has anything left in it, we MUST be able to match to
        # one of them.

        if len(temp) == 0:
            pass
        else:
            value = self.optionMap[name][value]
        # end if
            

        # If value is a string, put quotes around it and make it
        # the correct length, otherwise convert to string
        if isinstance(value,str): 
            spacesToAdd = self.optionMap[name]['len'] - len(value)
            value = '\'' + value + ' '*spacesToAdd + '\''
        else:
            value = str(value)
        # end if
      
        # Exec str is what is actually executed:
        exec_str = 'self.sumb.'+self.optionMap[name]['location'] + '=' + value

        exec(exec_str)

        return
        
    def _on_getOption(self, name):
        
        '''
        Get Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 21, 2008 - Ruben E. Perez
        '''
        
        pass
        
        
    def _on_getInform(self, infocode):
        
        '''
        Get Optimizer Result Information (Optimizer Specific Routine)
        
        Keyword arguments:
        -----------------
        id -> STRING: Option Name
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        # 
        return self.informs[infocode]


class SUmbDummyMesh(object):
    """
    Represents a dummy Multiblock structured Mesh for SUmb
    """
 
    def __init__(self):
        """Initialize the object."""

    def getSurfaceCoordinates(self,group_name):
        ''' 
        Returns a UNIQUE set of ALL surface points belonging to
        group "group_name", that belong to blocks on THIS processor
        '''

        return 

    def setSurfaceCoordinates(self,group_name,coordinates,reinitialize=True):
        ''' 
        Set the UNIQUE set of ALL surface points belonging to
        group "group_name", that belong to blocks on THIS processor
        This must be the same length as the list obtained from
        getSurfaceCoordiantesLocal
        '''

        return 

    def addFamilyGroup(self,group_name,families=None):

        return 
   
# =========================================================================
#                         Interface Functionality
# =========================================================================

    def setExternalMeshIndices(self,ind):
        ''' Take in a set of external indices from another flow solver
        and use this to setup the scatter contexts'''

        return 

    def setExternalForceIndices(self,ind):
        ''' Take in a set of external indices from another flow solver
        and use this to setup the scatter contexts'''

        return 

    def getSolverGrid(self):
        '''
        Return the grid as the external solver want it
        '''
        
        return

    def warp_to_solver_force(self,group_name,disp):

        return 

    def solver_to_warp_force(self,group_name,solver_forces):
    
        return

    def getdXs(self,group_name):

        return 

    def computeArea(self,axis):
        
        return 0.0

    def computeVolume(self):

        return 0.0

# ==========================================================================
#                        Output Functionality
# ==========================================================================

    def writeVolumeGrid(self,file_name):
        '''write volume mesh'''

        return

    def writeSurfaceGrid(self,file_name):
        '''write surface mesh'''

        return

    def writeFEGrid(self,file_name):
        ''' write the FE grid to file '''
        return

# =========================================================================
#                         Utiliy Functions
# =========================================================================

    def getMeshQuality(self,bins=None):
        
        return

# =========================================================================
#                      Mesh Warping Functionality
# =========================================================================

    def warpMesh(self):
        ''' 
        Run either the solid warping scheme or the algebraic
        warping scheme depending on the options
        '''

    def setupWarpDeriv(self):

        return 

    def WarpDeriv(self,solver_dxv):

        return 

    def verifySolidWarpDeriv(self):

        return

#==============================================================================
# SUmb Analysis Test
#==============================================================================
if __name__ == '__main__':
    
    # Test SUmb
    print 'Testing ...'
    sumb = SUMB()
    print sumb

