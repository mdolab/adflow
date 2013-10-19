#!/usr/bin/python
from __future__ import print_function
'''
pySUmb - A Python interface to SUmb.

Copyright (c) 2008 by Mr.C.A (Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 03/07/2008 11:00$

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)
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
import sys, string, os
import time
import copy

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from baseclasses import AeroSolver
from mdo_import_helper import mpiPrint, MExt, MPI

# =============================================================================
# SUMB Class
# =============================================================================
class SUMB(AeroSolver):
    '''
    SUmb Aerodynamic Analysis Class - Inherited from Solver Abstract Class
    '''
    def __init__(self, comm=None, options=None, mesh=None, **kwargs):
        '''
        SUMB Class Initialization
        '''
        
        name = 'SUMB'
        category = 'Three Dimensional CFD'
        defOpts = {
            # Common Paramters
            'gridfile':[str, 'default.cgns'],
            'restartfile':[str, 'default_restart.cgns'],
            'solrestart':[bool, False],

            # Output Parameters
            'storerindlayer':[bool, True],
            'probname':[str, 'defaultName'],
            'outputdir':[str, './'],
            'writesymmetry':[bool, True],
            'writefarfield':[bool, False],
            'writesurfacesolution':[bool,True],
            'writevolumesolution':[bool,True],
            'solutionprecision':[str,'single'],
            'gridprecision':[str,'double'],
            'isosurface':[dict, {}],
            'isovariables':[list, []],
            'viscoussurfacevelocities':[bool, True],

            # Physics Paramters
            'discretization':[str, 'central plus scalar dissipation'],
            'coarsediscretization':[str, 'central plus scalar dissipation'],
            'limiter':[str, 'vanalbeda'],
            'smoother':[str, 'runge kutta'],
            'equationtype': [str, 'euler'],
            'equationmode': [str, 'steady'],
            'flowtype':[str, 'external'],
            'turbulencemodel':[str, 'sa'], 
            'turbulenceorder':[str,'first order'],           
            'usewallfunctions':[bool, False],
            'useapproxwalldistance':[bool, True],
            'reynoldsnumber':[float, 1e6], 
            'reynoldslength':[float, 1.0], 
            'walltreatment':[str, 'linear pressure extrapolation'],
            'dissipationscalingexponent':[float, 0.67],
            'vis4':[float, 0.0156],
            'vis2':[float, 0.5],
            'vis2coarse':[float, 0.5], 
            'restrictionrelaxation':[float, 1.0],

            # Common Paramters
            'ncycles':[int, 500],
            'ncyclescoarse':[int, 500],
            'nsubiterturb':[int, 1],
            'cfl':[float, 1.7],
            'cflcoarse':[float, 1.0],
            'mgcycle':[str, '3w'],
            'mgstartlevel':[int, -1],
            'resaveraging':[str,'alternateresaveraging'],
            'smoothparameter':[float, 1.5],
            'cfllimit':[float, 1.5],

            # Unsteady Paramters           
            'timeintegrationscheme':[str, 'bdf'],
            'timeaccuracy':[int, 2],
            'ntimestepscoarse':[int, 48],
            'ntimestepsfine':[int, 400],
            'deltat':[float, .010],            

            # Time Spectral Paramters
            'timeintervals': [int, 1],
            'alphamode':[bool, False],
            'betamode':[bool, False],
            'machmode':[bool, False],
            'pmode':[bool, False],
            'qmode':[bool, False],
            'rmode':[bool, False],
            'altitudemode':[bool, False],
            'windaxis':[bool, False],
            'familyrot':[str, ''],
            'rotcenter':[list, [0.0,0.0, 0.0]],
            'rotrate':[list, [0.0,0.0, 0.0]],
            'tsstability': [bool, False],

            # Convergence Paramters
            'l2convergence':[float, 1e-6],
            'l2convergencerel':[float, 1e-16],
            'l2convergencecoarse':[float, 1e-2], 
            'maxl2deviationfactor':[float, 1.0],
            'coeffconvcheck':[bool, False],
            'miniterationnum':[int, 10],

            # Newton-Krylov Paramters
            'usenksolver':[bool, False],
            'nklinearsolver':[str, 'gmres'],
            'nkswitchtol':[float, 1e-2],
            'nksubspacesize':[int, 60],
            'nklinearsolvetol':[float, 1e-1],
            'nkpc':[str, 'additive schwartz'],
            'nkasmoverlap':[int, 3],
            'nkpcilufill':[int, 3],
            'nklocalpcordering':[str, 'rcm'],
            'nkjacobianlag':[int, 10],
            'rkreset':[bool, False],
            'nrkreset':[int, 5],
            'applypcsubspacesize':[int, 10],
            'nkinnerpreconits':[int, 1],
            'nkouterpreconits':[int, 1],

            # Load Balance Paramters
            'blocksplitting':[bool, True],
            'loadimbalance':[float, 0.1],
            'loadbalanceiter':[int, 10],
       
            # Misc Paramters
            'metricconversion':[float, 1.0],
            'autosolveretry':[bool, False],
            'autoadjointretry':[bool, False],
            'storehistory':[bool, False],
            'numbersolutions':[bool, False],
            'printiterations':[bool, False],
            'printtiming':[bool, True],
            'setmonitor':[bool, True],        
            'monitorvariables':[list, ['resrho','cl', 'cd']],
            'surfacevariables':[list, ['cp','vx', 'vy','vz', 'mach']],
            'volumevariables':[list, ['resrho']],

            # Multidisciplinary Coupling Parameters:
            'forcesastractions':[bool, True],

            # Adjoint Paramters
            'adjointl2convergence':[float, 1e-10],
            'adjointl2convergencerel':[float, 1e-16],
            'adjointl2convergenceabs':[float, 1e-16],
            'adjointdivtol':[float, 1e5],
            'approxpc': [bool, False],
            'adpc': [bool, False],
            'viscpc':[bool,False],
            'usediagtspc':[bool, False],
            'restartadjoint':[bool, False],
            'adjointsolver': [str, 'gmres'],
            'adjointmaxiter': [int, 500],
            'adjointsubspacesize' : [int, 80],
            'adjointmonitorstep': [int, 10],
            'dissipationlumpingparameter':[float, 6.0],
            'preconditionerside': [str, 'right'],
            'matrixordering': [str, 'rcm'],
            'globalpreconditioner': [str, 'additive schwartz'],
            'localpreconditioner' : [str, 'ilu'],
            'ilufill': [int, 1],
            'asmoverlap' : [int, 1],
            'innerpreconits':[int,1],
            'outerpreconits':[int,1],
            'usereversemodead':[bool, True],
            'applyadjointpcsubspacesize':[int, 20],
            'frozenturbulence':[bool, True],

            # ADjoint debugger
            'firstrun':[bool, True],
            'verifystate':[bool, True],
            'verifyspatial':[bool, True],
            'verifyextra':[bool, True],
            }

        informs = {
            }

        # Load the compiled module using MExt, which allow multiple imports
        try: 
            self.sumb
        except:
            curDir = os.path.dirname(os.path.realpath(__file__))
            # Explictly only search the local directory that this file
            # resides in
            self.sumb = MExt('sumb', [curDir])._module
        # end try
        
        # Next set the MPI Communicators and associated info
        if comm is None:
            comm = MPI.COMM_WORLD

        self.comm = comm
        self.sumb.communication.sumb_comm_world = self.comm.py2f()
        self.sumb.communication.sumb_comm_self  = MPI.COMM_SELF.py2f()
        self.sumb.communication.sendrequests = numpy.zeros(self.comm.size)
        self.sumb.communication.recvrequests = numpy.zeros(self.comm.size)
        self.myid = self.sumb.communication.myid = self.comm.rank
	self.nproc = self.sumb.communication.nproc = self.comm.size

        # Initialize petec in case the user has not already
        self.sumb.initializepetsc()

        # Set the stand-alone sumb flag to flase...this changes how
        # terminate calls are handled. 
        self.sumb.iteration.standalonemode = False

        # Set the frompython flag to true... this also changes how
        # terminate calls are handled
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
             'cBend':self.sumb.costfunctions.costfuncbendingcoef,
             }
        
        self.possibleObjectives = \
            { 'lift':'Lift',
              'drag':'Drag',
              'cl':'Cl', 'cd':'Cd',
              'fx':'Fx', 'fy':'Fy','fz':'Fz', 
              'cfx':'cFx','cfy':'cFy', 'cfz':'cFz',
              'mx':'Mx', 'my':'My', 'mz':'Mz', 
              'cmx':'cMx', 'cmy':'cMy', 'cmz':'cMz', 
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
              'cbend':'cBend', 
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
            'surfaceref':'adjointvars.ndesignsurfaceref',
            'disserror':'adjointvars.ndesigndisserror'
            }

        self.aeroDVs = []

        self.optionMap = {\
            # Common Paramters
            'gridfile':{'location':'inputio.gridfile',
                        'len':self.sumb.constants.maxstringlen},
            'restartfile':{'location':'inputio.restartfile',
                           'len':self.sumb.constants.maxstringlen},
            'solrestart':{'location':'inputio.restart'},
            'storerindlayer':{'location':'inputio.storerindlayer'},
            'writesymmetry':{'location':'inputio.writesymmetry'},
            'writefarfield':{'location':'inputio.writefarfield'},
            'viscoussurfacevelocities':{'location':'inputio.viscoussurfacevelocities'},

            'solutionprecision':{'single':
                                     self.sumb.inputio.precisionsingle,
                                 'double':
                                     self.sumb.inputio.precisiondouble,
                                 'location':
                                     'inputio.precisionsol'},            
            'gridprecision':{'single':
                                 self.sumb.inputio.precisionsingle,
                             'double':
                                 self.sumb.inputio.precisiondouble,
                             'location':
                                 'inputio.precisiongrid'},            

            # Physics Paramters
            'discretization':{'central plus scalar dissipation':
                                  self.sumb.inputdiscretization.dissscalar,
                              'central plus matrix dissipation':
                                  self.sumb.inputdiscretization.dissmatrix,
                              'central plus cusp dissipation':
                                  self.sumb.inputdiscretization.disscusp,
                              'upwind':
                                  self.sumb.inputdiscretization.upwind,
                              'location':
                                  'inputdiscretization.spacediscr'},
            'coarsediscretization':{'central plus scalar dissipation':
                                        self.sumb.inputdiscretization.dissscalar,
                                    'central plus matrix dissipation':
                                        self.sumb.inputdiscretization.dissmatrix,
                                    'central plus cusp dissipation':
                                        self.sumb.inputdiscretization.disscusp,
                                    'upwind':
                                        self.sumb.inputdiscretization.upwind,
                                    'location':
                                        'inputdiscretization.spacediscrcoarse'},
            'limiter':{'vanalbeda':
                           self.sumb.inputdiscretization.vanalbeda,
                       'minmod':
                           self.sumb.inputdiscretization.minmod,
                       'nolimiter':
                           self.sumb.inputdiscretization.nolimiter,
                       'location':
                           'inputdiscretization.limiter'},
            'smoother':{'runge kutta':
                            self.sumb.inputiteration.rungekutta,
                        'lu sgs':
                            self.sumb.inputiteration.nllusgs,
                        'lu sgs line':
                            self.sumb.inputiteration.nllusgsline,
                        'dadi':
                            self.sumb.inputiteration.dadi,
                        'location':
                            'inputiteration.smoother'},
            
            'equationtype':{'euler':
                                self.sumb.inputphysics.eulerequations,
                            'laminar ns':
                                self.sumb.inputphysics.nsequations,
                            'rans':
                                self.sumb.inputphysics.ransequations,
                            'location':
                                'inputphysics.equations'},
            'equationmode':{'steady':
                                self.sumb.inputphysics.steady,
                            'unsteady':
                                self.sumb.inputphysics.unsteady,
                            'time spectral':
                                self.sumb.inputphysics.timespectral,
                            'location':
                                'inputphysics.equationmode'},
            'flowtype':{'internal':
                            self.sumb.inputphysics.internalflow,
                        'external':
                            self.sumb.inputphysics.externalflow,
                        'location':
                            'inputphysics.flowtype'},
            'turbulencemodel':{'baldwin lomax':
                                   self.sumb.inputphysics.baldwinlomax,
                               'sa':
                                   self.sumb.inputphysics.spalartallmaras,
                               'sae':
                                   self.sumb.inputphysics.spalartallmarasedwards,
                               'k omega wilcox':
                                   self.sumb.inputphysics.komegawilcox,
                               'k omega modified':
                                   self.sumb.inputphysics.komegamodified,
                               'ktau':
                                   self.sumb.inputphysics.ktau,
                               'menter sst':
                                   self.sumb.inputphysics.mentersst,
                               'v2f':
                                   self.sumb.inputphysics.v2f,
                               'location':
                                   'inputphysics.turbmodel'},
            'turbulenceorder':{'first order':1,
                               'second order':2,
                               'location':
                                   'inputdiscretization.orderturb'},
            'usewallfunctions':{'location':'inputphysics.wallfunctions'},
            'useapproxwalldistance':{
                'location':'inputdiscretization.useapproxwalldistance'},
            'reynoldsnumber':{'location':'inputphysics.reynolds'},
            'reynoldslength':{'location':'inputphysics.reynoldslength'},
            'walltreatment':{'linear pressure extrapolation':
                                 self.sumb.inputdiscretization.linextrapolpressure,
                             'constant pressure extrapolation':
                                 self.sumb.inputdiscretization.constantpressure,
                             'quadratic pressure extrapolation':
                                 self.sumb.inputdiscretization.quadextrapolpressure,
                             'normal momentum':
                                 self.sumb.inputdiscretization.normalmomentum,
                             'location':
                                 'inputdiscretization.wallbctreatment'},
            'dissipationscalingexponent':{'location':'inputdiscretization.adis'},
            'vis4':{'location':'inputdiscretization.vis4'},
            'vis2':{'location':'inputdiscretization.vis2'},
            'vis2coarse':{'location':'inputdiscretization.vis2coarse'},
            'restrictionrelaxation':{'location':'inputiteration.fcoll'},
            'forcesastractions':{'location':'inputphysics.forcesastractions'},
            
            # Common Paramters
            'ncycles':{'location':'inputiteration.ncycles'},
            'ncyclescoarse':{'location':'inputiteration.ncyclescoarse'},
            'nsubiterturb':{'location':'inputiteration.nsubiterturb'},
            'cfl':{'location':'inputiteration.cfl'},        
            'cflcoarse':{'location':'inputiteration.cflcoarse'},        
            'mgcycle':{'location':'localmg.mgdescription',
                       'len':self.sumb.constants.maxstringlen},
            'mgstartlevel':{'location':'inputiteration.mgstartlevel'},
            'resaveraging':{'noresaveraging':
                                self.sumb.inputiteration.noresaveraging,
                            'alwaysresaveraging':
                                self.sumb.inputiteration.alwaysresaveraging,
                            'alternateresaveraging':
                                self.sumb.inputiteration.alternateresaveraging,
                            'location':
                                'inputiteration.resaveraging'},
            'smoothparameter':{'location':'inputiteration.smoop'},
            'cfllimit':{'location':'inputiteration.cfllimit'},
            
            # Unsteady Params
            'timeintegrationscheme':{'bdf':
                                         self.sumb.inputunsteady.bdf,
                                     'explicitrk':
                                         self.sumb.inputunsteady.explicitrk,
                                     'inplicitrk':
                                         self.sumb.inputunsteady.implicitrk,
                                     'md':
                                         self.sumb.inputunsteady.md,
                                     'location':
                                         'inputunsteady.timeintegrationscheme'},
            'timeaccuracy':{'location':'inputunsteady.timeaccuracy'},
            'ntimestepscoarse':{'location':'inputunsteady.ntimestepscoarse'},
            'ntimestepsfine':{'location':'inputunsteady.ntimestepsfine'},
            'deltat':{'location':'inputunsteady.deltat'},
            
            # Time Spectral Paramters
            'timeintervals':{'location':'inputtimespectral.ntimeintervalsspectral'},
            'alphamode':{'location':'inputtsstabderiv.tsalphamode'},
            'betamode':{'location':'inputtsstabderiv.tsbetamode'},
            'machmode':{'location':'inputtsstabderiv.tsmachmode'},
            'pmode':{'location':'inputtsstabderiv.tspmode'},
            'qmode':{'location':'inputtsstabderiv.tsqmode'},
            'rmode':{'location':'inputtsstabderiv.tsrmode'},
            'altitudemode':{'location':'inputtsstabderiv.tsaltitudemode'},
            'windaxis':{'location':'inputtsstabderiv.usewindaxis'},
            'rotcenter':{'location':'inputmotion.rotpoint'},
            'rotrate':{'location':'inputmotion.rotrate'},
            'tsstability':{'location':'inputtsstabderiv.tsstability'},
            
            # Convergence Paramters
            'l2convergence':{'location':'inputiteration.l2conv'},
            'l2convergencerel':{'location':'inputiteration.l2convrel'},
            'l2convergencecoarse':{'location':'inputiteration.l2convcoarse'},
            'maxl2deviationfactor':{'location':'inputiteration.maxl2deviationfactor'},
            'coeffconvcheck':{'location':'monitor.coeffconvcheck'},
            'miniterationnum':{'location':'inputiteration.miniternum'},
            
            # Newton-Krylov Paramters
            'usenksolver':{'location':'nksolvervars.usenksolver'},
            'nklinearsolver':{'gmres':'gmres',
                              'tfqmr':'tfqmr',
                              'location':
                                  'nksolvervars.ksp_solver_type',
                              'len':self.sumb.constants.maxstringlen},
            
            'nkswitchtol':{'location':'nksolvervars.nk_switch_tol'},
            'nksubspacesize':{'location':'nksolvervars.ksp_subspace'},
            'nklinearsolvetol':{'location':'nksolvervars.ksp_rtol'},
            'nkpc':{'additive schwartz':'asm',
                    'multigrid':'mg',
                    'location':
                        'nksolvervars.global_pc_type',
                    'len':self.sumb.constants.maxstringlen},
            'nkasmoverlap':{'location':'nksolvervars.asm_overlap'},
            'nkpcilufill':{'location':'nksolvervars.local_pc_ilu_level'},               
            'nklocalpcordering':{'natural':'natural',
                                 'rcm':'rcm',
                                 'nested dissection':'nd',
                                 'one way dissection':'1wd',
                                 'quotient minimum degree':'qmd',
                                 'location':
                                     'nksolvervars.local_pc_ordering',
                                 'len':self.sumb.constants.maxstringlen},
            'nkmaxlinearkspits':{'location':'nksolvervars.ksp_max_it'},
            'nkjacobianlag':{'location':'nksolvervars.jacobian_lag'},
            'rkreset':{'location':'nksolvervars.rkreset'},
            'nrkreset':{'location':'nksolvervars.nrkreset'},
            'nnkfnitedifferencepc':{'location':'nksolvervars.nkfinitedifferencepc'},
            'applypcsubspacesize':{'location':'nksolvervars.applypcsubspacesize'},
            'nkinnerpreconits':{'location':'nksolvervars.innerpreconits'},
            'nkouterpreconits':{'location':'nksolvervars.outerpreconits'},

            # Load Balance Paramters
            'blocksplitting':{'location':'inputparallel.splitblocks'},
            'loadimbalance':{'location':'inputparallel.loadimbalance'},
            'loadbalanceiter':{'location':'inputparallel.loadbalanceiter'},

            # Misc Paramters
            'printiterations':{'location':'inputiteration.printiterations'},
            'printtiming':{'location':'inputadjoint.printtiming'},
            'setmonitor':{'location':'inputadjoint.setmonitor'},
            
            # Adjoint Params
            'adjointl2convergence':{'location':'inputadjoint.adjreltol'},
            'adjointl2convergencerel':{'location':'inputadjoint.adjreltolrel'},
            'adjointl2convergenceabs':{'location':'inputadjoint.adjabstol'},
            'adjointdivtol':{'location':'inputadjoint.adjdivtol'},
            'approxpc':{'location':'inputadjoint.approxpc'},
            'adpc':{'location':'inputadjoint.adpc'},
            'viscpc':{'location':'inputadjoint.viscpc'},
            'frozenturbulence':{'location':'inputadjoint.frozenturbulence'},
            'usediagtspc':{'location':'inputadjoint.usediagtspc'},
            'restartadjoint':{'location':'inputadjoint.restartadjoint'},
            'adjointsolver':{'gmres':'gmres',
                             'tfqmr':'tfqmr',
                             'richardson':'richardson',
                             'bcgs':'bcgs',
                             'ibcgs':'ibcgs',
                             'location':
                                 'inputadjoint.adjointsolvertype',
                             'len':self.sumb.constants.maxstringlen},
            'adjointmaxiter':{'location':'inputadjoint.adjmaxiter'},
            'adjointsubspacesize':{'location':'inputadjoint.adjrestart'},
            'adjointmonitorstep':{'location':'inputadjoint.adjmonstep'},
            'dissipationlumpingparameter':{'location':'inputdiscretization.sigma'},
            'preconditionerside':{'left':'left',
                                  'right':'right',
                                  'location':
                                      'inputadjoint.adjointpcside',
                                  'len':self.sumb.constants.maxstringlen},
            'matrixordering':{'natural':'natural',
                              'rcm':'rcm',
                              'nested dissection':'nd',
                              'one way dissection':'1wd',
                              'quotient minimum degree':'qmd',
                              'location':
                                  'inputadjoint.matrixordering',
                              'len':self.sumb.constants.maxstringlen},
            'globalpreconditioner':{'additive schwartz':'asm',
                                    'multigrid':'mg',
                                    'location':
                                        'inputadjoint.precondtype',
                                    'len':self.sumb.constants.maxstringlen},
            'localpreconditioner':{'ilu':'ilu',
                                   'location':
                                       'inputadjoint.localpctype',
                                   'len':self.sumb.constants.maxstringlen},
            'ilufill':{'location':'inputadjoint.filllevel'},
            'applyadjointpcsubspacesize':{
                'location':'inputadjoint.applyadjointpcsubspacesize'},
            'asmoverlap':{'location':'inputadjoint.overlap'},
            'innerpreconits':{'location':'inputadjoint.innerpreconits'},
            'outerpreconits':{'location':'inputadjoint.outerpreconits'},
            'firstrun':{'location':'inputadjoint.firstrun'},
            'verifystate':{'location':'inputadjoint.verifystate'},
            'verifyspatial':{'location':'inputadjoint.verifyspatial'},
            'verifyextra':{'location':'inputadjoint.verifyextra'},
            }                

        # These "ignore_options" are NOT actually, ignored, rather,
        # they DO NOT GET SET IN THE FORTRAN CODE. Rather, they are
        # used strictly in Python

        self.ignoreOptions = [
            'defaults',
            'storehistory',
            'numbersolutions',
            'writesurfacesolution',
            'writevolumesolution',
            'familyrot',  # -> Not sure how to do
            'autosolveretry',
            'autoadjointretry',
            'usereversemodead'
            ]

        # Deprecated options. These should not be used, but old
        # scripts can continue to run
        self.deprecatedOptions = {'finitedifferencepc':
                                      'Use the ADPC option.',
                                  'writesolution':
                                      'Use writeSurfaceSolution and writeVolumeSolution options instead.'}

        self.specialOptions = ['surfacevariables',
                               'volumevariables',
                               'monitorvariables',
                               'metricconversion',
                               'outputdir',
                               'probname',
                               'isovariables',
                               'isosurface',
                               ]

        # Info for flowCases
        self.curFlowCase = None
        self.flowCases = {}
        
        self.updateTime = 0.0
        self.nSlice = 0
        self.nLiftDist = 0

        # Set default values --- actual options will be set when
        # aero_solver is initialized
        self.sumb.setdefaultvalues()

        # If 'options' is not None, go through and make sure all keys
        # are lower case:
        if options is not None:
            for key in options.keys():
                options[key.lower()] = options.pop(key)
        else:
            options = {}
        # end if

        # Initialize the inherited aerosolver
        AeroSolver.__init__(\
            self, name, category, defOpts, informs, options=options)
        self.sumb.inputio.autoparameterupdate = False

        # Set the external Mesh Warping is provided
        self._updateGeomInfo = False

        # Sumb can be used without a external mesh warping
        # object, however, geometric sensitivities cannot be computed. 
        if mesh is None:
            self.mesh = SUmbDummyMesh()
        else:
            self.mesh = mesh
        # end if

        # Set Flags that are used to keep of track of what is "done"
        # in fortran
        self.allInitialized = False    # All flow solver initialization   

        # Matrix Setup Flags
        self.adjointSetup = False
        self._updateGeomInfo = False
        self._updatePeriodInfo = True
        self._updateVelInfo = True
        self.fatalFail = False
        self.solveFailed = False
        self.adjointFailed = False
        self.dtype = 'd'

        # Write the intro message
        self.sumb.writeintromessage()

        return

    def initialize(self, aeroProblem, partitionOnly=False):
        '''
        Run High Level Initialization 
        
        Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
        '''
        
        if self.allInitialized == True:
            return
        
        # Set periodic paramters
        self._setPeriodicParams(aeroProblem)
  
        # Make sure all the params are ok
        for option in self.options:
            if option != 'defaults':
                self.setOption(option.lower(), self.options[option][1])
            # end if
        # end for
      
        # Do the remainder of the operations that would have been done
        # had we read in a param file
        self.sumb.iteration.deforming_grid = True
        self._setMachNumber(aeroProblem)
        self._setRefState(aeroProblem)
        self._setPeriodicParams(aeroProblem)
        self.sumb.dummyreadparamfile()

        mpiPrint(' -> Partitioning and Reading Grid', comm=self.comm)
        self.sumb.partitionandreadgrid()
        if partitionOnly:
            return
        # end if

        mpiPrint(' -> Preprocessing', comm=self.comm)
        self.sumb.preprocessing()

        mpiPrint(' -> Initializing flow', comm=self.comm)
        self.sumb.initflow()
        self._setInflowAngle(aeroProblem)

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

        # Setup Surface/Force info
        npatch = self.sumb.getnpatches()
        patchnames = []
        patchsizes = []
        for i in xrange(npatch):
            tmp = numpy.zeros(256,'c')
            self.sumb.getpatchname(i+1, tmp)
            patchnames.append(
                ''.join([tmp[j] for j in range(256)]).lower().strip())
            patchsizes.append(self.sumb.getpatchsize(i+1))
        # end for
 
        conn = self.getForceConnectivity()
        pts  = self.getForcePoints()

        self.mesh.setExternalSurface(patchnames, patchsizes, conn, pts)
        
        # Solver is initialize
        self.allInitialized = True

        # 
        self.initAdjoint()

        # Setup a default flowCase 
        self.addFlowCase('default')
        self.curFlowCase = 'default'

        return

    def addFlowCase(self, flowCaseName):
        '''Add a flowCase called 'flowCaseName' to SUmb. These flow
        cases should be setup after the problem has been initialized'''

        # First check that the problem is initialized:
        if not self.allInitialized:
            mpiPrint('Error: SUmb must be initialized before\
 flowCases can be added', comm=self.comm)
            return 
        # end if

        # Special treatment of first user-supplied flowCase; overwrite
        # the default case that was automatically added:
        if self.flowCases.keys() == ['default']:
            del self.flowCases['default']

        # Add the new case
        self.flowCases[flowCaseName] = {}
        self.flowCases[flowCaseName]['states'] = None
        self.flowCases[flowCaseName]['adjoints'] = {}
        self.flowCases[flowCaseName]['surfMesh'] = \
            self.getSurfaceCoordinates('all')
        self.flowCases[flowCaseName]['aeroProblem'] = None
        self.flowCases[flowCaseName]['callCounter'] = -1
        # Only set the current flowCase name IF it is the first case
        # added. 

        if len(self.flowCases.keys()) == 1:
            self.curFlowCase = flowCaseName

        return
    
    def getFlowCase(self):
        '''Return to the user the current flowCaseName'''
        return self.curFlowCase

    def setFlowCase(self, flowCaseName, doNotSwitch=False, funcName=None):
        '''Set 'flowCaseName' in SUmb'''

        # Do nothing if the name is None. or it is the same. 
        if flowCaseName in [None, self.curFlowCase]:
            return

        # Next make sure that flowCaseName exists:
        if not flowCaseName in self.flowCases.keys():
            mpiPrint('Error: %s has not been added using \
\'addFlowCase()\''%flowCaseName, comm=self.comm)
            return 
        # end if

        # We now know the flowCaseName exists. Check that doNotSwitch
        # is not True since this is not allow. Print an error message
        # including where it happened
        if doNotSwitch:
            if funcName is not None:
                mpiPrint('Error: Trying to switch to flowCase %s from a \
function where this operation is not allowed. The offending function \
is %s'%(flowCaseName, funcName), comm=self.comm)
                sys.exit(1)
            else:
                mpiPrint('Error: Trying to switch to flowCase %s from a \
function where this operation is not allowed. The offending function \
name is unavailable.'%(flowCase), comm=self.comm)
                sys.exit(1)
            # end if
        # end if

        mpiPrint('+'+'-'*70+'+',comm=self.comm)
        mpiPrint('|  Switching to flowCase: %-45s|'%flowCaseName,comm=self.comm)
        mpiPrint('+'+'-'*70+'+',comm=self.comm)

        # Now, store the data for the current flowCase:
        self.flowCases[self.curFlowCase]['states'] = \
            self.getStates()
        self.flowCases[self.curFlowCase]['surfMesh'] = \
            self.getSurfaceCoordinates('all')

        # Set the stored data for the new flowCase:
        if self.flowCases[flowCaseName]['states'] is not None:
            self.setStates(self.flowCases[flowCaseName]['states'])
        coords = self.flowCases[flowCaseName]['surfMesh']
        if coords is not None:
            self.setSurfaceCoordinates('all', coords)
        # end if
    
        # Now we have to do a bunch of updates. This is fairly
        # expensive and flow cases should only be switched when
        # required.
        self._updateGeomInfo = True
        self._updateGeometryInfo()
        
        self.curFlowCase = flowCaseName
        self.flowCases[self.curFlowCase]['adjointRHS'] = None

        # Destroy the NK solver and the adjoint memory
        self.sumb.destroynksolver()
        self.releaseAdjointMemory()

        return

    def addLiftDistribution(self, nSegments, direction,
                            groupName=None, description=''):
        '''
        Add a lift distribution to the surface output. 
        nSegments: Number of slices to use for the distribution. Typically 150-250 is sufficient
        direction: str, one of 'x', 'y', or 'z'. Auto bases the direction
                   on the liftDir given in aeroproblem. 
        groupName: The family (as defined in pyWarp) to use for the lift distribution
        description: An additional string that can be used to destingush between 
                     multiple lift distributions in the output
        '''
        direction=direction.lower()
        if direction not in ['x','y','z']:
            mpiPrint(' Error: \'direction\' must be one of \'x\', \
\'y\', \'z\'', comm=self.comm)
            groupTag = '%s: '%groupName
            return
        else:
            groupTag = ''
        # end if

        if groupName is not None:
            mpiPrint(' Error: lift distributions by group is not yet supported')
            return
        # end if

        if direction == 'x':
            dirVec = [1.0, 0.0, 0.0]
            dirInd = 1
        elif direction == 'y':
            dirVec = [0.0, 1.0, 0.0]
            dirInd = 2
        else:
            dirVec = [0.0, 0.0, 1.0]
            dirInd = 3
        # end if
        
        distName = 'LiftDist_%2.2d %s: %s normal'%(self.nLiftDist + 1, groupTag, direction)
        self.nLiftDist += 1

        self.sumb.addliftdistribution(nSegments, dirVec, dirInd, distName)

    def addSlices(self, direction, positions, sliceType='relative', groupName=None):
        '''
        Add parametric slice positions. Slices are taken of the wing
        at time addParaSlices() is called and the parametric positions
        of the intersections on the surface mesh are stored. On
        subsequent output, the position of the slice moves as the mesh
        moves/deforms. This effectively tracks the same location on
        the wing.
        
        direction: one of 'x', 'y', 'z' 
        positions: scalar or list or array: List of slice positions 
        sliceType: One of 'relative' or 'absolute'
        groupName: The family (as defined in pyWarp) to use for the slices
        '''

        if groupName is not None:
            mpiPrint(' Error: slices by group is not yet supported')
            groupTag = '%s: '%groupName
            return
        else:
            groupTag = ''
        # end if

        direction = direction.lower()
        if direction not in ['x','y','z']:
            mpiPrint(' Error: \'direction\' must be one of \'x\', \
\'y\', \'z\'', comm=self.comm)
            return
        # end if

        sliceType = sliceType.lower()
        if sliceType not in ['relative', 'absolute']:
            mpiPrint(' Error: \'sliceType\' must be \'relative\' or \
\'absolute\'', comm=self.comm)
            return
        # end if

        positions = numpy.atleast_1d(positions)
        N = len(positions)
        tmp = numpy.zeros((N, 3),self.dtype)
        if direction == 'x':
            tmp[:, 0] = positions
            dirVec = [1.0, 0.0, 0.0]
        elif direction == 'y':
            tmp[:, 1] = positions
            dirVec = [0.0, 1.0, 0.0]
        elif direction == 'z':
            tmp[:, 2] = positions
            dirVec = [0.0, 0.0, 1.0]
        # end if

        for i in xrange(len(positions)):
            # It is important to ensure each slice get a unique
            # name...so we will number sequentially from pythhon
            j = self.nSlice + i + 1
            if sliceType == 'relative':
                sliceName = 'Slice_%4.4d %s Para Init %s=%7.3f'%(j, groupTag, direction, positions[i])
                self.sumb.addparaslice(sliceName, tmp[i], dirVec)
            else:
                sliceName = 'Slice_%4.4d %s Absolute %s=%7.3f'%(j, groupTag, direction, positions[i])
                self.sumb.addabsslice(sliceName, tmp[i], dirVec)
            # end if
        # end for
        self.nSlice += N
        return
        
    def _setInflowAngle(self, aeroProblem):
        '''
        Set the alpha and beta fromthe desiggn variables
        '''
        
        [velDir, liftDir, dragDir] = self.sumb.adjustinflowangleadjts(\
            (aeroProblem._flows.alpha*(numpy.pi/180.0)),
            (aeroProblem._flows.beta*(numpy.pi/180.0)),
            aeroProblem._flows.liftIndex)
        self.sumb.inputphysics.veldirfreestream = velDir
        self.sumb.inputphysics.liftdirection = liftDir
        self.sumb.inputphysics.dragdirection = dragDir

        if self.sumb.inputiteration.printiterations:
            mpiPrint('-> Alpha... %f %f'%(
                    numpy.real(aeroProblem._flows.alpha*(numpy.pi/180.0)),
                    numpy.real(aeroProblem._flows.alpha)), comm=self.comm)

        #update the flow vars
        self.sumb.updateflow()
        self._updateVelInfo = True

        return

    def _setElasticCenter(self, aeroProblem):
        '''
        set the value of pointRefEC for the bending moment calculation
        '''
        
        self.sumb.inputphysics.pointrefec[0] = aeroProblem._geometry.xRootec\
            *self.metricConversion
        self.sumb.inputphysics.pointrefec[1] = aeroProblem._geometry.yRootec\
            *self.metricConversion
        self.sumb.inputphysics.pointrefec[2] = aeroProblem._geometry.zRootec\
            *self.metricConversion
    
    def _setReferencePoint(self, aeroProblem):
        '''
        Set the reference point for rotations and moment calculations
        '''
        self.sumb.inputphysics.pointref[0] = aeroProblem._refs.xref\
            *self.metricConversion
        self.sumb.inputphysics.pointref[1] = aeroProblem._refs.yref\
            *self.metricConversion
        self.sumb.inputphysics.pointref[2] = aeroProblem._refs.zref\
            *self.metricConversion
        self.sumb.inputmotion.rotpoint[0] = aeroProblem._refs.xrot\
            *self.metricConversion
        self.sumb.inputmotion.rotpoint[1] = aeroProblem._refs.yrot\
            *self.metricConversion
        self.sumb.inputmotion.rotpoint[2] = aeroProblem._refs.zrot\
            *self.metricConversion
        #update the flow vars
        self.sumb.updatereferencepoint()
        self._updateVelInfo = True

        return

    def _setRotationRate(self, aeroProblem):
        '''
        Set the rotational rate for the grid
        '''
        a  = numpy.sqrt(self.sumb.flowvarrefstate.gammainf*\
                      self.sumb.flowvarrefstate.pinfdim/ \
                      self.sumb.flowvarrefstate.rhoinfdim)
        V = (self.sumb.inputphysics.machgrid+self.sumb.inputphysics.mach)*a
        
        p = aeroProblem._flows.phat*V/aeroProblem._refs.bref
        q = aeroProblem._flows.qhat*2*V/aeroProblem._refs.cref
        r = aeroProblem._flows.rhat*V/aeroProblem._refs.bref

        self.sumb.updaterotationrate(p, r, q)
        self._updateVelInfo = True

        return
    
    def _setRefArea(self, aeroProblem):
        self.sumb.inputphysics.surfaceref = aeroProblem._refs.sref*self.metricConversion**2
        self.sumb.inputphysics.lengthref = aeroProblem._refs.cref*self.metricConversion
        
        return

    def _setPeriodicParams(self, aeroProblem):
        '''
        Set the frequecy and amplitude of the oscillations
        '''
        if  self.getOption('alphaMode'):
            self.sumb.inputmotion.degreepolalpha = int(aeroProblem._flows.degreePol)
            self.sumb.inputmotion.coefpolalpha = aeroProblem._flows.coefPol
            self.sumb.inputmotion.omegafouralpha   = aeroProblem._flows.omegaFourier
            self.sumb.inputmotion.degreefouralpha  = aeroProblem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffouralpha = aeroProblem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffouralpha = aeroProblem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif  self.getOption('betaMode'):
            self.sumb.inputmotion.degreepolmach = int(aeroProblem._flows.degreePol)
            self.sumb.inputmotion.coefpolmach = aeroProblem._flows.coefPol
            self.sumb.inputmotion.omegafourbeta   = aeroProblem._flows.omegaFourier
            self.sumb.inputmotion.degreefourbeta  = aeroProblem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourbeta = aeroProblem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourbeta = aeroProblem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('machMode'):
            self.sumb.inputmotion.degreepolmach = int(aeroProblem._flows.degreePol)
            self.sumb.inputmotion.coefpolmach = aeroProblem._flows.coefPol
            self.sumb.inputmotion.omegafourmach   = aeroProblem._flows.omegaFourier
            self.sumb.inputmotion.degreefourmach  = aeroProblem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourmach = aeroProblem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourmach = aeroProblem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif  self.getOption('pMode'):
            ### add in lift axis dependence
            self.sumb.inputmotion.degreepolxrot = int(aeroProblem._flows.degreePol)
            self.sumb.inputmotion.coefpolxrot = aeroProblem._flows.coefPol
            self.sumb.inputmotion.omegafourxrot = aeroProblem._flows.omegaFourier
            self.sumb.inputmotion.degreefourxrot  = aeroProblem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourxrot = aeroProblem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourxrot = aeroProblem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('qMode'):
            self.sumb.inputmotion.degreepolzrot = int(aeroProblem._flows.degreePol)
            self.sumb.inputmotion.coefpolzrot = aeroProblem._flows.coefPol
            self.sumb.inputmotion.omegafourzrot = aeroProblem._flows.omegaFourier
            self.sumb.inputmotion.degreefourzrot  = aeroProblem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffourzrot = aeroProblem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffourzrot = aeroProblem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('rMode'):
            self.sumb.inputmotion.degreepolyrot = int(aeroProblem._flows.degreePol)
            self.sumb.inputmotion.coefpolyrot = aeroProblem._flows.coefPol
            self.sumb.inputmotion.omegafouryrot = aeroProblem._flows.omegaFourier
            self.sumb.inputmotion.degreefouryrot  = aeroProblem._flows.degreeFourier
            self.sumb.inputmotion.coscoeffouryrot = aeroProblem._flows.cosCoefFourier
            self.sumb.inputmotion.sincoeffouryrot = aeroProblem._flows.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        # end if

        self._update_period_info = True
        self._update_geom_info = True
        self._update_vel_info = True
 
        return

    def _setMachNumber(self, aeroProblem):
        '''
        Set the mach number for the problem...
        '''
        if self.getOption('familyRot') != '':
            Rotating = True
        else:
            Rotating = False
        #endif

        if Rotating or self.getOption('equationMode').lower()=='time spectral':
            self.sumb.inputphysics.mach = 0.0
            self.sumb.inputphysics.machcoef = aeroProblem._flows.mach
            self.sumb.inputphysics.machgrid = aeroProblem._flows.mach
            self.sumb.inputmotion.gridmotionspecified = True
        else:
            self.sumb.inputphysics.mach = aeroProblem._flows.mach 
            self.sumb.inputphysics.machcoef = aeroProblem._flows.mach
            self.sumb.inputphysics.machgrid = 0.0
        # end if

        return

    def _setRefState(self, aeroProblem):
        ''' Set the Pressure, density and viscosity/reynolds number
        from the aeroProblem
        '''
        self.sumb.flowvarrefstate.pref = aeroProblem._flows.P
        self.sumb.flowvarrefstate.rhoref = aeroProblem._flows.rho
        self.sumb.flowvarrefstate.tref = aeroProblem._flows.T

        # Reynolds number info not setup yet...

        return

    def _updatePeriodicInfo(self):
        """Update the SUmb TS period info"""
        if (self._update_period_info):
            self.sumb.updateperiodicinfoalllevels()
            self._updatePeriodInfo = False
        # end if

        return 

    def _updateVelocityInfo(self):
        if (self._update_vel_info):
            self.sumb.updategridvelocitiesalllevels()
            self._updateVelInfo = False
        # end if
        
        return 
    
    def resetAdjoint(self, obj, flowCase=None):
        '''
        Reset a possible stored adjoint 'obj'
        '''
        self.setFlowCase(flowCase)

        if obj in self.flowCases[self.curFlowCase]['adjoints'].keys():
            self.flowCases[self.curFlowCase]['adjoints'][obj][:] = 0.0
        # end if

        return

    def resetFlow(self, aeroProblem=None, flowCase=None):
        '''
        Reset the flow for the complex derivative calculation
        '''
        self.setFlowCase(flowCase)

        if aeroProblem is not None:
            self._setInflowAngle(aeroProblem)
            self._setMachNumber(aeroProblem)
            self._setRefState(aeroProblem)
            self.sumb.referencestate()
            self.sumb.setflowinfinitystate()
            if self.myid == 0:
                print ('alpha:',aeroProblem._flows.alpha*(numpy.pi/180))
                print ('Mach:',aeroProblem._flows.mach)

        #mgLvlSave =  self.sumb.inputiteration.mgstartlevel
        #self.sumb.inputiteration.mgstartlevel = 1
        strLvl =  self.getOption('MGStartLevel')
        nLevels = self.sumb.inputiteration.nmglevels
        if strLvl < 0 or strLvl > nLevels :
            strLvl = nLevels
        # end if
        self.sumb.inputiteration.mgstartlevel = strLvl
        self.sumb.inputiteration.groundlevel = strLvl
        self.sumb.inputiteration.currentlevel = strLvl
        self.sumb.monitor.niterold = 0
        self.sumb.monitor.nitercur = 0
        self.sumb.iteration.itertot = 0
        self.sumb.setuniformflow()
        self.sumb.nksolvervars.nksolvecount = 0
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
        a  = numpy.sqrt(self.sumb.flowvarrefstate.gammainf*\
                      self.sumb.flowvarrefstate.pinfdim/ \
                      self.sumb.flowvarrefstate.rhoinfdim)
        V = (self.sumb.inputphysics.machgrid+self.sumb.inputphysics.mach)*a

        return V

    def __solve__(self, aeroProblem, nIterations=500, flowCase=None, 
                  MDCallBack=None, writeSolution=True):
        
        '''
        Run Analyzer (Analyzer Specific Routine)
        
        Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
        '''
        # Set the desired flow Case
        self.setFlowCase(flowCase)
        
        # Possibly release adjoint memory (if flowCase is unchanged,
        # this would not have been  done in the above caall
        self.releaseAdjointMemory()

        # Save aeroProblem, and other information into the current flow case
        self.flowCases[self.curFlowCase]['aeroProblem'] = aeroProblem
        self.flowCases[self.curFlowCase]['adjointRHS'] = None
        self.flowCases[self.curFlowCase]['callCounter'] += 1

        # Run Initialize, if already run it just returns.
        self.initialize(aeroProblem)

        # Set all the infor contained in the aeroProblem object
        self._setMachNumber(aeroProblem)
        self._setPeriodicParams(aeroProblem)
        self._setInflowAngle(aeroProblem)
        self._setReferencePoint(aeroProblem)
        #self._setElasticCenter(aeroProblem)
        self._setRotationRate(aeroProblem)
        self._setRefArea(aeroProblem)
        self._setRefState(aeroProblem)

        # Run Solver
        t0 = time.time()

        # set the number of cycles for this call
        self.sumb.inputiteration.ncycles = nIterations

        # Cold Start:
        if self.sumb.monitor.niterold == 0 and \
            self.sumb.monitor.nitercur == 0 and \
            self.sumb.iteration.itertot == 0:
            if self.myid == 0:
                desiredSize = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles
                self.sumb.allocconvarrays(desiredSize)
            # end if
        else:
            # More Time Steps / Iterations OR a restart
            # Reallocate convergence history array and time array
            # with new size,  storing old values from previous runs
            if self.getOption('storeHistory'):
                currentSize = len(self.sumb.monitor.convarray)
                desiredSize = currentSize + self.sumb.inputiteration.ncycles+1
                self.sumb.monitor.niterold  = self.sumb.monitor.nitercur+1
            else:
                self.sumb.monitor.nitercur  = 0
                self.sumb.monitor.niterold  = 1
                desiredSize = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles +1
            # end if
            # Allocate Arrays
            if self.myid == 0:
                self.sumb.allocconvarrays(desiredSize)

            self.sumb.inputiteration.mgstartlevel = 1
            self.sumb.iteration.itertot = 0
        # end if

        if self.getOption('equationMode') == 'unsteady':
            self.sumb.alloctimearrays(self.getOption('nTimeStepsFine'))

        # Reset Fail Flags
        self.sumb.killsignals.routinefailed =  False
        self.sumb.killsignals.fatalfail = False
        self.solveFailed =  self.fatalFail = False

        self._updatePeriodicInfo()
        if (self.getOption('equationMode').lower() == 'steady' or 
            self.getOption('equationMode').lower() == 'time spectral'):
            self._updateGeometryInfo()
        self._updateVelocityInfo()

        # Check to see if the above update routines failed.
        self.sumb.killsignals.routinefailed = \
            self.comm.allreduce(
            bool(self.sumb.killsignals.routinefailed), op=MPI.LOR)

        if self.sumb.killsignals.routinefailed:
            mpiPrint('Fatal failure during mesh warp', comm=self.comm)
            self.fatalFail = True
            self.solveFailed = True
            return
     
        t1 = time.time()

        # Call the Solver or the MD callback solver
        if MDCallBack is None:
            self.sumb.solver()
        else:
            self.sumb.solverunsteadymd(MDCallBack)
        # end if

        # Save the states into the flowCase
        self.flowCases[self.curFlowCase]['states'] = \
            self.getStates()
            
        # Assign Fail Flags
        self.solveFailed = self.sumb.killsignals.routinefailed
        self.fatalFail = self.sumb.killsignals.fatalfail

        # Reset Flow if there's a fatal fail reset and return;
        # --> Do not write solution
        if self.fatalFail:
            self.resetFlow()
            return
        # end if

        t2 = time.time()
        solTime = t2 - t0
        self.updateTime = t1 - t0

        if self.getOption('printTiming'):
            mpiPrint('Solution Time: %10.3f sec'%solTime, comm=self.comm)
        # end if

        # Post-Processing -- Write Solutions is requested
        if writeSolution:
            self.writeSolution()
                    
        if self.getOption('TSStability'):
            self.computeStabilityParameters()
        # end if
        
        return

    def solveCL(self, aeroProblem, CLStar, nIterations=500, alpha0=0, 
                delta=0.5, tol=1e-3, autoReset=True, flowCase=None):
        '''This is a simple secant method search for solving for a
        fixed CL. This really should only be used to determine the
        starting alpha for a lift constraint in an optimization.

        Input:  aeroProblem -> aerodynamic problem definition
                CLStar     -> Target CL
                nIterations -> Number of CFD iterations to run (same as 
                               input to __solve__
                alpha0      -> Initial guess for secant search (deg)
                delta       -> Initial step direction (deg)
                tol         -> Absolute tolerance for CL convergence
        Output: aeroProblem._flows.alpha is updated with correct alpha
        '''
        self.setFlowCase(flowCase)

        anm2 = alpha0
        anm1 = alpha0 + delta

        # Solve for the n-2 value:
        aeroProblem._flows.alpha = anm2
        self.__solve__(aeroProblem, nIterations=nIterations)
        sol = self.getSolution()
        fnm2 =  sol['cl'] - CLStar

        L2ConvSave = self.getOption('l2convergence')
        for iIter in xrange(20):
            # We need to reset the flow since changing the alpha leads
            # to problems with the NK solver
            if autoReset:
                self.resetFlow()

            # Sometimes with the RKSolver, the residual doesn't spike
            # immediately due to the alpha, and the solver will
            # convergnce after 2 iterations. We slightly lower the
            # tolerance at each iteration to prevent this. 
            self.setOption('l2convergence', 0.95*self.getOption('l2convergence'))

            # Set current alpha
            aeroProblem._flows.alpha = anm1

            # Solve for n-1 value (anm1)
            self.__solve__(aeroProblem, nIterations=nIterations, writeSolution=False)
            sol = self.getSolution()
            fnm1 =  sol['cl'] - CLStar
            
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
        
        # Restore the initial tolerance so the user isn't confused
        self.setOption('l2convergence', L2ConvSave)
        return

    def getSurfaceCoordinates(self, groupName):
        ''' 
        See MultiBlockMesh.py for more info
        '''

        return self.mesh.getSurfaceCoordinates(groupName)

    def setSurfaceCoordinates(self, groupName, coordinates, flowCase=None):
        ''' 
        See MultiBlockMesh.py for more info
        '''

        self._updateGeomInfo = True
        self.mesh.setSurfaceCoordinates(groupName, coordinates)
        
        if flowCase is None:
            for case in self.flowCases.keys():
                self.flowCases[case]['surfMesh'] = \
                    coordinates.copy()
            # end for
        else:
            self.flowCases[flowCase]['surfMesh'] = \
                coordinates.copy()
        # end if

        return 

    def getSurfaceConnectivity(self, groupName):
        '''
        See MultiBlockMesh.py for more info
        '''
        return self.mesh.getSurfaceConnectivity(groupName)
        
    def writeSolution(self, outputDir=None, baseName=None, number=None):
        '''This is a generic shell function that potentially writes
        the various output files. The intent is that the user or
        calling program can call this file and SUmb write all the
        files that the user has defined. It is recommneded that this
        function is used along with the associated logical flags in
        the options to determine the desired writing procedure

        Optional arguments 

        outputDir: Use the supplied output directory

        baseName: Use this supplied string for the base filename. Typically
                  only used from an external solver.
        number: Use the user spplied number to index solutino. Again, only
                typically used from an external solver.
                '''
        if outputDir is None:
            outputDir = self.getOption('outputDir')

        if baseName is None:
            baseName = self.getOption('probName')

        # Join to get the base filename
        base = os.path.join(outputDir, baseName)

        # If we have flow cases, add the flow case name:
        if self.curFlowCase <> "default":
            base = base + '_%s'%self.curFlowCase
        # end if

        # If we are numbering solution, it saving the sequence of
        # calls, add the call number
        if number is not None:
            # We need number based on the provided number:
            base = base + '_%3.3d'%number
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
            if self.getOption('numberSolutions'):            
                base = base + '_%3.3d'%self.flowCases[self.curFlowCase]['callCounter']
            # end if
        # end if
        
        # Now call each of the 4 routines with the appropriate file name:
        if self.getOption('writevolumesolution'):
            self.writeVolumeSolutionFile(base + '_vol.cgns')
        if self.getOption('writesurfacesolution'):
            self.writeSurfaceSolutionFile(base + '_surf.cgns')

        self.writeLiftDistributionFile(base + '_lift.dat')
        self.writeSlicesFile(base + '_slices.dat')
        
        return

    def writeMeshFile(self, fileName):
        """Write the current mesh to a CGNS file. This call isn't used
                normally since the volume solution usually contains the grid
        
        fileName -- the name of the file 
        """
        
        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        if ext <> 'cgns':
            fileName += '.cgns'
        # end if

        # Set Flags for writing
        self.sumb.monitor.writegrid = True
        self.sumb.monitor.writevolume = False
        self.sumb.monitor.writesurface = False

        # Set filename in sumb
        self.sumb.inputio.solfile[:] = ''
        self.sumb.inputio.solfile[0:len(filename)] = filename
        
        self.sumb.inputio.newgridfile[:] = ''
        self.sumb.inputio.newgridfile[0:len(filename)] = filename

        # Actual fortran write call
        self.sumb.writesol()

        return

    def writeVolumeSolutionFile(self, fileName, writeGrid=True):
        """Write the current state of the volume flow solution to a CGNS file.
                Keyword arguments:
        
        filename -- the name of the file 
        writeGrid -- Include the grid or use links. Always 
                     writing the grid is recommended, even in cases
                     where it is not strictly necessary

        """
        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        if ext <> 'cgns':
            fileName += '.cgns'
        # end if

        # Set Flags for writing
        self.sumb.monitor.writegrid = writeGrid
        self.sumb.monitor.writevolume = True
        self.sumb.monitor.writesurface = False

        # Set filename in sumb
        self.sumb.inputio.solfile[:] = ''
        self.sumb.inputio.solfile[0:len(fileName)] = fileName

        self.sumb.inputio.newgridfile[:] = ''
        self.sumb.inputio.newgridfile[0:len(fileName)] = fileName

        # Actual fortran write call
        self.sumb.writesol()

        return

    def writeSurfaceSolutionFile(self, fileName):
        '''Write the current state of the surface flow solution to a CGNS file.
        Keyword arguments:
        fileName -- the name of the file
        '''
        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        if ext <> 'cgns':
            fileName += '.cgns'
        # end if

        # Set Flags for writing
        self.sumb.monitor.writegrid=False
        self.sumb.monitor.writevolume=False
        self.sumb.monitor.writesurface=True

        # Set filename in sumb
        self.sumb.inputio.surfacesolfile[:] = ''
        self.sumb.inputio.surfacesolfile[0:len(fileName)] = fileName

        # Actual fortran write call
        self.sumb.writesol()

        return

    def writeLiftDistributionFile(self, fileName):
        '''Evaluate and write the lift distibution to a tecplot file.
        fileName: Filename of the lift file.
        '''
        
        # Ensure filename is .dat even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        if ext <> 'dat':
            fileName += '.dat'
        # end if
            
        # Actual write command
        self.sumb.writeliftdistributionfile(fileName)

        return

    def writeSlicesFile(self, fileName):
        '''Evaluate and write the defined slice information to a
        tecplot file.
        fileName: Filename of the slice file
        '''

        # Ensure filename is .dat even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        if ext <> 'dat':
            fileName += '.dat'
        # end if
            
        # Actual write command
        self.sumb.writeslicesfile(fileName)

        return

    def getTriangulatedMeshSurface(self, flowCase=None):
        '''
        This function returns a trianguled verision of the surface
        mesh on all processors. The intent is to use this for doing
        constraints in DVConstraints.
        '''

        self.setFlowCase(flowCase, True, 'getTriangulatedMeshSurface')

        # Use first spectral instance
        pts = self.comm.allgather(self.getForcePoints(0))
        conn = self.comm.allgather(self.mesh.getSurfaceConnectivity('all'))

        # Triangle info...point and two vectors 
        p0 = []
        v1 = []
        v2 = []

        for iProc in xrange(len(conn)):
            for i in xrange(len(conn[iProc])/4):
                i0 = conn[iProc][4*i+0]
                i1 = conn[iProc][4*i+1]
                i2 = conn[iProc][4*i+2]
                i3 = conn[iProc][4*i+3]

                p0.append(pts[iProc][i0])
                v1.append(pts[iProc][i1]-pts[iProc][i0])
                v2.append(pts[iProc][i3]-pts[iProc][i0])

                p0.append(pts[iProc][i2])
                v1.append(pts[iProc][i1]-pts[iProc][i2])
                v2.append(pts[iProc][i3]-pts[iProc][i2])

            # end for
        # end for
                
        return [p0, v1, v2]

    def writeForceFile(self, fileName, TS=0, groupName='all', 
                       cfdForcePts=None, flowCase=None):
        '''This function collects all the forces and locations and
        writes them to a file with each line having: X Y Z Fx Fy Fz.
        This can then be used to set a set of structural loads in TACS
        for structural only optimization 

        Like the getForces() routine, an external set of forces may be
        passed in on which to evaluate the forces. This is only
        typically used in an aerostructural case. 

        '''
        self.setFlowCase(flowCase, True, 'writeForceFile')

        if self.mesh is None:
            mpiPrint('Error: A pyWarp mesh be specified to use writeForceFile',
                     comm=self.comm)
            return

        # Now we need to gather the data:
        if cfdForcePts is None:
            pts = self.comm.gather(self.getForcePoints(TS), root=0)
        else:
            pts = self.comm.gather(cfdForcePts)
        # end if
            
        # Forces are still evaluated on the displaced surface so do NOT pass in pts.
        forces = self.comm.gather(self.getForces(groupName, TS=TS), root=0)
        conn   = self.comm.gather(self.mesh.getSurfaceConnectivity(groupName),
                                  root=0)

        # Write out Data only on root proc:
        if self.myid == 0:
            # First sum up the total number of nodes and elements:
            nPt = 0
            nCell = 0
            for iProc in xrange(len(pts)):
                nPt += len(pts[iProc])
                nCell += len(conn[iProc])/4
            # end for
     
            # Open output file
            f = open(fileName, 'w')
            
            # Write header with number of nodes and number of cells
            f.write("%d %d\n"%(nPt, nCell))

            # Now write out all the Nodes and Forces (or tractions)
            for iProc in xrange(len(pts)):
                for i in xrange(len(pts[iProc])):
                    f.write('%15.8g %15.8g %15.8g '%(
                            numpy.real(pts[iProc][i,0]),
                            numpy.real(pts[iProc][i,1]),
                            numpy.real(pts[iProc][i,2])))
                    f.write('%15.8g %15.8g %15.8g\n'%(
                            numpy.real(forces[iProc][i,0]),
                            numpy.real(forces[iProc][i,1]),


                            numpy.real(forces[iProc][i,2])))
                # end for
            # end for

            # Now write out the connectivity information. We have to
            # be a little careful, since the connectivitiy is given
            # locally per proc. As we loop over the data from each
            # proc, we need to increment the connectivity by the numer
            # of nodes we've "used up" so far

            nodeOffset = 0
            for iProc in xrange(len(conn)):
                for i in xrange(len(conn[iProc])/4):
                    f.write('%d %d %d %d\n'%(
                            conn[iProc][4*i+0]+nodeOffset,
                            conn[iProc][4*i+1]+nodeOffset,
                            conn[iProc][4*i+2]+nodeOffset,
                            conn[iProc][4*i+3]+nodeOffset))
                # end for
                nodeOffset += len(pts[iProc])
            # end for

            f.close()
        # end if (root proc )

        return 

    def getForces(self, groupName=None, TS=0, pressure=True, viscous=True, 
                  flowCase=None):
        ''' Return the forces on this processor.
        '''
        self.setFlowCase(flowCase, True, 'getForces')
        [npts, ncell] = self.sumb.getforcesize()

        if npts > 0:
            forcesp, forcesv = self.sumb.getforces(npts, TS+1)
            forcesp = forcesp.T
            forcesv = forcesv.T

            forces = numpy.zeros_like(forcesp)
            if pressure:
                forces += forcesp
            if viscous:
                forces += forcesv
        else:
            forces = numpy.zeros((0,3),self.dtype)
        # end if

        if groupName is not None:
            forces = self.mesh.sectionVectorByFamily(groupName, forces)

        return forces

    def getForcePoints(self, TS=0, flowCase=None):
        self.setFlowCase(flowCase, True, 'getForcePoints')

        [npts, ncell] = self.sumb.getforcesize()
        pts = numpy.zeros((npts, 3),self.dtype)
        if npts > 0:
            self.sumb.getforcepoints(pts.T, TS+1)
        
        return pts

    def getForceConnectivity(self):
        [npts, ncells] = self.sumb.getforcesize()
        conn =  numpy.zeros((ncells, 4), dtype='intc')
        self.sumb.getforceconnectivity(numpy.ravel(conn))

        return conn

    def verifyBendingPartial(self):
        self.sumb.verifybendingderivatives()
        
        return

    def globalNKPreCon(self, inVec, outVec, flowCase=None):
        '''This function is ONLY used as a preconditioner to the
        global Aero-Structural system'''
        self.setFlowCase(flowCase, True, 'globalNKPreCon')
        outVec = self.sumb.applypc(inVec, outVec)
        
        return outVec

    def globalAdjointPreCon(self, inVec, outVec, flowCase=None):
        ''' This function is ONLY used as a preconditioner for the
        global Aero-Structural Adjoint system'''
        self.setFlowCase(flowCase, True, 'globalAdjointPreCon')

        outVec = self.sumb.applyadjointpc(inVec, outVec)

        return outVec

    def verifyAD(self):
        '''
        Use Tapenade TGT debugger to verify AD
        '''
        self.sumb.verifyad()

        return

    def initAdjoint(self):
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
        self.sumb.adjointvars.ndesigndisserror = -1
        
        # Set the required paramters for the aero-Only design vars:
        self.nDVAero = len(self.aeroDVs)#for debuggin with check all...

        self.sumb.adjointvars.ndesignextra = self.nDVAero
        
        if self.nDVAero >0:           
            self.sumb.adjointvars.dida = numpy.zeros(self.nDVAero)
            for i in xrange(self.nDVAero):
                execStr = 'self.sumb.' + self.possibleAeroDVs[self.aeroDVs[i]] + \
                           '= %d'%(i)
                # Leave this zero-based since we only need to use it in petsc
                exec(execStr)
            # end for
        # end if

        # Run the small amount of fortran code for adjoint initialization
        self.sumb.preprocessingadjoint()

        return

    def addAeroDV(self, *dvs):
        '''Take in a list of DVs that the flow solver will use in
        addition to shape-type design variables'''
        for i in xrange(len(dvs)):
            if dvs[i] in self.possibleAeroDVs and not dvs[i] in self.aeroDVs:
                self.aeroDVs.append(dvs[i])
            else:
                print('Warning: %s was not one of the possible AeroDVs'%(dvs[i]))
                print('Full AeroDV list is:')
                print(self.possibleAeroDVs)
            # end if
        # end for

        # Need to run initAdjoint() to update the aeroDVs in fortran
        self.initAdjoint()

        return

    def setupAdjoint(self, reform=False, flowCase=None):
        '''
        Setup the data structures required to solve the adjoint problem
        '''

        # Set the flow Case
        self.setFlowCase(flowCase)

        # Destroy the NKsolver to free memory -- Call this even if the
        # solver is not used...a safeguard check is done in Fortran
        self.sumb.destroynksolver()

        # For now, just create all the petsc variables
        if not self.adjointSetup or reform:
            self.sumb.createpetscvars()

            if self.getOption('useReverseModeAD'):
                self.sumb.setupallresidualmatrices()
            else:
                self.sumb.setupallresidualmatricesfwd()
            # end if
                
            # Create coupling matrix struct whether we need it or not
            [npts, ncells] = self.sumb.getforcesize()
            nTS  = self.sumb.inputtimespectral.ntimeintervalsspectral
            forcePoints = numpy.zeros((nTS, npts, 3),self.dtype)
            for i in xrange(nTS):
                forcePoints[i] = self.getForcePoints(TS=i)

            self.sumb.setupcouplingmatrixstruct(forcePoints.T)
            
            # Setup the KSP object
            self.sumb.setuppetscksp()

            # Set the flag
            self.adjointSetup = True
        # end if

        return

    def printMatrixInfo(self, dRdwT=False, dRdwPre=False, dRdx=False,
                        dRda=False, dSdw=False, dSdx=False,
                        printLocal=False, printSum=False, printMax=False):
        
        # Call sumb matrixinfo function
        self.sumb.matrixinfo(dRdwT, dRdwPre, dRdx, dRda, dSdw, dSdx, 
                             printLocal, printSum, printMax)

        return

    def checkPartitioning(self, nprocs):
        '''This function determine the potential load balancing for
        nprocs. The intent is this function can be run in serial
        to determine the best number of procs for load balancing. The
        grid is never actually loaded so this function can be run with
        VERY large grids without issue.'''
  
        loadInbalance, faceInbalance = self.sumb.checkpartitioning(nprocs)
                
        return loadInbalance, faceInbalance
    
    def releaseAdjointMemory(self):
        '''
        release the PETSc Memory that have been allocated
        '''
        if self.adjointSetup:
            self.sumb.destroypetscvars()
            self.adjointSetup = False

        return

    def _on_adjoint(self, objective, forcePoints=None, structAdjoint=None, 
                    groupName=None, flowCase=None):

        # Try to see if obj is an aerodynamic objective. If it is, we
        # will have a non-zero RHS, otherwise its an objective with a
        # zero aerodynamic RHS
       
        self.setFlowCase(flowCase)

        # We need to reset some of the flow condition because the flow
        # case may have changed.
        aeroProblem = self.flowCases[self.curFlowCase]['aeroProblem']
        self._setMachNumber(aeroProblem)
        self._setPeriodicParams(aeroProblem)
        self._setInflowAngle(aeroProblem)
        self._setReferencePoint(aeroProblem)
        self._setRotationRate(aeroProblem)
        self._setRefArea(aeroProblem)
        self._setRefState(aeroProblem)

        obj, aeroObj = self._getObjective(objective)

        # Setup adjoint matrices/vector as required
        self.setupAdjoint()

        # Check to see if the RHS Partials have been computed
        if not self.flowCases[self.curFlowCase]['adjointRHS'] == obj:
            self.computeObjPartials(objective, forcePoints)
        # end if

        # Check to see if we need to agument the RHS with a structural
        # adjoint:
        if structAdjoint is not None and groupName is not None:
            if self.getOption('usereversemodead'):
                print('Reverse mode AD no longer supported with \
aerostructural analysis. Use Forward mode AD for the adjoint')
                sys.exit(0)
            # end if

            phi = self.mesh.expandVectorByFamily(groupName, structAdjoint)
            self.sumb.agumentrhs(numpy.ravel(phi))
        # end if

        # Check if objective is allocated:
        if obj not in self.flowCases[self.curFlowCase]['adjoints'].keys():
            self.flowCases[self.curFlowCase]['adjoints'][obj] = \
                numpy.zeros(self.getAdjointStateSize(), float)
        self.sumb.setadjoint(self.flowCases[self.curFlowCase]['adjoints'][obj])

        # Actually Solve the adjoint system
        self.sumb.solveadjointtransposepetsc()

        # Possibly try another solve
        if self.sumb.killsignals.adjointfailed and self.getOption('restartAdjoint'):
            # Only retry if the following conditions are met:

            # 1. restartAdjoint is true -> that is we were starting
            # from a non-zero starting point

            # 2. The stored adjoint must have been already set at
            # least once; that is we've already tried one solve
            
            self.flowCases[self.curFlowCase]['adjoints'][obj][:] = 0.0
            if self.getOption('autoAdjointRetry'):
                self.sumb.solveadjointtransposepetsc()
            # end if
        # end if

        # Now set the flags and possibly reset adjoint
        if self.sumb.killsignals.adjointfailed == False:
            self.flowCases[self.curFlowCase]['adjoints'][obj] = \
                self.sumb.getadjoint(self.getAdjointStateSize())
            self.adjointFailed = False
        else:
            self.adjointFailed = True

            # Reset stored adjoint
            self.flowCases[self.curFlowCase]['adjoints'][obj][:] = 0.0
        # end if
       
        return

    def totalSurfaceDerivative(self, objective, flowCase=None):
        # The adjoint vector is now calculated so perform the
        # following operation to produce dI/dX_surf:
        # (p represents partial, d total)
        # dI/dX_s = pI/pX_s - (dXv/dXs)^T * ( dRdX_v^T * psi)
        # 
        # The derivative wrt the surface captures the effect of ALL
        # GLOBAL Multidisciplinary variables -- any DV that changes
        # the surface. 
        self.setFlowCase(flowCase, True, 'totalSurfaceDerivative')

        obj, aeroObj = self._getObjective(objective)

        # NOTE: do dRdxvPsi MUST be done first since this
        # allocates spatial memory if required.
        dIdxs_2 = self.getdRdXvPsi('all', objective)
          
        # Direct partial derivative contibution 
        dIdxs_1 = self.getdIdx(objective, groupName='all')

        # Total derivative of the obective with surface coordinates
        dIdXs = dIdxs_1 - dIdxs_2

        return dIdXs

    def totalAeroDerivative(self, objective, flowCase=None):
        # The adjoint vector is now calculated. This function as above
        # computes dI/dX_aero = pI/pX_aero - dR/dX_aero^T * psi. The
        # "aero" variables are intrinsic ONLY to the aero
        # discipline. Nothing in the structural process should depend
        # on these functions directly. 
        self.setFlowCase(flowCase, True, 'totalAeroDerivative')

        obj, aeroObj = self._getObjective(objective)

        if obj in self.flowCases[self.curFlowCase]['adjoints'].keys():
            psi = self.flowCases[self.curFlowCase]['adjoints'][obj]
        else:
            mpiPrint('%s adjoint for flowCase %s is not computed.'%(
                    obj, self.curFlowCase), comm=self.comm)
            sys.exit(1)
        # end if

        # Direct partial derivative contibution 
        dIda_1 = self.getdIda(objective)

        # dIda contribution for drda^T * psi
        dIda_2 = self.getdRdaPsi(psi)

        # Total derivative of the obective wrt aero-only DVs
        dIda = dIda_1 - dIda_2

        return dIda

    def saveAdjointMatrix(self, fileName):
        ''' Save the adjoint matrix to a binary petsc file for
        possible future resue'''
        if self.adjointSetup:
            self.sumb.saveadjointmatrix(fileName)
        else:
            mpiPrint('Cannot save matrix since adjoint not setup.',
                     comm=self.comm)
        # end if

        return

    def saveAdjointPC(self, fileName):
        ''' Save the adjoint preconditioning matrix to a binary petsc
        file for possible future resue'''
        if self.adjointSetup and self.getOption('approxpc'):
            self.sumb.saveadjointpc(fileName)
        else:
            mpiPrint('Cannot save PC matrix since adjoint not setup.',
                     comm=self.comm)
        # end if

        return

    def saveAdjointRHS(self, fileName):
        ''' Save the current adjoint RHS to a binary petsc file for
        possible future resue'''
        ''' Save the adjoint matrix to a binary petsc file for
        possible future resue'''
        if self.adjointSetup:
            self.sumb.saveadjointrhs(fileName)
        else:
            mpiPrint('Cannot save RHS since adjoint not setup.',
                     comm=self.comm)
        # end if

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

        if self._updateGeomInfo:
            self.mesh.warpMesh()
            newGrid = self.mesh.getSolverGrid()

            if newGrid is not None:
                self.sumb.setgrid(newGrid)
            # end if
                
            self.sumb.updatecoordinatesalllevels()
            self.sumb.updatewalldistancealllevels()
            self.sumb.updateslidingalllevels()
            self.sumb.updatemetricsalllevels()
            self.sumb.updategridvelocitiesalllevels()
            self._updateGeomInfo = False
        # end if

        return 
        
    def getMonitoringVariables(self):
        """Return a list of the text strings describing the variables being
        monitored.

        """
        return self.monnames.keys()
    
    def getConvergenceHistory(self, name):
        """Return an array of the convergence history for a particular
        quantity.

        Keyword arguments:

        name -- the text string for a particular quantity

        """
        try:
            index = self.monnames[name]
        except KeyError:
            mpiPrint('Error: No such quantity %s'%name, comm=self.comm)
            return None
        if (self.myid == 0):
            if (self.sumb.monitor.niterold == 0 and
                self.sumb.monitor.nitercur == 0 and
                self.sumb.iteration.itertot == 0):
                history = None
            elif (self.sumb.monitor.nitercur == 0 and
                  self.sumb.iteration.itertot == 0):
                niterold = self.sumb.monitor.niterold[0]	    
                history = self.sumb.monitor.convarray[:niterold+1, index]
            else:
                history = self.sumb.monitor.convarray[:, index]
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
        finalRes = self.sumb.adjointpetsc.adjreshist[finalIt-1]
        fail = self.sumb.killsignals.adjointfailed

        return startRes, finalRes, fail

    def getResNorms(self):
        '''Return the initial, starting and final Res Norms'''
        return \
            numpy.real(self.sumb.nksolvervars.totalr0), \
            numpy.real(self.sumb.nksolvervars.totalrstart),\
            numpy.real(self.sumb.nksolvervars.totalrfinal)

    def setResNorms(self, initNorm=None, startNorm=None, finalNorm=None):
        ''' Set one of these norms if not none'''
        if initNorm is not None:
            self.sumb.nksolvervars.totalr0 = initNorm
        if startNorm is not None:
            self.sumb.nksolvervars.totalrstart = startNorm
        if finalNorm is not None:
            self.sumb.nksolvervars.finalNorm = finalNorm

        return 
    
    def getMeshIndices(self):
        ndof_1_instance = self.sumb.adjointvars.nnodeslocal[0]*3
        indices = self.sumb.getcgnsmeshindices(ndof_1_instance)

        return indices
    
    def getdRdXvPsi(self, groupName=None, objective=None, flowCase=None):

        self.setFlowCase(flowCase, True, 'getdRdXvpsi') 

        # Get objective
        obj, aeroObj = self._getObjective(objective)

        if obj in self.flowCases[self.curFlowCase]['adjoints'].keys():
            psi = self.flowCases[self.curFlowCase]['adjoints'][obj]
        else:
            mpiPrint('%s adjoint for flowCase %s is not computed.'%(
                    obj, self.curFlowCase), comm=self.comm)
            sys.exit(1)
        # end if

        # Now call getdrdxvpsi WITH the psi vector:
        dxvSolver = self.sumb.getdrdxvpsi(self.getSpatialSize(), psi)

        # If we are doing a prescribed motion TS motion, we need to
        # convert this back to a single instance 
        if self._prescribedTSMotion():
            ndof_1_instance = self.sumb.adjointvars.nnodeslocal[0]*3
            dxvSolver = self.sumb.spectralprecscribedmotion(
                dxvSolver, ndof_1_instance)
        # end if

        if groupName is not None:
            self.mesh.warpDeriv(dxvSolver)
            dxs = self.mesh.getdXs(groupName)
            return dxs
        else:
            return dxvSolver
        # end if

    def _prescribedTSMotion(self):
    
        if self.getOption('alphamode') or self.getOption('betamode') or \
                self.getOption('machmode') or self.getOption('pmode') or \
                self.getOption('qmode') or self.getOption('rmode') or \
                self.getOption('altitudemode'):

            return True
        else:
            return False
        # end if

    def getdRdXvVec(self, inVec, groupName, flowCase=None):
        self.setFlowCase(flowCase, True, 'getdRdXvVec')

        ndof = self.sumb.adjointvars.nnodeslocal[0]*3

        # Now call getdrdxvpsi WITH the psi vector:
        dxvSolver = self.sumb.getdrdxvpsi(self.getSpatialSize(), inVec)
        self.mesh.warpDeriv(dxvSolver)
        dxs = self.mesh.getdXs(groupName)

        return dxs

    def getdRdaPsi(self,  psi):

        if self.nDVAero > 0:
            dIda = self.sumb.getdrdapsi(self.nDVAero, psi)
        else:
            dIda = numpy.zeros((0))
        # end if

        return dIda

    def getdRdwTVec(self, inVec, outVec, flowCase=None):
        ''' Compute the result: outVec = dRdw^T * inVec'''
        self.setFlowCase(flowCase, True, 'getdRdwTVec')

        outVec = self.sumb.getdrdwtvec(inVec, outVec)
        
        return outVec

    def getdFdxVec(self, groupName, vec, flowCase=None):
        # Calculate dFdx * vec and return the result
        self.setFlowCase(flowCase, True, 'getdFdxVec')

        vec = self.mesh.expandVectorByFamily(groupName, vec)
        if len(vec) > 0:
            vec = self.sumb.getdfdxvec(numpy.ravel(vec))
        vec = self.mesh.sectionVectorByFamily(groupName, vec)

        return vec

    def getdFdxTVec(self, groupName, vec, flowCase=None):
        # Calculate dFdx^T * vec and return the result
        self.setFlowCase(flowCase, True, 'getdFdxTVec')

        vec = self.mesh.expandVectorByFamily(groupName, vec)
        if len(vec) > 0:
            vec = self.sumb.getdfdxtvec(numpy.ravel(vec))
        vec = self.mesh.sectionVectorByFamily(groupName, vec)

        return vec

    def computeObjPartials(self, objective, forcePoints=None, flowCase=None):
        self.setFlowCase(flowCase)
        obj, aeroObj = self._getObjective(objective)

        if aeroObj:
            objNum = self.SUmbCostfunctions[obj]

            if self.getOption('useReverseModeAD'):

                # Note: Computeobjective partials MUST be called with the full
                # force pt list.
                if forcePoints is None:
                    [npts, ncells] = self.sumb.getforcesize()
                    nTS  = self.sumb.inputtimespectral.ntimeintervalsspectral
                    forcePoints = numpy.zeros((nTS, npts, 3),self.dtype)
                    for i in xrange(nTS):
                        forcePoints[i] = self.getForcePoints(TS=i)
                    # end force
                # end if

                self.sumb.computeobjpartials(
                    objNum, forcePoints.T, True, True)
            else:
                self.sumb.computeobjectivepartialsfwd(objNum)
            # end if

            # Store the current RHS
                self.flowCases[self.curFlowCase]['adjointRHS'] = obj
        else:
            self.sumb.zeroobjpartials(True, True)
        # end if

        return 

    def getdIdx(self, objective, forcePoints=None, TS=0, groupName=None,
                flowCase=None):
        self.setFlowCase(flowCase, True, 'getdIdx')
        obj, aeroObj = self._getObjective(objective)

        # Compute the partials
        self.computeObjPartials(objective, forcePoints)
        dXv = numpy.zeros(self.getSpatialSize())
        self.sumb.getdidx(dXv)

        # If we are doing a prescribed motion TS motion, we need to
        # convert this back to a single instance 
        if self._prescribedTSMotion():
            ndof_1_instance = self.sumb.adjointvars.nnodeslocal[0]*3
            dXv = self.sumb.spectralprecscribedmotion(dXv, ndof_1_instance)
        # end if

        if groupName is not None:
            # We have a decision to make here: If we have euler
            # analysis, we can do a "surfOnly" meshDerivative since
            # there is no information on the interior anyway. However,
            # if we have a viscous analysis, then we DO have to do a
            # proper mesh warp, its fairly costly, but worth it.

            if self.getOption('equationType') == 'euler':
               self.mesh.warpDeriv(dXv, surfOnly=True)
            else:
                self.mesh.warpDeriv(dXv, surfOnly=False)
            # end if

            dxs = self.mesh.getdXs(groupName)
            return dxs
        else:
            return dXv
        # end if
        
    def getdIda(self, objective, forcePoints=None, flowCase=None):

        self.setFlowCase(flowCase, True, 'getdIda')

        obj, aeroObj = self._getObjective(objective)

        if self.nDVAero > 0:
            
            self.computeObjPartials(objective, forcePoints)
            if aeroObj:
                dIdaLocal = self.sumb.adjointvars.dida
            else:
                dIdaLocal = numpy.zeros_like(self.sumb.adjointvars.dida)
            # end if

            # We must MPI all reuduce
            dIda = self.comm.allreduce(dIdaLocal,  op=MPI.SUM)
        else:
            dIda = numpy.zeros((0))
        # end if

        return dIda
        
    def getdIdw(self, dIdw, objective, forcePoints=None, flowCase=None):
        self.setFlowCase(flowCase, True, 'getdIdw')
        obj, aeroObj = self._getObjective(objective)

        if aeroObj:
            self.computeObjPartials(objective, forcePoints)
            dIdw = self.sumb.getdidw(dIdw)
        # end if

        return dIdw

    def finalizeAdjoint(self):
        '''
        destroy the PESTcKSP context
        '''
        self.releaseAdjointMemory()
        self.sumb.releasememadjoint()
        
        return

    def getStateSize(self):
        '''Return the number of degrees of freedom (states) that are
        on this processor'''

        nstate = self.sumb.flowvarrefstate.nw
        ncells = self.sumb.adjointvars.ncellslocal[0]
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return nstate*ncells*ntime

    def getAdjointStateSize(self):
        '''Return the number of degrees of freedom (states) that are
        on this processor'''
        if self.getOption('frozenTurbulence'):
            nstate = self.sumb.flowvarrefstate.nwf
        else:
            nstate = self.sumb.flowvarrefstate.nw
        #end if

        ncells = self.sumb.adjointvars.ncellslocal[0]
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return nstate*ncells*ntime

    def getSpatialSize(self):
        '''Return the number of degrees of spatial degrees of freedom on this processor.'''

        nnodes = self.sumb.adjointvars.nnodeslocal[0]
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return 3*nnodes*ntime

    def getStates(self):
        '''Return the states on this processor. Used in aerostructural
        analysis'''

        states = self.sumb.getstates(self.getStateSize())

        return states

    def setStates(self, states):
        ''' Set the states on this processor. Used in aerostructural
        analysis'''

        self.sumb.setstates(states)

        return 

    def setAdjoint(self, adjoint, objective=None, flowCase=None):
        '''Sets the adjoint vector externally. Used in coupled solver'''
        self.setFlowCase(flowCase)

        self.sumb.setadjoint(adjoint)
        if objective is not None:
            obj, aeroObj = self._getObjective(objective)
            self.flowCases[self.curFlowCase]['adjoints'][obj] = adjoint.copy()
        # end if

        return

    def getResidual(self, res=None, flowCase=None):
        '''Return the residual on this processor. Used in aerostructural
        analysis'''
        self.setFlowCase(flowCase)
        if res is None:
            res = numpy.zeros(self.getStateSize())
        res = self.sumb.getres(res)
        
        return res

    def computedSdwTVec(self, inVec, outVec, groupName, flowCase=None):
        '''This function computes: outVec = outVec + dFdw^T*inVec'''
        self.setFlowCase(flowCase, True, 'computedSdwTVec')
        phi = self.mesh.expandVectorByFamily(groupName, inVec)
        outVec = self.sumb.getdfdwtvec(numpy.ravel(phi), outVec)

        return outVec

    def getSolution(self, sps=1, flowCase=None):
        ''' Retrieve the solution variables from the solver. Note this
        is a collective function and must be called on all processors
        '''
        self.setFlowCase(flowCase, True, 'getSolution')
        # We should return the list of results that is the same as the
        # possibleObjectives list
        self.sumb.getsolution(sps)

        funcVals = self.sumb.costfunctions.functionvalue
        SUmbsolution = {
            'lift':funcVals[self.sumb.costfunctions.costfunclift-1],
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
            'clq'        :funcVals[self.sumb.costfunctions.costfuncclq-1],
            'cbend'      :funcVals[self.sumb.costfunctions.costfuncbendingcoef-1]
            }
        
        return SUmbsolution

    def computeArea(self, axis, groupName=None, TS=0, flowCase=None):
        """
        Compute the projected area of the surface mesh

        Input Arguments:
           axis, numpy array, size(3): The projection vector
               along which to determine the shadow area
           groupName, str: The group from which to obtain the coordinates.
               This name must have been obtained from addFamilyGroup() or 
               be the default 'all' which contains all surface coordiantes 

        Output Arguments: 
            Area: The resulting area  
            """
        self.setFlowCase(flowCase)

        cfdForcePts = self.getForcePoints(TS)
        if len(cfdForcepts) > 0:
            areas = self.sumb.getareas(cfdForcePts.T, TS+1, axis).T
        else:
            areas = numpy.zeros((0,3), self.dtype)
        # end if

        if groupName is not None:
            areas = self.mesh.sectionVectorByFamily(groupName, areas)
        # end if

        # Now we do an mpiallreduce with sum:
        area = self.comm.allreduce(numpy.sum(areas), op=MPI.SUM)
        
        return area

    def computeAreaSensitivity(self, axis, groupName=None, TS=0, flowCase=None):
        """ 
        Compute the projected area of the surface mesh

        Input Arguments:
           axis, numpy array, size(3): The projection vector
               along which to determine the shadow area  
           groupName, str: The group from which to obtain the coordinates.
               This name must have been obtained from addFamilyGroup() or 
               be the default 'all' which contains all surface coordiantes 

        Output Arguments:
            Area: The resulting area    
            """
        self.setFlowCase(flowCase)
        cfdForcePts = self.getForcePoints(TS)

        if len(cfdForcePts) > 0:
            da = self.sumb.getareasensitivity(cfdForcePts.T, TS+1, axis).T
        else:
            da = numpy.zeros((0,3), self.dtype)
        # end if

        if groupName is not None:
            da = self.mesh.sectionVectorByFamily(groupName, da)
        # end if

        return da

    def _getObjective(self, objective):
        '''Check to see if objective is one of the possible
        objective. If it is, return the obj value for SUmb and
        True. Otherwise simply return the objective string and
        False'''
        if objective in self.possibleObjectives.keys():
            obj = self.possibleObjectives[objective]
            aeroObj = True
        else:
            obj = objective 
            aeroObj = False
        # end try

        return obj, aeroObj

    def setOption(self, name, value):
        '''
        Set Solver Option Value 
        '''
        name = name.lower()

        # Check to see if we have a deprecated option. Print a useful
        # warning that this is deprecated.
        if name in self.deprecatedOptions.keys():
            mpiPrint('+'+'-'*78+'+',comm=self.comm)
            mpiPrint('| WARNING: Option: \'%-29s\' is a deprecated SUmb Option |'%name,comm=self.comm)
            mpiPrint('| %-78s'%self.deprecatedOptions[name],comm=self.comm)
            mpiPrint('+'+'-'*78+'+',comm=self.comm)
            return
        # end if

        # Try to the option in the option dictionary
        defOptions = self.options['defaults']
        try: 
            defOptions[name]
        except: 
            mpiPrint('+'+'-'*78+'+',comm=self.comm)
            mpiPrint('| WARNING: Option: \'%-30s\' is not a valid SUmb Option |'%name,comm=self.comm)
            mpiPrint('+'+'-'*78+'+',comm=self.comm)
            return
        # end try
        
        # Now we know the option exists, lets check if the type is ok:
        if type(value) == self.options[name][0]:
            # Just set:
            self.options[name] = [type(value),value]
        else:
            mpiPrint('+'+'-'*78+'+',comm=self.comm)
            mpiPrint('| ERROR: Datatype for Option %-35s was not valid |'%name,comm=self.comm)
            mpiPrint('|        Expected data type is %-47s |'%self.options[name][0],comm=self.comm)
            mpiPrint('|        Received data type is %-47s |'%type(value),comm=self.comm)
            mpiPrint('+'+'-'*78+'+',comm=self.comm)
            sys.exit(1)
        # end if
            
        # If the option is in the ignoredOption list, we just return. 
        if name in self.ignoreOptions:
            return
        
        # Do special Options individually
        if name in self.specialOptions:
            if name in ['monitorvariables',
                        'surfacevariables',
                        'volumevariables',
                        'isovariables']:
                varStr = ''
                for i in xrange(len(value)):
                    varStr = varStr + value[i] + '_'
                # end if
                varStr = varStr[0:-1] # Get rid of last '_'
                if name == 'monitorvariables':
                    self.sumb.monitorvariables(varStr)
                if name == 'surfacevariables':
                    self.sumb.surfacevariables(varStr)
                if name == 'volumevariables':
                    self.sumb.volumevariables(varStr)
                if name == 'isovariables':
                    self.sumb.isovariables(varStr)

            if name == 'isosurface':
                # We have a bit of work to do...extract out the
                # names, and there can be more than 1 value per variables
                var = []
                val = []
                isoDict = value
                for key in isoDict.keys():
                    
                    isoVals = numpy.atleast_1d(isoDict[key])
                    for i in xrange(len(isoVals)):
                        var.append(key)
                        val.append(isoVals[i])
                    # end for
                # end for
                val = numpy.array(val)

                self.sumb.initializeisosurfacevariables(val)
                for i in xrange(len(val)):
                    self.sumb.setisosurfacevariable(var[i], i+1)
                # end for

            # end if
            if name == 'metricconversion':
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
            #Convert the value to lower case:
            value = self.optionMap[name][value.lower()]
        # end if

        # If value is a string, put quotes around it and make it
        # the correct length, otherwise convert to string
        if isinstance(value, str): 
            spacesToAdd = self.optionMap[name]['len'] - len(value)
            value = '\'' + value + ' '*spacesToAdd + '\''
        else:
            value = str(value)
        # end if
      
        # Exec str is what is actually executed:
        execStr = 'self.sumb.'+self.optionMap[name]['location'] + '=' + value

        exec(execStr)

        return
        
    def getOption(self, name):

        # Redefine the getOption def from the base class so we can
        # mane sure the name is lowercase

        defOptions = self.options['defaults']
        if defOptions.has_key(name.lower()):
            return self.options[name.lower()][1]
        else:    
            raise InputError(repr(name) + ' is not a valid option name')
        # end if
        
        return

class SUmbDummyMesh(object):
    """
    Represents a dummy Multiblock structured Mesh for SUmb
    """
 
    def __init__(self):
        return
        
    def getSurfaceCoordinates(self, groupName):
    
        return 

    def getSurfaceConnectivity(self, groupName):

        return

    def setSurfaceCoordinates(self, groupName, coordinates):
       
        return 

    def addFamilyGroup(self, groupName, families=None):

        return 
   
# =========================================================================
#                         Interface Functionality
# =========================================================================

    def setExternalMeshIndices(self, ind):

        return 

    def setExternalSurface(self, patchNames, patchSizes, conn, pts):

        return

    def getSolverGrid(self):
      
        return

    def getdXs(self, groupName):

        return 

# ==========================================================================
#                        Output Functionality
# ==========================================================================

    def writeVolumeGrid(self, fileName):
        '''write volume mesh'''

        return

    def writeSurfaceGrid(self, fileName):
        '''write surface mesh'''

        return

    def writeFEGrid(self, fileName):
        ''' write the FE grid to file '''
        return

# =========================================================================
#                         Utiliy Functions
# =========================================================================

    def getMeshQuality(self, bins=None):
        
        return

# =========================================================================
#                      Mesh Warping Functionality
# =========================================================================

    def warpMesh(self):
        ''' 
        Run either the solid warping scheme or the algebraic
        warping scheme depending on the options
        '''

    def WarpDeriv(self, solverdXdv):

        return

#==============================================================================
# SUmb Analysis Test
#==============================================================================
if __name__ == '__main__':
    
    # Test SUmb
    mpiPrint('Testing ...')
    sumb = SUMB()
    print(sumb)

