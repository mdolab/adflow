#!/usr/bin/python
from __future__ import print_function
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
import sys, string
import time
import copy

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from mdo_import_helper import MPI, import_modules, mpiPrint, MExt
exec(import_modules('pyAero_solver'))

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
        
        Documentation last updated:  July. 03, 2008 - C.A.(Sandy) Mader
        '''
        
        name = 'SUMB'
        category = 'Three Dimensional CFD'
        def_opts = {
            # Common Paramters
            'gridfile':[str, 'default.cgns'],
            'restartfile':[str, 'default_restart.cgns'],
            'probname':[str, ''],
            'outputdir':[str, './'],
            'solrestart':[bool, False],
            'writesolution':[bool, True],
            'writemesh':[bool, False],
            'storerindlayer':[bool, True],

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
            'blocksplitting':[bool, False],
            'loadimbalance':[float, 0.1],
       
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
            'finitedifferencepc':[bool, True],
            'usereversemodead':[bool, True],
            'lowmemory':[bool, True],
            'applyadjointpcsubspacesize':[int, 20]
            }

        informs = {
            }

        # Load the compiled module using MExt, which allow multiple imports
        try: 
            self.sumb
        except:
            self.sumb = MExt('sumb')._module
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
            'finitedifferencepc':{'location':'inputadjoint.finitedifferencepc'},
            }                

        # These "ignore_options" are NOT actually, ignored, rather,
        # they DO NOT GET SET IN THE FORTRAN CODE. Rather, they are
        # used strictly in Python

        self.ignore_options = [
            'defaults',
            'storehistory',
            'numbersolutions',
            'writesolution',
            'writemesh',
            'familyrot',  # -> Not sure how to do
            'lowmemory',
            'autosolveretry',
            'autoadjointretry',
            'usereversemodead'
            ]

        self.special_options = ['surfacevariables',
                                'volumevariables',
                                'monitorvariables',
                                'metricconversion',
                                'outputdir',
                                'probname']

        self.storedADjoints = {}

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
            self, name, category, def_opts, informs, options=options)
        self.sumb.inputio.autoparameterupdate = False

        # Set the external Mesh Warping is provided
        self._update_geom_info = False
        if mesh is None:
            self.mesh = SUmbDummyMesh()
        else:
            self.mesh = mesh
        # end if

        # Set Flags that are used to keep of track of what is "done"
        # in fortran
        self.allInitialized = False    # All flow solver initialization   
        self.adjointPreprocessed = False

        # Matrix Setup Flags
        self.adjointSetup = False
        self.adjointRHS = None # When this is setup, it has
                               # the current objective
        
        self._update_geom_info = False
        self._update_period_info = True
        self._update_vel_info = True
        self.fatalFail = False
        self.solve_failed = False
        self.adjoint_failed = False
        self.dtype = 'd'

        # Write the intro message
        self.sumb.writeintromessage()

        return

    def initialize(self, aero_problem, partitionOnly=False):
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
                self.setOption(option.lower(), self.options[option][1])
            # end if
        # end for
      
        # Do the remainder of the operations that would have been done
        # had we read in a param file
        self.sumb.iteration.deforming_grid = True

        self.setMachNumber(aero_problem)
        self.setRefState(aero_problem)
        self.setPeriodicParams(aero_problem)

        self.sumb.dummyreadparamfile()

        #This is just to flip the -1 to 1 possibly a memory issue?
        self.sumb.inputio.storeconvinneriter = \
            abs(self.sumb.inputio.storeconvinneriter)

    
        mpiPrint(' -> Partitioning and Reading Grid', comm=self.comm)
        self.sumb.partitionandreadgrid()
        if partitionOnly:
            return
        # end if

        mpiPrint(' -> Preprocessing', comm=self.comm)
        self.sumb.preprocessing()

        mpiPrint(' -> Initializing flow', comm=self.comm)
        self.sumb.initflow()
        self.setInflowAngle(aero_problem)

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
            self.sumb.getpatchname(i, tmp)
            patchnames.append(
                ''.join([tmp[j] for j in range(256)]).lower().strip())
            patchsizes.append(self.sumb.getpatchsize(i))
        # end for
 
        conn = self.getForceConnectivity()
        pts  = self.getForcePoints()

        self.mesh.setExternalSurface(patchnames, patchsizes, conn, pts)
        
        # Solver is initialize
        self.allInitialized = True
        self.initAdjoint()
        
        return

    def setInflowAngle(self, aero_problem):
        '''
        Set the alpha and beta fromthe desiggn variables
        '''
        
        [velDir, liftDir, dragDir] = self.sumb.adjustinflowangleadjts(\
            (aero_problem._flows.alpha*(numpy.pi/180.0)),
            (aero_problem._flows.beta*(numpy.pi/180.0)),
            aero_problem._flows.liftIndex)
        self.sumb.inputphysics.veldirfreestream = velDir
        self.sumb.inputphysics.liftdirection = liftDir
        self.sumb.inputphysics.dragdirection = dragDir

        if self.sumb.inputiteration.printiterations:
            mpiPrint('-> Alpha... %f %f'%(
                    numpy.real(aero_problem._flows.alpha*(numpy.pi/180.0)),
                    numpy.real(aero_problem._flows.alpha)), comm=self.comm)

        #update the flow vars
        self.sumb.updateflow()
        self._update_vel_info = True

        return

    def setElasticCenter(self, aero_problem):
        '''
        set the value of pointRefEC for the bending moment calculation
        '''
        
        self.sumb.inputphysics.pointrefec[0] = aero_problem._geometry.xRootec\
            *self.metricConversion
        self.sumb.inputphysics.pointrefec[1] = aero_problem._geometry.yRootec\
            *self.metricConversion
        self.sumb.inputphysics.pointrefec[2] = aero_problem._geometry.zRootec\
            *self.metricConversion
    
    def setReferencePoint(self, aero_problem):
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

    def setRotationRate(self, aero_problem):
        '''
        Set the rotational rate for the grid
        '''
        a  = numpy.sqrt(self.sumb.flowvarrefstate.gammainf*\
                      self.sumb.flowvarrefstate.pinfdim/ \
                      self.sumb.flowvarrefstate.rhoinfdim)
        V = (self.sumb.inputphysics.machgrid+self.sumb.inputphysics.mach)*a
        
        p = aero_problem._flows.phat*V/aero_problem._refs.bref
        q = aero_problem._flows.qhat*2*V/aero_problem._refs.cref
        r = aero_problem._flows.rhat*V/aero_problem._refs.bref

        self.sumb.updaterotationrate(p, r, q)
        self._update_vel_info = True

        return
    
    def setRefArea(self, aero_problem):
        self.sumb.inputphysics.surfaceref = aero_problem._refs.sref*self.metricConversion**2
        self.sumb.inputphysics.lengthref = aero_problem._refs.cref*self.metricConversion
        
        return

    def setPeriodicParams(self, aero_problem):
        '''
        Set the frequecy and amplitude of the oscillations
        '''
        if  self.getOption('alphaMode'):
            self.sumb.inputmotion.degreepolalpha = int(aero_problem._flows.degreePol)
            self.sumb.inputmotion.coefpolalpha = aero_problem._flows.coefPol
            self.sumb.inputmotion.omegafouralpha   = aero_problem._flows.omegaFourier
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
        # end if

        self._update_period_info = True
        self._update_geom_info = True
        self._update_vel_info = True
 
        return

    def setMachNumber(self, aero_problem):
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
            self.sumb.inputphysics.machcoef = aero_problem._flows.mach
            self.sumb.inputphysics.machgrid = aero_problem._flows.mach
            self.sumb.inputmotion.gridmotionspecified = True
        else:
            self.sumb.inputphysics.mach = aero_problem._flows.mach 
            self.sumb.inputphysics.machcoef = aero_problem._flows.mach
            self.sumb.inputphysics.machgrid = 0.0
        # end if

        return

    def setRefState(self, aero_problem):
        ''' Set the Pressure, density and viscosity/reynolds number
        from the aero_problem
        '''
        self.sumb.flowvarrefstate.pref = aero_problem._flows.P
        self.sumb.flowvarrefstate.rhoref = aero_problem._flows.rho
        self.sumb.flowvarrefstate.tref = aero_problem._flows.T

        # Reynolds number info not setup yet...

        return

    def resetAdjoint(self, obj):
        '''
        Reset a possible stored adjoint 'obj'
        '''

        if obj in self.storedADjoints.keys():
            self.storedADjoints[obj][:] = 0.0
        # end if

        return

    def resetFlow(self, aeroProblem=None):
        '''
        Reset the flow for the complex derivative calculation
        '''

        if aeroProblem is not None:
            self.setInflowAngle(aeroProblem)
            self.setMachNumber(aeroProblem)
            self.setRefState(aeroProblem)
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
        self.sumb.inputiteration.currentlevle = strLvl
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

    def __solve__(self, aero_problem, nIterations=500, MDCallBack=None):
        
        '''
        Run Analyzer (Analyzer Specific Routine)
        
        Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
        '''
        # Release adjoint memory in case an adjoint was previously solved
        self.releaseAdjointMemory()
        self.adjointRHS         = None
        self.callCounter += 1

        # Run Initialize, if already run it just returns.
        self.initialize(aero_problem)

        #set inflow angle, refpoint etc.
        self.setMachNumber(aero_problem)
        self.setPeriodicParams(aero_problem)
        self.setInflowAngle(aero_problem)
        self.setReferencePoint(aero_problem)
        #self.setElasticCenter(aero_problem)
        self.setRotationRate(aero_problem)
        self.setRefArea(aero_problem)
        self.setRefState(aero_problem)

        # Run Solver
        t0 = time.time()

        # set the number of cycles for this call
        self.sumb.inputiteration.ncycles = nIterations

        # Cold Start:
        if self.sumb.monitor.niterold == 0 and \
            self.sumb.monitor.nitercur == 0 and \
            self.sumb.iteration.itertot == 0:
            if self.myid == 0:
                desired_size = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles
                self.sumb.allocconvarrays(desired_size)
            # end if
        else:
            # More Time Steps / Iterations OR a restart
            # Reallocate convergence history array and time array
            # with new size,  storing old values from previous runs
            if self.getOption('storeHistory'):
                current_size = len(self.sumb.monitor.convarray)
                desired_size = current_size + self.sumb.inputiteration.ncycles+1
                self.sumb.monitor.niterold  = self.sumb.monitor.nitercur+1
            else:
                self.sumb.monitor.nitercur  = 0
                self.sumb.monitor.niterold  = 1
                desired_size = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles +1
            # end if
            # Allocate Arrays
            if self.myid == 0:
                self.sumb.allocconvarrays(desired_size)

            self.sumb.inputiteration.mgstartlevel = 1
            self.sumb.iteration.itertot = 0
        # end if

        if self.getOption('equationMode') == 'unsteady':
            self.sumb.alloctimearrays(self.getOption('nTimeStepsFine'))

        # Reset Fail Flags
        self.sumb.killsignals.routinefailed =  False
        self.sumb.killsignals.fatalfail = False
        self.solve_failed =  self.fatalFail = False

        self._updatePeriodInfo()
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
            self.solve_failed = True
            return
     
        t1 = time.time()
        # Call the Solver
        if MDCallBack is None:
            self.sumb.solver()
        else:
            self.sumb.solverunsteadymd(MDCallBack)
        # end if

        # Assign Fail Flags
        self.solve_failed = self.sumb.killsignals.routinefailed
        self.fatalFail = self.sumb.killsignals.fatalfail

        # Reset Flow if there's a fatal fail reset and return;
        # --> Do not write solution
        if self.fatalFail:
            self.resetFlow()
            return
        # end if

        t2 = time.time()
        sol_time = t2 - t0
        self.update_time = t1 - t0

        if self.getOption('printTiming'):
            mpiPrint('Solution Time: %10.3f sec'%sol_time, comm=self.comm)
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
        # end if
        
        return

    def solveCL(self, aeroProblem, CL_star, nIterations=500, alpha0=0, 
                delta=0.5, tol=1e-3):
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
        self.__solve__(aeroProblem, nIterations=nIterations)
        sol = self.getSolution()
        fnm2 =  sol['cl'] - CL_star

        for iIter in xrange(20):
            # We need to reset the flow since changing the alpha leads
            # to problems with the NK solver
            self.resetFlow()

            # Set current alpha
            aeroProblem._flows.alpha = anm1

            # Solve for n-1 value (anm1)
            self.__solve__(aeroProblem, nIterations=nIterations)
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

    def getSurfaceCoordinates(self, group_name):
        ''' 
        See MultiBlockMesh.py for more info
        '''
        return self.mesh.getSurfaceCoordinates(group_name)

    def setSurfaceCoordinates(self, group_name, coordinates):
        ''' 
        See MultiBlockMesh.py for more info
        '''
        self._update_geom_info = True
        self.mesh.setSurfaceCoordinates(group_name, coordinates)

        return 

    def getSurfaceConnectivity(self, group_name):
        '''
        See MultiBlockMesh.py for more info
        '''
        return self.mesh.getSurfaceConnectivity(group_name)
        
    def writeMeshFile(self, filename=None):
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

    def writeVolumeSolutionFile(self, filename=None, writeGrid=True):
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

    def writeSurfaceSolutionFile(self, *filename):
        '''Write the current state of the surface flow solution to a CGNS file.
        Keyword arguments:
        filename -- the name of the file (optional)
        '''
        if (filename):
            self.sumb.inputio.surfacesolfile[:] = ''
            self.sumb.inputio.surfacesolfile[0:len(filename[0])] = filename[0]
        # end if
        self.sumb.monitor.writegrid=False
        self.sumb.monitor.writevolume=False
        self.sumb.monitor.writesurface=True
        self.sumb.writesol()

        return

    def writeForceFile(self, file_name, TS=0, group_name='all', cfd_force_pts=None):
        '''This function collects all the forces and locations and
        writes them to a file with each line having: X Y Z Fx Fy Fz.
        This can then be used to set a set of structural loads in TACS
        for structural only optimization 

        Like the getForces() routine, an external set of forces may be
        passed in on which to evaluate the forces. This is only
        typically used in an aerostructural case. 

        '''
      
        if self.mesh is None:
            mpiPrint('Error: A pyWarp mesh be specified to use writeForceFile',
                     comm=self.comm)
            return

        # Now we need to gather the data:
        if cfd_force_pts is None:
            pts = self.comm.gather(self.getForcePoints(TS), root=0)
        else:
            pts = self.comm.gather(cfd_force_pts)
        # end if
            
        # Forces are still evaluated on the displaced surface so do NOT pass in pts.
        forces = self.comm.gather(self.getForces(group_name, TS=TS), root=0)
        conn   = self.comm.gather(self.mesh.getSurfaceConnectivity(group_name),
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
            f = open(file_name, 'w')
            
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

    def getForces(self, group_name=None, cfd_force_pts=None, TS=0):
        ''' Return the forces on this processor. Use
        cfd_force_pts to compute the forces if given
        '''

        if cfd_force_pts is None:
            cfd_force_pts = self.getForcePoints(TS)
        # end if
        
        if len(cfd_force_pts) > 0:
            forces = self.sumb.getforces(cfd_force_pts.T, TS).T
        else:
            forces = numpy.zeros((0,3),self.dtype)
        # end if

        if group_name is not None:
            forces = self.mesh.sectionVectorByFamily(group_name, forces)

        return forces

    def getForcePoints(self, TS=0):
        [npts, ncell, nTS] = self.sumb.getforcesize()
        pts = numpy.zeros((nTS, npts, 3),self.dtype)
        self.sumb.getforcepoints(pts.T)
        
        return pts[TS]

    def getForceConnectivity(self):
        conn_size = self.sumb.getforceconnectivitysize()
        conn =  numpy.zeros((conn_size, 4), dtype='intc')
        self.sumb.getforceconnectivity(numpy.ravel(conn))

        return conn

    def verifyForces(self, cfd_force_pts=None):

        # Adjoint must be initialized for force verification

        self.initAdjoint()
        if cfd_force_pts is None:
            cfd_force_pts = self.getForcePoints()
        # end if

        self.sumb.verifyforces(cfd_force_pts.T)

        return

    def verifyBendingPartial(self):
        self.sumb.verifybendingderivatives()
        
        return

    def globalNKPreCon(self, in_vec, out_vec):
        '''This function is ONLY used as a preconditioner to the
        global Aero-Structural system'''

        out_vec = self.sumb.applypc(in_vec, out_vec)
        
        return out_vec

    def globalAdjointPreCon(self, in_vec, out_vec):
        ''' This function is ONLY used as a preconditioner for the
        global Aero-Structural Adjoint system'''

        out_vec = self.sumb.applyadjointpc(in_vec, out_vec)

        return out_vec

    def verifydCdx(self, objective):
        '''
        call the reouttine to compare the partial dIda
        against FD
        '''

        self.setupAdjoint()

        # Short form of objective--easier code reading
        obj = self.possibleObjectives[objective.lower()]
        costFunc =  self.SUmbCostfunctions[obj]
        self.sumb.verifydcfdx(1, costFunc)
        
        return

    def verifydCdw(self, objective):
	
	self.sumb.verifydcdwfile(1)

	return

    def verifydIdw(self, objective):
        '''
        run compute obj partials, then write to a file...
        '''
        self.setupAdjoint()
        self.computeObjPartials(objective)
        obj = self.possibleObjectives[objective.lower()]
        filename= self.getOption('outputDir') + '/' +'ADw%s'%(obj)
        costFunc =  self.SUmbCostfunctions[obj]
        level = 1

        self.sumb.verifydidwfile(level, costFunc, filename)
        
        return

    def verifydIdx(self, objective):
        '''
        run compute obj partials, then write to a file...
        '''
        self.setupAdjoint()
        self.computeObjPartials(objective)
        obj = self.possibleObjectives[objective.lower()]
        filename= self.getOption('outputDir') + '/' +'ADx%s'%(obj)

        costFunc =  self.SUmbCostfunctions[obj]
        level = 1

        self.sumb.verifydidxfile(level, costFunc, filename)

        return

    def verifydRdw(self):
        ''' run the verify drdw scripts in fortran'''
        # Make sure adjoint is initialize
        self.initAdjoint()
        self.setupAdjoint()
        # end if
	level = 1
        self.sumb.iteration.currentlevel=level
        self.sumb.iteration.groundlevel=level
	self.sumb.verifydrdwfile(1)

	return

    def verifydRdx(self):
        ''' run the verify drdw scripts in fortran'''
        # Make sure adjoint is initialize
        self.initAdjoint()
	level = 1
        self.sumb.iteration.currentlevel=level
        self.sumb.iteration.groundlevel=level
	self.sumb.verifydrdxfile(1)
        
	return

    def verifydRda(self):
        ''' run the verify drdw scripts in fortran'''
        # Make sure adjoint is initialize
        self.initAdjoint()
        #if not self.adjointMatrixSetup:
        #    #self.sumb.createpetscvars()
        #    self.setupAdjoint(None)#(forcePoints)
        # end if
        level = 1
        self.sumb.iteration.currentlevel=level
        self.sumb.iteration.groundlevel=level
	self.sumb.verifydrdextrafile(1)

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
                exec_str = 'self.sumb.' + self.possibleAeroDVs[self.aeroDVs[i]] + \
                           '= %d'%(i)
                # Leave this zero-based since we only need to use it in petsc
                exec(exec_str)
            # end for
        # end if

        #Set the mesh level and timespectral instance for this
        #computation
                
        if not self.adjointPreprocessed:
            self.sumb.preprocessingadjoint()
            self.adjointPreprocessed = True
        # end if

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

        return

    def setupAdjoint(self, reform=False):
        '''
        Setup the data structures required to solve the adjoint problem
        '''
        # Destroy the NKsolver to free memory -- Call this even if the
        # solver is not used...a safeguard check is done in Fortran
        self.sumb.destroynksolver()

        # Run initAdjoint in case this is the first adjoint solve
        self.initAdjoint()

        # For now, just create all the petsc variables
        if not self.adjointSetup or reform:
            self.sumb.createpetscvars()

            if self.getOption('useReverseModeAD'):
                self.sumb.setupallresidualmatrices()
            else:
                self.sumb.setupallresidualmatricesfwd()
            # end if
                
            # Create coupling matrix struct whether we need it or not
            self.sumb.setupcouplingmatrixstruct(self.getForcePoints().T)
            
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
  
        load_inbalance, face_inbalance = self.sumb.checkpartitioning(nprocs)
                
        return load_inbalance, face_inbalance
    
    def releaseAdjointMemory(self):
        '''
        release the PETSc Memory that have been allocated
        '''
        if self.adjointSetup:
            self.sumb.destroypetscvars()
            self.adjointSetup = False

        return

    def _on_adjoint(self, objective, forcePoints=None, structAdjoint=None, 
                    group_name=None):

        # Try to see if obj is an aerodynamic objective. If it is, we
        # will have a non-zero RHS, otherwise its an objective with a
        # zero aerodynamic RHS

        obj, aeroObj = self._getObjective(objective)

        # Setup adjoint matrices/vector as required
        self.setupAdjoint()

        # Check to see if the RHS Partials have been computed
        if not self.adjointRHS == obj:
            self.computeObjPartials(objective, forcePoints)
        # end if

        # Check to see if we need to agument the RHS with a structural
        # adjoint:
        if structAdjoint is not None and group_name is not None:
            phi = self.mesh.expandVectorByFamily(group_name, structAdjoint)
            self.sumb.agumentrhs(numpy.ravel(phi))
        # end if

        # If we have saved adjoints, 
        if self.getOption('restartAdjoint') or self.getOption('lowMemory'):
            # Objective is already stored, so just set it
            if obj in self.storedADjoints.keys():
                self.sumb.setadjoint(self.storedADjoints[obj])
            else:
                # Objective is not yet run,  allocated zeros and set
                self.storedADjoints[obj]= numpy.zeros(self.getStateSize(),float)
                self.sumb.setadjoint(self.storedADjoints[obj])
            # end if
        # end if

        # Actually Solve the adjoint system
        self.sumb.solveadjointtransposepetsc()

        # Possibly try another solve
        if self.sumb.killsignals.adjointfailed and self.getOption('restartAdjoint'):
            # Only retry if the following conditions are met:

            # 1. restartAdjoint is true -> that is we were starting
            # from a non-zero starting point

            # 2. The stored adjoint must have been already set at
            # least once; that is we've already tried one solve
            
            if  obj in self.storedADjoints.keys():
                self.storedADjoints[obj][:] = 0.0 # Always reset a
                                                  # stored adjoint

                if self.getOption('autoAdjointRetry'):
                    self.sumb.solveadjointtransposepetsc()
                # end if
            # end if
        # end if

        # Now set the flags and possibly reset adjoint
        if self.sumb.killsignals.adjointfailed == False:
            self.adjoint_failed = False
            # Copy out the adjoint to store
            if self.getOption('restartAdjoint') or self.getOption('lowMemory'):
                self.storedADjoints[obj] =  \
                    self.sumb.getadjoint(self.getStateSize())
            # end if
        else:
            self.adjoint_failed = True
            # Reset stored adjoint
            if self.getOption('restartAdjoint') or self.getOption('lowMemory'):
                self.storedADjoints[obj][:] = 0.0
            # end if
        # end if
       
        return

    def totalSurfaceDerivative(self, objective):
        # The adjoint vector is now calculated so perform the
        # following operation to produce dI/dX_surf:
        # (p represents partial, d total)
        # dI/dX_s = pI/pX_s - (dXv/dXs)^T * ( dRdX_v^T * psi)
        # 
        # The derivative wrt the surface captures the effect of ALL
        # GLOBAL Multidisciplinary variables -- any DV that changes
        # the surface. 

        obj, aeroObj = self._getObjective(objective)

        # NOTE: do dRdxvPsi MUST be done first since this
        # allocates spatial memory if required.
        dIdxs_2 = self.getdRdXvPsi('all', objective)
          
        # Direct partial derivative contibution 
        dIdxs_1 = self.getdIdx(objective, group_name='all')

        # Total derivative of the obective with surface coordinates
        dIdXs = dIdxs_1 - dIdxs_2

        return dIdXs

    def totalAeroDerivative(self, objective):
        # The adjoint vector is now calculated. This function as above
        # computes dI/dX_aero = pI/pX_aero - dR/dX_aero^T * psi. The
        # "aero" variables are intrinsic ONLY to the aero
        # discipline. Nothing in the structural process should depend
        # on these functions directly. 

        obj, aeroObj = self._getObjective(objective)

        if self.getOption('lowMemory') or self.getOption('restartAdjoint'):
            if obj in self.storedADjoints.keys():
                psi = self.storedADjoints[obj]
            else:
                mpiPrint('%s adjoint is not computed.'%(obj), comm=self.comm)
                sys.exit(1)
            # end if
        else:
            psi = self.sumb.getadjoint(self.getStateSize())
        # end if

        # Direct partial derivative contibution 
        dIda_1 = self.getdIda(objective)

        # dIda contribution for drda^T * psi
        dIda_2 = self.getdRdaPsi(psi)

        # Total derivative of the obective wrt aero-only DVs
        dIda = dIda_1 - dIda_2

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

        if self._update_geom_info:
            self.mesh.warpMesh()
            #self.mesh.writeVolumeGrid('/scratch/j/jmartins/kenway/warped_grid.cgns')
            #self.mesh.writeSurfaceGrid('/scratch/j/jmartins/kenway/warped_surf.cgns')
            newGrid = self.mesh.getSolverGrid()

            if newGrid is not None:
                self.sumb.setgrid(newGrid)
            # end if
                
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
    
    def getdRdXvPsi(self, group_name=None, objective=None):

        self.setupAdjoint()

        # Get objective
        obj, aeroObj = self._getObjective(objective)

        if self.getOption('lowMemory') or self.getOption('restartAdjoint'):
            if obj in self.storedADjoints.keys():
                psi = self.storedADjoints[obj]
            else:
                mpiPrint('%s adjoint is not computed.'%(obj), comm=self.comm)
                sys.exit(1)
            # end if
        else:
            psi = self.sumb.getadjoint(self.getStateSize())
        # end if

        # Now call getdrdxvpsi WITH the psi vector:
        dxv_solver = self.sumb.getdrdxvpsi(self.getSpatialSize(), psi)

        # If we are doing a prescribed motion TS motion, we need to
        # convert this back to a single instance 
        if self._prescribedTSMotion():
            ndof_1_instance = self.sumb.adjointvars.nnodeslocal[0]*3
            dxv_solver = self.sumb.spectralprecscribedmotion(
                dxv_solver, ndof_1_instance)
        # end if

        if group_name is not None:
            self.mesh.warpDeriv(dxv_solver)
            dxs = self.mesh.getdXs(group_name)
            return dxs
        else:
            return dxv_solver
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

    def getdRdXvVec(self, in_vec, group_name):

        ndof = self.sumb.adjointvars.nnodeslocal[0]*3

        # Now call getdrdxvpsi WITH the psi vector:
        dxv_solver = self.sumb.getdrdxvpsi(ndof, in_vec)
        self.mesh.warpDeriv(dxv_solver)
        dxs = self.mesh.getdXs(group_name)

        return dxs

    def getdRdaPsi(self,  psi):

        if self.nDVAero > 0:
            dIda = self.sumb.getdrdapsi(self.nDVAero, psi)
        else:
            dIda = numpy.zeros((0))
        # end if

        return dIda

    def getdRdwTVec(self, in_vec, out_vec):
        ''' Compute the result: out_vec = dRdw^T * in_vec'''

        out_vec = self.sumb.getdrdwtvec(in_vec, out_vec)
        
        return out_vec

    def getdFdxVec(self, group_name, vec):
        # Calculate dFdx * vec and return the result
        vec = self.mesh.expandVectorByFamily(group_name, vec)
        if len(vec) > 0:
            vec = self.sumb.getdfdxvec(numpy.ravel(vec))
        vec = self.mesh.sectionVectorByFamily(group_name, vec)

        return vec

    def getdFdxTVec(self, group_name, vec):
        # Calculate dFdx^T * vec and return the result
        vec = self.mesh.expandVectorByFamily(group_name, vec)
        if len(vec) > 0:
            vec = self.sumb.getdfdxtvec(numpy.ravel(vec))
        vec = self.mesh.sectionVectorByFamily(group_name, vec)

        return vec

    def computeObjPartials(self, objective, forcePoints=None):

        obj, aeroObj = self._getObjective(objective)

        # Note: Computeobjective partials MUST be called with the full
        # force pt list.
        if forcePoints is None:
            [npts, ncells, nTS] = self.sumb.getforcesize()
            forcePoints = numpy.zeros((nTS, npts, 3),self.dtype)
            self.sumb.getforcepoints(forcePoints.T)
        # end if

        if aeroObj:
            obj_num = self.SUmbCostfunctions[obj]

            if self.getOption('useReverseModeAD'):
                self.sumb.computeobjpartials(
                    obj_num, forcePoints.T, True, True)
            else:
                self.sumb.computeobjectivepartialsfwd(obj_num)
                

            self.adjointRHS = obj


        else:
            self.sumb.zeroobjpartials(True, True)
        # end if

        return 

    def getdIdx(self, objective, forcePoints=None, TS=0, group_name=None):

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

        if group_name is not None:
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

            dxs = self.mesh.getdXs(group_name)
            return dxs
        else:
            return dXv
        # end if
        
    def getdIda(self, objective, forcePoints=None):

        obj, aeroObj = self._getObjective(objective)

        if self.nDVAero > 0:
            
            self.computeObjPartials(objective, forcePoints)
            if aeroObj:
                dIda_local = self.sumb.adjointvars.dida
            else:
                dIda_local = numpy.zeros_like(self.sumb.adjointvars.dida)
            # end if

            # We must MPI all reuduce
            dIda = self.comm.allreduce(dIda_local,  op=MPI.SUM)
        else:
            dIda = numpy.zeros((0))
        # end if

        return dIda
        
    def getdIdw(self, dIdw, objective, forcePoints=None):

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

        nw     = self.sumb.flowvarrefstate.nw
        ncells = self.sumb.adjointvars.ncellslocal[0]
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return nw*ncells*ntime

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

    def setAdjoint(self, adjoint, objective=None):
        '''Sets the adjoint vector externally. Used in coupled solver'''
        self.sumb.setadjoint(adjoint)

        if objective is not None:
            obj, aeroObj = self._getObjective(objective)
            self.storedADjoints[obj]= adjoint.copy()
        
        return

    def getResidual(self, res=None):
        '''Return the residual on this processor. Used in aerostructural
        analysis'''
        if res is None:
            res = numpy.zeros(self.getStateSize())
        res = self.sumb.getres(res)
        
        return res

    def computedSdwTVec(self, in_vec, out_vec, group_name):
        '''This function computes: out_vec = out_vec + dFdw^T*in_vec'''
        phi = self.mesh.expandVectorByFamily(group_name, in_vec)
        out_vec = self.sumb.getdfdwtvec(numpy.ravel(phi), out_vec)

        return out_vec

    def getSolution(self, sps=1):
        ''' Retrieve the solution variables from the solver. Note this
        is a collective function and must be called on all processors
        '''

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

    def computeArea(self, axis, group_name=None, cfd_force_pts=None, TS=0):
        """
        Compute the projected area of the surface mesh

        Input Arguments:
           axis, numpy array, size(3): The projection vector
               along which to determine the shadow area
           group_name, str: The group from which to obtain the coordinates.
               This name must have been obtained from addFamilyGroup() or 
               be the default 'all' which contains all surface coordiantes 

        Output Arguments: 
            Area: The resulting area  
            """

        if cfd_force_pts is None:
            cfd_force_pts = self.getForcePoints(TS)
        # end if
        if len(cfd_force_pts) > 0:
            areas = self.sumb.getareas(cfd_force_pts.T, TS, axis).T
        else:
            areas = numpy.zeros((0,3), self.dtype)
        # end if

        if group_name is not None:
            areas = self.mesh.sectionVectorByFamily(group_name, areas)
        # end if

        # Now we do an mpiallreduce with sum:
        area = self.comm.allreduce(numpy.sum(areas), op=MPI.SUM)
        
        return area

    def computeAreaSensitivity(self, axis, group_name=None, 
                               cfd_force_pts=None, TS=0):
        """ 
        Compute the projected area of the surface mesh

        Input Arguments:
           axis, numpy array, size(3): The projection vector
               along which to determine the shadow area  
           group_name, str: The group from which to obtain the coordinates.
               This name must have been obtained from addFamilyGroup() or 
               be the default 'all' which contains all surface coordiantes 

        Output Arguments:
            Area: The resulting area    
            """

        if cfd_force_pts is None:
            cfd_force_pts = self.getForcePoints(TS)
        # end if
        if len(cfd_force_pts) > 0:
            da = self.sumb.getareasensitivity(cfd_force_pts.T, TS, axis).T
        else:
            da = numpy.zeros((0,3), self.dtype)
        # end if

        if group_name is not None:
            da = self.mesh.sectionVectorByFamily(group_name, da)
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

    def _on_setOption(self, name, value):
        '''
        Set Solver Option Value 
        '''
        # Ignored options do NOT get set in solver
        if name in self.ignore_options:
            return

        # Do special Options individually
        if name in self.special_options:
            if name in ['monitorvariables',
                        'surfacevariables',
                        'volumevariables']:
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
        exec_str = 'self.sumb.'+self.optionMap[name]['location'] + '=' + value

        exec(exec_str)

        return
        
    def getOption(self, name):

        # Redefine the getOption def from the base class so we can
        # mane sure the name is lowercase

        def_options = self.options['defaults']
        if def_options.has_key(name.lower()):
            return self.options[name.lower()][1]
        else:    
            raise InputError(repr(name) + ' is not a valid option name')
        #end
        
        return
  

class SUmbDummyMesh(object):
    """
    Represents a dummy Multiblock structured Mesh for SUmb
    """
 
    def __init__(self):
        return
        
    def getSurfaceCoordinates(self, group_name):
    
        return 

    def getSurfaceConnectivity(self, group_name):

        return

    def setSurfaceCoordinates(self, group_name, coordinates):
       
        return 

    def addFamilyGroup(self, group_name, families=None):

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

    def getdXs(self, group_name):

        return 

# ==========================================================================
#                        Output Functionality
# ==========================================================================

    def writeVolumeGrid(self, file_name):
        '''write volume mesh'''

        return

    def writeSurfaceGrid(self, file_name):
        '''write surface mesh'''

        return

    def writeFEGrid(self, file_name):
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

    def WarpDeriv(self, solver_dxv):

        return

#==============================================================================
# SUmb Analysis Test
#==============================================================================
if __name__ == '__main__':
    
    # Test SUmb
    mpiPrint('Testing ...')
    sumb = SUMB()
    print(sumb)

