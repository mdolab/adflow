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
import pickle

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
        
        # 
        name = 'SUMB'
        category = 'Three Dimensional CFD'
        def_opts = {
            'gridFile':[str,'default.cgns'],
            'restartFile':[str,'default_restart.cgns'],
            'Discretization':[str,'Central plus scalar dissipation'],
            'Smoother':[str,'Runge-Kutta'],
            'turblenceModel':[str,'SA'], 
            'useWallFunctions':[bool,False],
            'reynoldsNumber':[float,1e6], 
            'reynoldsLength':[float,1.0], 
            'wallTreatment':[str,'Linear Pressure Extrapolation'],
            'nCycles':[int,500],
            'CFL':[float,1.7],
            'CFLCoarse':[float,1.0], # New
            'equationType': [str,'Euler'],
            'equationMode': [str,'Steady'],
            'flowType':[str,'External'],
            'Mach':[float,0.6],
            'machCoef':[float,0.6],
            'machGrid':[float,0.0],
            'L2Convergence':[float,1e-6],
            'L2ConvergenceRel':[float,1e-16],
            'L2ConvergenceCoarse':[float,1e-2], 
            'MGCycle':[str,'3w'],
            'metricConversion':[float,1.0],
            'storeHistory':[bool,False],
            'numberSolutions':[bool,False],
            'refTemp.':[float,398.0],
            'refPressure':[float,17654.0],
            'refDensity':[float,.28837],
            'blockSplitting':[bool,False],
            'loadImbalance':[float,0.1],
            'dissipationScalingExponent':[float,0.67],
            'vis4':[float,0.0156],
            'vis2':[float,0.25],
            'vis2Coarse':[float,0.5], 
            'restrictionRelaxation':[float,1.0],
            'printIterations':[bool,False],
            'printTiming':[bool,True],
            'adjointL2Convergence':[float,1e-10],
            'adjointL2ConvergenceRel':[float,1e-16],
            'adjointL2ConvergenceAbs':[float,1e-16],
            'adjointDivTol':[float,1e5],
            'monitorVariables':[list,['cl','cd','cmz']],
            'probName':[str,''],
            'outputDir':[str,'./'],
            'solRestart':[bool,False],
            'solveAdjoint':[bool,False],
            'setMonitor':[bool,True],
            'writeSolution':[bool,True],
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
            'TSStability': [bool,False],
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
            'omegaFourier': [float,3.14],
            'fourierSineCoefficient':[float,3.52e-4],
            'surfaceVariables':[list,['cp','vx','vy','vz','mach']],
            'volumeVariables':[list,['resrho']]
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

        self.callCounter = 0
        self.meshWarpingInitialized = False
        self.allInitialized = False
        
        self.sumb.iteration.standalonemode = False
        self.sumb.iteration.deforming_grid = deforming_mesh = False
        # Write the intro message
        self.sumb.writeintromessage()

        # Set the frompython flag to true
        self.sumb.killsignals.frompython = True

        self.SUmbCostfunctions = \
            {'cl':self.sumb.adjointvars.costfuncliftcoef,
             'cd':self.sumb.adjointvars.costfuncdragcoef,
             'cFx':self.sumb.adjointvars.costfuncforcexcoef,
             'cFy':self.sumb.adjointvars.costfuncforceycoef,
             'cFz':self.sumb.adjointvars.costfuncforcezcoef,
             'cMx':self.sumb.adjointvars.costfuncmomxcoef,
             'cMy':self.sumb.adjointvars.costfuncmomycoef,
             'cMz':self.sumb.adjointvars.costfuncmomzcoef,
             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,
             'cM0':self.sumb.adjointvars.costfunccm0,
             'clAlpha':self.sumb.adjointvars.costfuncclalpha,
             'cl0':self.sumb.adjointvars.costfunccl0,
             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,
             'cd0':self.sumb.adjointvars.costfunccd0
             }
        
        self.possibleObjectives = \
            { 'lift':'cl','cl':'cl','Cl':'cl',
              'drag':'cd','cd':'cd','Cd':'Cd',
              'forcex':'cFx','xforce':'cFx','cfx':'cFx',
              'forcey':'cFy','yForce':'cFy','cfy':'cFy',
              'forcez':'cFz','zForce':'cFz','cfz':'cFz',
              'momentx':'cMx','xmoment':'cMx','cmx':'cMx',
              'momenty':'cMy','ymoment':'cMy','cmy':'cMy',
              'momentz':'cMz','zmoment':'cMz','cmz':'cMz',
              'cmzalpha':'cMzAlpha',
              'cm0':'cM0',
              'clalpha':'clAlpha',
              'cl0':'cl0',
              'cdalpha':'cdAlpha',
              'cd0':'cd0'
              }
        if 'optionMap' in kwargs:
            self.optionMap = kwargs['optionMap']
        else:
            
            self.optionMap = \
                {'Discretization':{'Central plus scalar dissipation':
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
                                 'Laminar NS':
                                     self.sumb.inputphysics.unsteady,
                                 'RANS':
                                     self.sumb.inputphysics.timespectral,
                                 'location':
                                     'inputphysics.equationmode'},
                 'flowType':{'Internal':
                                 self.sumb.inputphysics.internalflow,
                             'External':
                                 self.sumb.inputphysics.externalflow,
                             'location':
                                 'inputphysics.flowtype'},
                 
                 'MGCycle':{'location':'localmg.mgdescription',
                            'len':self.sumb.constants.maxstringlen},
                 'gridFile':{'location':'inputio.gridfile',
                             'len':self.sumb.constants.maxstringlen},
                 'restartFile':{'location':'inputio.restartfile',
                                'len':self.sumb.constants.maxstringlen},
                 'useWallFunctions':{'location':'inputphysics.wallfunctions'},
                 'reynoldsNumber':{'location':'inputphysics.reynolds'},
                 'reynoldsLength':{'location':'inputphysics.reynoldslength'},
                 'solRestart':{'location':'inputio.restart'},
                 'nCycles':{'location':'inputiteration.ncycles'},
                 'CFL':{'location':'inputiteration.cfl'},        
                 'CFLCoarse':{'location':'inputiteration.cflcoarse'},        
                 'Mach':{'location':'inputphysics.mach'},
                 'machCoef':{'location':'inputphysics.machcoef'},
                 'machGrid':{'location':'inputphysics.machgrid'},
                 'L2Convergence':{'location':'inputiteration.l2conv'},
                 'L2ConvergenceRel':{'location':'inputiteration.l2convrel'},
                 'L2ConvergenceCoarse':{'location':'inputiteration.l2convcoarse'},
                 'refTemp.':{'location':'flowvarrefstate.tref'},
                 'refPressure':{'location':'flowvarrefstate.pref'},
                 'refDensity':{'location':'flowvarrefstate.rhoref'},
                 'blockSplitting':{'location':'inputparallel.splitblocks'},
                 'loadImbalance':{'location':'inputparallel.loadimbalance'},
                 'dissipationScalingExponent':{'location':'inputdiscretization.adis'},
                 'vis4':{'location':'inputdiscretization.vis4'},
                 'vis2':{'location':'inputdiscretization.vis2'},
                 'vis2Coarse':{'location':'inputdiscretization.vis2coarse'},
                 'restrictionRelaxation':{'location':'inputiteration.fcoll'},
                 'printIterations':{'location':'inputiteration.printiterations'},
                 'adjointL2Convergence':{'location':'inputadjoint.adjreltol'},
                 'adjointL2ConvergenceRel':{'location':'inputadjoint.adjreltolrel'},
                 'adjointL2ConvergenceAbs':{'location':'inputadjoint.adjabstol'},
                 'adjointDivTol':{'location':'inputadjoint.adjdivtol'},
                 'solveAdjoint':{'location':'inputadjoint.solveadjoint'},
                 'setMonitor':{'location':'inputadjoint.setmonitor'},
                 'approxPC':{'location':'inputadjoint.approxpc'},
                 'restartAdjoint':{'location':'inputadjoint.restartadjoint'},
                 'ILUFill':{'location':'inputadjoint.filllevel'},
                 'ASMOverlap':{'location':'inputadjoint.overlap'},
                 'dissipationLumpingParameter':{'location':'inputadjoint.sigma'},
                 'adjointMaxIter':{'location':'inputadjoint.adjmaxiter'},
                 'adjointMonitorStep':{'location':'inputadjoint.adjmonstep'},
                 'adjointSubspaceSize':{'location':'inputadjoint.adjrestart'},
                 'printTiming':{'location':'inputadjoint.printtiming'},
                 'TSStability':{'location':'inputtsstabderiv.tsstability'},
                 'alphaMode':{'location':'inputtsstabderiv.tsalphamode'},
                 'betaMode':{'location':'inputtsstabderiv.tsbetamode'},
                 'machMode':{'location':'inputtsstabderiv.tsmachmode'},
                 'pMode':{'location':'inputtsstabderiv.tspmode'},
                 'qMode':{'location':'inputtsstabderiv.tsqmode'},
                 'rMode':{'location':'inputtsstabderiv.tsrmode'},
                 'altitudeMode':{'location':'inputtsstabderiv.tsaltitudemode'},
                 'windAxis':{'location':'inputtsstabderiv.usewindaxis'},
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
                 'rotCenter':{'location':'inputmotion.rotPoint'},
                 'timeIntervals':{'location':'inputtimespectral.ntimeintervalsspectral'},
                 }
        # end if

        if 'ignore_options' in kwargs:
            self.ignore_options = kwargs['ignore_options']
        else:
            self.ignore_options = ['defaults',
                                   'storeHistory',
                                   'numberSolutions',
                                   'writeSolution',
                                   'familyRot',  # -> Not sure how to do
                                   'rotRate',    # -> Not sure how to do
                                   'omegaFourier', # -> Not sure how to do
                                   'fourierSineCoefficient' # -> Not sure how to do
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
        if 'mesh' in kwargs:
            self.mesh = kwargs['mesh']
        else:
            print 'Mesh must be specified'
            sys.exit(1)
        # end if
            
            
        # Write the intro message
        self.sumb.writeintromessage()
        return

    def initialize(self, aero_problem, sol_type, grid_file='default',
                   *args,**kwargs):
        '''
        Run High Level Initialization 
        
        Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
        '''

        if self.allInitialized == True:
            return
        
        # Make sure all the params are ok
        for option in self.options:
            if option != 'defaults':
                self.setOption(option, self.options[option][1])

        self.sumb.dummyreadparamfile()

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
        

        # Create Surface Node list
        if(self.myid==0): print ' -> Creating Surface Node List'
        self.sumb.mdcreatensurfnodes()

        # Setup External Warping
        meshInd = self.getMeshIndices()
        self.mesh.setExternalMeshIndices(meshInd)

        forceInd = self.getForceIndices()
        self.mesh.setExternalForceIndices(forceInd)
        
        #Set flags for ADjoint initialization
        self.adjointInitialized = False
        self.allInitialized = True

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
        self.sumb.inputmotion.rotpoint[0] = aero_problem._refs.xref\
            *self.metricConversion
        self.sumb.inputmotion.rotpoint[1] = aero_problem._refs.yref\
            *self.metricConversion
        self.sumb.inputmotion.rotpoint[2] = aero_problem._refs.zref\
            *self.metricConversion
        #update the flow vars
        self.sumb.updatereferencepoint()

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
        q = aero_problem._flows.qhat*V/aero_problem._refs.cref
        r = aero_problem._flows.rhat*V/aero_problem._refs.bref

        self.sumb.updaterotationrate(p,r,q)

        return
    
    def setRefArea(self,aero_problem):
        self.sumb.inputphysics.surfaceref = aero_problem._refs.sref*self.metricConversion**2
        self.sumb.inputphysics.lengthref = aero_problem._refs.cref*self.metricConversion
        
        return

    def resetFlow(self):
        '''
        Reset the flow for the complex derivative calculation
        '''

        self.sumb.setuniformflow()

        return

    def __solve__(self, aero_problem, niterations, sol_type, 
                  grid_file='default', *args, **kwargs):
        
        '''
        Run Analyzer (Analyzer Specific Routine)
        
        Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
        '''


        self.initialize(aero_problem,sol_type,grid_file,*args,**kwargs)
        #set inflow angle
        self.setInflowAngle(aero_problem)
        self.setReferencePoint(aero_problem)
        self.setRotationRate(aero_problem)
        self.setRefArea(aero_problem)
        # Run Solver
    
        t0 = time.time()
        self.sumb.monitor.niterold = self.sumb_comm_world.bcast(
            self.sumb.monitor.niterold, root=0)

        ncycles = niterations
        storeHistory = self.getOption('storeHistory')

        #reset python failure check to false
        self.sumb.killsignals.routinefailed=False

        if sol_type.lower() in ['steady', 'time spectral']:
            
            #set the number of cycles for this call
            self.sumb.inputiteration.ncycles = niterations
            
            if (self.sumb.monitor.niterold == 0 and 
                self.sumb.monitor.nitercur == 0 and 
                self.sumb.iteration.itertot == 0):

                # No iterations have been done
                if (self.sumb.inputio.storeconvinneriter):
                    nn = self.sumb.inputiteration.nsgstartup + \
                        self.sumb.inputiteration.ncycles
                    if(self.myid==0):
                        self.sumb.deallocconvarrays()
                        self.sumb.allocconvarrays(nn)
                    # end if
                    # end if

            elif(self.sumb.monitor.nitercur == 0 and  
                 self.sumb.iteration.itertot == 0):
                # Reallocate convergence history array and
                # time array with new size, storing old values from restart
                if (self.myid == 0):
                    # number of time steps from restart
                    ntimestepsrestart = self.sumb.monitor.ntimestepsrestart
                    
                    if (self.sumb.inputio.storeconvinneriter):
                        # number of iterations from restart
                        niterold = self.sumb.monitor.niterold#[0]
                        if storeHistory:
                            # store restart convergence history and
                            # deallocate array
                            temp = copy.deepcopy(
                                self.sumb.monitor.convarray[:niterold+1,:])
                            self.sumb.deallocconvarrays()
                            # allocate convergence history array with
                            # new extended size
                            self.sumb.allocconvarrays(
                                temp.shape[0]+self.sumb.inputiteration.ncycles-1)
                            # recover values from restart and
                            # deallocate temporary array
                            self.sumb.monitor.convarray[\
                                :temp.shape[0],:temp.shape[1]] = temp
                            temp = None
                        else:
                            temp=copy.deepcopy(
                                self.sumb.monitor.convarray[0,:,:])
                            self.sumb.deallocconvarrays()
                            # allocate convergence history array with
                            # new extended size
                            self.sumb.allocconvarrays(
                                self.sumb.inputiteration.ncycles+1+niterold)
                            self.sumb.monitor.convarray[0,:,:] = temp
                            self.sumb.monitor.convarray[1,:,:] = temp
                            temp = None
                        #endif
                    #endif
                #endif
            else:

                # More Time Steps / Iterations in the same session
                # Reallocate convergence history array and time array
                # with new size, storing old values from previous runs
                if (self.myid == 0):
                    if (self.sumb.inputio.storeconvinneriter):
                        if storeHistory:
                            # store previous convergence history and
                            # deallocate array
                            temp = copy.deepcopy(self.sumb.monitor.convarray)

                            self.sumb.deallocconvarrays()
                            # allocate convergence history array with
                            # new extended size
                            nn = self.sumb.inputiteration.nsgstartup + \
                                self.sumb.inputiteration.ncycles
                            
                            self.sumb.allocconvarrays(temp.shape[0]+nn-1)
                        
                            # recover values from previous runs and
                            # deallocate temporary array
                            self.sumb.monitor.convarray[:temp.shape[0],:] = \
                                copy.deepcopy(temp)
                            
                            temp = None
                        else:
                            temp=copy.deepcopy(
                                self.sumb.monitor.convarray[0,:,:])
                            self.sumb.deallocconvarrays()
                            self.sumb.allocconvarrays(
                                self.sumb.inputiteration.ncycles+1)
                            self.sumb.monitor.convarray[0,:,:] = temp
                            temp = None
                        #endif
                    #endif
                #endif

                # re-initialize iteration variables
                self.sumb.inputiteration.mgstartlevel = 1

                if storeHistory:
                    self.sumb.monitor.niterold  = self.sumb.monitor.nitercur
                else:
                    self.sumb.monitor.niterold  = 1#0#self.sumb.monitor.nitercur
                #endif

                self.sumb.monitor.nitercur  = 0#1
                self.sumb.iteration.itertot = 0#1
                
                # update number of time steps from restart
                self.sumb.monitor.ntimestepsrestart = \
                    self.sumb.monitor.ntimestepsrestart + \
                    self.sumb.monitor.timestepunsteady

                # re-initialize number of time steps previously run
                # (excluding restart)
                self.sumb.monitor.timestepunsteady = 0

                # update time previously run
                self.sumb.monitor.timeunsteadyrestart = \
                    self.sumb.monitor.timeunsteadyrestart + \
                    self.sumb.monitor.timeunsteady

                # re-initialize time run
                self.sumb.monitor.timeunsteady = 0.0
                
            #endif
        elif sol_type.lower()=='unsteady':
            print 'unsteady not implemented yet...'
            sys.exit(0)
        #endif

        self._UpdateGeometryInfo()
        self.routineFailed = self.sumb_comm_world.allreduce(
            self.sumb.killsignals.routinefailed,MPI.MIN)
        
        if (abs(self.routineFailed)==True):
            if self.myid ==0:print 'Error raise in updateGeometry'
            raise ValueError
        #endif

        # Now call the solver
        
        self.sumb.solver()

        #Check to see whether we have a valid solution
        #in this case routineFailed will be triggered on all processors
        self.routineFailed = self.sumb_comm_world.allreduce(
            abs(self.sumb.killsignals.routinefailed),MPI.SUM)/self.comm.size
        #therefore sum up and devide by nProc
        if (abs(self.routineFailed)==True):
            if self.myid ==0: print 'Error Raised in solver'
            raise ValueError
        #endif

        sol_time = time.time() - t0

        if self.getOption('printTiming'):
            if(self.myid==0):
                print 'Solution Time',sol_time
        
        # Post-Processing
        #Write solutions
        if self.getOption('writeSolution'):
            volname=self.OutputDir+self.probName+self.filename+'_vol.cgns'
            surfname=self.OutputDir+self.probName+self.filename+'_surf.cgns'

            if self.getOption('numberSolutions'):
                volname=self.OutputDir+self.probName+self.filename+'_vol%d.cgns'%(self.callCounter)
                surfname=self.OutputDir+self.probName+self.filename+'_surf%d.cgns'%(self.callCounter)
            #endif
        
            if(self.myid==0):print volname,surfname
            self.WriteVolumeSolutionFile(volname)
            self.WriteSurfaceSolutionFile(surfname)
        # end if
        # get forces? from SUmb attributes
        if self.getOption('TSStability'):
            self.computeStabilityParameters()
        #endif
        
        self.callCounter+=1

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

    def getResiduals(self):
        if self.myid == 0:
            return self.sumb.monitor.convarray[0,0,0],\
                self.sumb.monitor.convarray[1,0,0],\
                self.sumb.monitor.convarray[self.sumb.monitor.nitercur,0,0]
        else:
            return None,None,None
        # end if

        
    def WriteVolumeSolutionFile(self,filename=None,gridname=None):
        """Write the current state of the volume flow solution to a CGNS file.
        
        Keyword arguments:
        
        filename -- the name of the file (optional)

        """

        self.sumb.monitor.writegrid = False
        if (filename):
            self.sumb.inputio.solfile[:] = ''
            self.sumb.inputio.solfile[0:len(filename)] = filename
        if (gridname):
            self.sumb.inputio.newgridfile[:] = ''
            self.sumb.inputio.newgridfile[0:len(gridname)] = gridname
            self.sumb.monitor.writegrid = True


        self.sumb.monitor.writevolume=True
        self.sumb.monitor.writesurface=False
        self.sumb.writesol()

    def WriteSurfaceSolutionFile(self,*filename):
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

    def getForces(self,group_name,cfd_force_pts=None):
        ''' Return the forces on this processor. Use
        cfd_force_pts to compute the forces if given
        
        '''
        if cfd_force_pts==None:
            pts = self.getForcePoints()
        # end if
        if len(cfd_force_pts) > 0:
            forces = self.sumb.getforces(cfd_force_pts.T).T
        else:
            forces = numpy.empty([0],dtype=self.dtype)
        # end if

        return self.mesh.solver_to_warp_force(group_name,forces)

    def getForcePoints(self):
        npts = self.sumb.getforcesize()
        if npts > 0:
            return self.sumb.getforcepoints(npts).T
        else:
            return numpy.empty([0],dtype=self.dtype)
        # end if


    def initAdjoint(self, *args, **kwargs):
        '''
        Initialize the Ajoint problem for this test case
        in SUMB
        '''

        
        #Set the mesh level and timespectral instance for this
        #computation
        self.level = 1
        self.sps = 1
        
        self.sumb.iteration.currentlevel=1
        self.sumb.iteration.groundlevel=1

        #Check to see if initialization has already been performed
        if(self.adjointInitialized):
            return
        
        #Run the preprocessing routine. Sets the node numbering and
        #allocates memory.
        self.sumb.preprocessingadjoint(self.level)
        
        #Initialize the design variable and function storage
        self.sumb.designinit()
        self.sumb.initializepetsc()
        self.sumb.createpetscvars()
        self.nSpatial = self.sumb.adjointvars.ndesignspatial

        self.adjointInitialized = True
        if(self.myid==0):
            print 'ADjoint Initialized Succesfully...'
        #endif

        return

    def setupAdjoint(self, forcePoints=None, **kwargs):
        '''
        Setup the adjoint matrix for the current solution
        '''
        self.sumb.setupallresidualmatrices(self.level)
        if forcePoints is None:
            forcePoints = self.getForcePoints()
        # end if

        self.sumb.setupcouplingmatrixstruct(forcePoints.T)
        self.sumb.setuppetscksp(self.level)
        self.mesh.setupWarpDeriv()

        return


    def releaseAdjointMemeory(self):
        '''
        release the KSP memory...
        '''

        return

    def _on_adjoint(self,objective,*args,**kwargs):

        # Setup the RHS
        obj = self.SUmbCostfunctions[self.possibleObjectives[objective.lower()]]
        self.sumb.setupadjointrhs(self.level,obj)

        if 'structAdjoint' in kwargs and 'group_name' in kwargs:
            group_name = kwargs['group_name']
            phi = kwargs['structAdjoint']
            solver_phi = self.mesh.warp_to_solver_force(group_name,phi)
            self.sumb.agumentrhs(phi)
        # end if

        restart = self.getOption('restartAdjoint')

        nw = self.sumb.flowvarrefstate.nw
        obj = self.possibleObjectives[objective.lower()]

        # If we have saved adjoints, 
        if restart:
            # Objective is already stored, so just set it
            if obj in self.storedADjoints.keys():
                self.sumb.setadjoint(self.storedADjoints[obj])
            else:
                # Objective is not yet run, allocated zeros and set
                self.storedADjoints[obj]=numpy.zeros([self.sumb.adjointvars.ncellslocal*nw],float)
                self.sumb.setadjoint(self.storedADjoints[obj])
            # end if
            # end if

        self.sumb.solveadjointtransposepetsc()

        if restart:
            self.storedADjoints[obj] =  self.sumb.getadjoint(self.sumb.adjointvars.ncellslocal*nw)
        # end if

        return

    def totalSurfaceDerivative(self,objective):
        # The adjoint vector is now calculated so perform the
        # following operatino to produce dI/dX_surf:
        # p represents partial, d total
        # dI/dX_s = (dXv/dXs)^T * (pI/pX_volume + dRdX_v^T * psi)

        restart = self.getOption('restartAdjoint')
        obj = self.possibleObjectives[objective.lower()]
        
        if restart: # Selected stored adjoint
            self.sumb.setadjoint(self.storedADjoints[obj])

        ndof = self.sumb.adjointvars.nnodeslocal*3
        dxv_solver0 = self.sumb.getdrdxvpsi(ndof)

        obj = self.SUmbCostfunctions[obj]
        self.sumb.setupgradientrhsvolume(obj)
        dxv_solver1 = self.sumb.getdjdx(ndof)
        self.mesh.WarpDeriv(dxv_solver0+dxv_solver1)
        dIdXs = self.mesh.getdXs('all')

        return dIdXs

    def verifyPartials(self):
        '''
        Run solverADjoint to verify the partial derivatives in the ADjoint
        '''
        self.sumb.verifydcfdx(1)

        return
    
    def computeStabilityParameters(self):
        '''
        run the stability derivative driver to compute the stability parameters
        from the time spectral solution
        '''

        self.sumb.stabilityderivativedriver()

        return

    def getFunctionValues(self):
        '''
        retrieve the solution values from SUmb
        '''
        if self.myid==0:print 'interface getting solution'
        # Map cost functions
        self.sumb.getsolution()

        SUmbsolutions = {'cl':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncliftcoef-1],\
                             'cd':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncdragcoef-1],\
                             'cFx':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncforcexcoef-1],\
                             'cFy':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncforceycoef-1],\
                             'cFz':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncforcezcoef-1],\
                             'cMx':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncmomxcoef-1],\
                             'cMy':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncmomycoef-1],\
                             'cMz':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncmomzcoef-1],\
                             'cMzAlpha':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccmzalpha-1],\
                             'cM0':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccm0-1],\
                             'clAlpha':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncclalpha-1],\
                             'cl0':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccl0-1],\
                             'cdAlpha':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccdalpha-1],\
                             'cd0':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccd0-1]
                         }

        return SUmbsolutions
    
    def _UpdateGeometryInfo(self):
        """Update the SUmb internal geometry info, if necessary."""
        if (self._update_geom_info):
            self.mesh.warpMesh()
            self.mesh.writeVolumeGrid('debug.cgns')
            self.sumb.setgrid(self.mesh.getSolverGrid())
            self.sumb.updatecoordinatesalllevels()
            self.sumb.updatewalldistancealllevels()
            self.sumb.updateslidingalllevels()
            self.sumb.updatemetricsalllevels()
            self.sumb.updategridvelocitiesalllevels()
            self._update_geom_info = False


    def GetMonitoringVariables(self):
        """Return a list of the text strings describing the variables being
        monitored.

        """
        return self.monnames.keys()
    

    def GetConvergenceHistory(self,name):
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

        if (finalIt < 0):
            fail = True
        else:
            fail = False
        # end if

        return startRes,finalRes,fail

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

    def finalizeAdjoint(self):
        '''
        destroy the PESTcKSP context
        '''
        self.meshDerivatives =[]
        self.releaseAdjointMemeory()
        
        return

    def getSolution(self):
        '''
        retrieve the solution variables from the solver.
        '''
        solution = self.getFunctionValues()
        
        return solution

#   #====================
#   #MD Coupling routines
#   #=====================
    
    def _on_setOption(self, name, value):
        
        '''
        Set Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 21, 2008 - Ruben E. Perez
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
            return
        # end if

        # All other options do genericaly by setting value in module:
        # Check if there is an additional mapping to what actually
        # has to be set in the solver
        try: 
            value = self.optionMap[name][value]
        except:
            pass

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

#==============================================================================
# SUmb Analysis Test
#==============================================================================
if __name__ == '__main__':
    
    # Test SUmb
    print 'Testing ...'
    sumb = SUMB()
    print sumb

