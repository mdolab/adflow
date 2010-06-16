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
	v. 1.0	- Original pyAero Framework Implementation (RP,SM 2008)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# SUmb Library
# =============================================================================
try:
	from sumbInterface import SUmbInterface,SUmbMesh
except ImportError:
	print 'Error: SUmbInterface failed to import'
#end

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
from numpy import real
# =============================================================================
# Extension modules
# =============================================================================
sys.path.append(os.path.abspath('../../../../pyACDT/Aerodynamics'))
from pyAero_solver import AeroSolver

# =============================================================================
# Misc Definitions
# =============================================================================


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
		}
		informs = {
		}
		AeroSolver.__init__(self, name, category, def_opts, informs, *args, **kwargs)
		
		#self.sumb = SUmbInterface()
		self.interface = SUmbInterface(*args,**kwargs)
		self.Mesh = self.interface.Mesh
		self.myid = self.interface.myid
		self.callCounter = 0
		self.volumeMatrixInitialized = False
		self.flowMatrixInitialized=False
		self.meshWarpingInitialized=False
		self.allInitialized = False
		# If a new option is added, make sure to add a default here!!
		self.solver_options_default={
			'probName':'',
			'OutputDir':'./',
			'Equation Type': 'Euler',\
			'reinitialize':False,
			'CFL':1.0,
			'L2Convergence':1e-6,
			'L2ConvergenceRel':0.0,
			'MGCycle':'sg',
			'MetricConversion':1.0,
			'Discretization':'Central plus scalar dissipation',
			'Dissipation Scaling Exponent':0.67,\
			'Dissipation Coefficients':[0.5,0.015625],\
			'sol_restart':'no',
			'Allow block splitting':'no',\
			'solveADjoint':'no',
			'set Monitor':'Yes',
			'printIterations':True,
			'printSolTime':True,
			'writeSolution':True,
			'Approx PC': 'no',
			'Adjoint solver type': 'GMRES',
			'adjoint relative tolerance':1e-10,
			'adjoint absolute tolerance':1e-16,
			'adjoint max iterations': 500,
			'adjoint restart iteration' : 80,
			'adjoint monitor step': 10,
			'dissipation lumping parameter':6,
			'Preconditioner Side': 'LEFT',
			'Matrix Ordering': 'NestedDissection',
			'Global Preconditioner Type': 'Additive Schwartz',
			'Local Preconditioner Type' : 'ILU',
			'ILU Fill Levels': 2,
			'ASM Overlap' : 5,
			'TS Stability': 'no',
			'Reference Temp.':398,
			'Reference Pressure':101325.0}

		return

	def initialize(self,aero_problem,sol_type,grid_file='default',*args,**kwargs):
		'''
		Run High Level Initialization 
		
		Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
		'''
	
		
		try:
			kwargs['solver_options'] = self._checkOptions(kwargs['solver_options'])
		except:
			kwargs['solver_options'] = self.solver_options_default
		# end try
		if self.allInitialized==True:return
		
		self.interface.initializeFlow(aero_problem,sol_type,grid_file, *args, **kwargs)
		self.filename=grid_file
		self.solver_options_default['reinitialize']=False

		#check if meshwarping is initialize. If not, do so
		if (not self.meshWarpingInitialized):
			if(self.interface.myid==0):print ' -> Initializing Internal Warping'
			self.interface.Mesh.initializeInternalWarping()
			self.meshWarpingInitialized=True
		#endif

		self.allInitialized = True

		return
		

	def __solve__(self, aero_problem, sol_type,grid_file='default', *args, **kwargs):
		
		'''
		Run Analyzer (Analyzer Specific Routine)
		
		Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
		'''

		self.initialize(aero_problem,sol_type,grid_file,*args,**kwargs)
		

			
		if kwargs['solver_options']['reinitialize'] == True:
			print 'reinitialization not yet implemented...'
			sys.exit(0)
			#self.interface.reinitialize()
			self.solver_options_default['reinitialize']=False

		# end if
		

		try:
			niterations = kwargs['niterations']
		except KeyError:
			test=1
		#endtry
	
		if(self.interface.myid==0):print ' ->Setting inflowangle'
		#set inflow angle
		self.interface.setInflowAngle(aero_problem)
		
		# Run Solver
		if(self.interface.myid==0):print ' ->Running iterations'
		# get flow and ref from aero_problem

		t0 = time.time()
		self.interface.RunIterations(sol_type=sol_type,\
					     ncycles=niterations,\
					     *args, **kwargs)
		sol_time = time.time() - t0
		if kwargs['solver_options']['printSolTime']:
			if(self.interface.myid==0):
				print 'Solution Time',sol_time
		
		# Post-Processing
		#Write solutions
## 		try:  kwargs['solver_options']['OutputDir']
## 		except KeyError:
## 			#volname=self.filename+'%d'%(self.callCounter)+'vol.cgns'
##                         #surfname=self.filename+'%d'%(self.callCounter)+'surf.cgns'
## 			volname=self.filename+'vol.cgns'
## 			surfname=self.filename+'surf.cgns'
## 		else:
## 			volname=kwargs['solver_options']['OutputDir']+self.filename+'vol.cgns'
## 			surfname=kwargs['solver_options']['OutputDir']+self.filename+'surf.cgns'
## 		#endif
		if kwargs['solver_options']['writeSolution']:
			volname=self.interface.OutputDir+self.interface.probName+self.filename+'_vol.cgns'
			surfname=self.interface.OutputDir+self.interface.probName+self.filename+'_surf.cgns'
		
			if(self.interface.myid==0):print volname,surfname
			self.interface.WriteVolumeSolutionFile(volname)
			self.interface.WriteSurfaceSolutionFile(surfname)
		# end if
		# get forces? from SUmb attributes
		try:  kwargs['solver_options']['TS Stability']
		except KeyError:
			test=1
		else:
			if (kwargs['solver_options']['TS Stability'])=='yes':
				self.interface.computeStabilityParameters()
 		        #endif
		#end
		
		# Store Results
		#aero_problem.addSol(self.__class__.__name__, sol_name, sol_time, sol_inform, 
		#	sol_geom, sol_flows, sol_options, display_opts=disp_opts,
		#	#Lift,Drag,CL,CD,CDi,CDp,CDw,Distribs,etc...
		#	arguments=args, **kwargs)
		
		self.callCounter+=1
		return

# 	def updateMesh(self,xyz,mapping={},meshwarping={}, *args, **kwargs):
# 		'''
# 		Take in a new surface and update the mesh
# 		'''
# 		# Compute the new CFD coordinates for each set of local CFD
# 		# Surface points
# 		#if mapping.cfd_surf_orig.shape[1]!=0:
# 		#	new_cfd_surf = mapping.getMappedSurface(xyz)
# 		##endif
# 		#Global surface defined on all processors
# 		new_cfd_surf = mapping.getMappedSurface(xyz)
		
	

# #		indices = self.Mesh.GetSurfaceIndicesLocal()
#                 #print 'surface indices',indices.shape,indices

# 		## # Translate new coordinates into perturbations
# ## 		if self.xyz_cfd_orig.shape[1]!=0:
# ## 			cfd_dispts = xyz_cfd_new - self.xyz_cfd_orig
# ## 		else:
# ## 			cfd_dispts = 0.0*self.xyz_cfd_orig
# ##                 #endif
# ##                 # what about case where initial block has no surface points but
# ##                 # association with another block causes perturbations? Handled
# ##                 # in warp!

# ##                 #return from perturbations to actual coordinates
# ## 	        cfd_xyz = self.xyz_cfd_orig + cfd_dispts
# ## 	        #print 'cfd',cfd_xyz
# ## 		# do we need this, can't we just pass the coordinates?
# 		[_parallel,mpi]=self.interface.CheckIfParallel()
		
# 		if _parallel :

# 			try:
				
# 				if(self.interface.myid==0):print 'Trying pySUmb internal meshwarping...'
# 				#Set the internal coordinates
# 				self.interface.Mesh.SetGlobalSurfaceCoordinates(new_cfd_surf)
# 				if(self.interface.myid==0):print 'Internal Surfaces set...'
				
# 				#now Warp the mesh
# 				self.interface.Mesh.warpMesh()
# 				if(self.interface.myid==0):print 'Internal warping Finished...'
				
# 			except:
# 				if(self.interface.myid==0):print 'pySUmb Internal meshwarping failed, using external...'
# 				sys.exit(0)
# 				# Update the Local blocks based on the perturbed points
# 				# on the local processor
# 				if mapping.cfd_surf_orig.shape[1]!=0:
# 					meshwarping.updateLocalBlockFaces(new_cfd_surf,indices)
# 				#endif
# 				mpi.WORLD.Barrier()
				
# 				# Propogate perturbations to all blocks based on
# 				# flow solver information
# 				meshwarping.synchronizeBlockFaces()
			
#                                 #print 'Block faces synchronized...'
# 				#Update the local blocks with the fortran routine WARPBLK
# 				meshwarping.updateLocalBlockCoords()

# 			#endtry
# ##	            #repeat
# ##             self.meshwarping.synchronizeBlockFaces()

# ##             print 'Block faces synchronized...'
# ##             #Update the local blocks with the fortran routine WARPBLK
# ##             self.meshwarping.updateLocalBlockCoords()
	

# 		else:
# 			ijk_blnum = self.Mesh.GetSurfaceIndices()
# 			print mpi.rank, "Done with self.Mesh.GetSurfaceIndices('')"
# 			if mpi.rank == 0:
# 				if(self.interface.myid==0):print "calling self.meshwarping.Warp ...."
# 				meshwarping.Warp(new_cfd_surf, ijk_blnum[3,:], ijk_blnum[0:3,:])
# 				if(self.interface.myid==0):print "done with self.meshwarping.Warp"
				
# 				for n in range(1, self.Mesh.GetNumberBlocks()+1):
# 					#print "In n loop:", mpi.rank, n
				 
# 					if mpi.rank == 0:
# 						# Get the new coordinates from warp
# 						[blocknum,ijkrange,xyz_new] = self.meshwarping.GetCoordinates(n)
# 					# end if
				
# 					if mpi.rank != 0:
# 						blocknum = None
# 						ijkrange = None
# 						xyz_new = None
# 					# end if
				
# 					blocknum = self.comm_world.bcast(blocknum) 
# 					ijkrange = self.comm_world.bcast(ijkrange) 
# 					xyz_new = self.comm_world.bcast(xyz_new) 
                    
# 					# Set new mesh coordinates in SUmb
# 					self.Mesh.SetCoordinates(blocknum, ijkrange, xyz_new.real)	
# 				# end for
# 			#endif
# 		#endif
# 		if self.myid==0: 
# 			print " New coordinates set in SUmb ... \n"
# 		#endif
# 		return

	def initAdjoint(self, *args, **kwargs):
		'''
		Initialize the Ajoint problem for this test case
		in SUMB
		'''
		self.interface.initializeADjoint()

		return

	def setupAdjointMatrix(self, *args, **kwargs):
		'''
		Setup the adjoint matrix for the current solution
		'''
		self.interface.setupADjointMatrix()

		return
		

	def _on_adjoint(self,objective,*args,**kwargs):

		if(self.interface.myid==0):print 'running cfd adjoint',objective

		self.interface.setupADjointRHS(objective)

		try:  kwargs['structAdjoint']
		except KeyError:
			test=1
		else:
			self.interface.augmentADjointRHS(objective,kwargs['structAdjoint'])
		#endtry

		self.interface.solveADjointPETSc()

		return

	def computeTotalFlowDerivatives(self,objective, *args, **kwargs):
		''' compute derivatives with respect to flow variables like alpha,Mach, beta....
		'''
		#Setup the partial derivative of the objective in sumb
		self.interface.setupGradientRHSFlow(objective)

		#setup the partial derivative of the volume coords. in sumb
		if not self.flowMatrixInitialized:
			self.interface.setupGradientMatrixFlow()
			self.flowMatrixInitialized=True
		#endif

		#compute and store the volume derivatives
		self.interface.computeTotalFlowDerivative(objective)

		#Retrieve a vector of the volume derivatives
		flowDerivative=self.interface.getTotalFlowDerivatives(objective)

		return flowDerivative

	def computeSurfaceDerivative(self, objective, *args, **kwargs):
		#def computeTotalSurfaceDerivative(self, objective,surface={},mapping={},meshwarping={}, *args, **kwargs):
		'''
		Compute the derivative of the objective function wrt the
		surface.
		'''
		#Setup the partial derivative of the objective in sumb
		self.interface.setupGradientRHSVolume(objective)

		#setup the partial derivative of the volume coords. in sumb
		if not self.volumeMatrixInitialized:
			self.interface.setupGradientMatrixVolume()
			self.volumeMatrixInitialized = True
		#endif

		#setup the derivative of the vol. coords. wrt. the surf.
		#coords. in SUmb (meshwarping...)
		self.interface.setupVolumeSurfaceDerivatives()

		self.interface.computeTotalVolumeDerivative(objective)
		#compute and store the surface derivatives
		self.interface.computeTotalSurfaceDerivative(objective)

		#Retrieve a vector of the volume derivatives
		surfaceDerivative=self.interface.getTotalSurfaceDerivatives(objective)
	## 	print surfaceDerivative.shape
## 		for i in xrange(len(surfaceDerivative[:])):
## 			print 'surfaceDerivative:',i,surfaceDerivative[i],i/3
## 			#endfor
## 		#endfor
		#print surfaceDerivative
		#stop
		#sys.exit(0)

		return surfaceDerivative*self.interface.Mesh.metricConversion

#*****************
# Depricated!!
#*****************
## 	def computeSurfaceDerivative(self, objective,surface={},mapping={},meshwarping={}, *args, **kwargs):
## 		#def computeTotalSurfaceDerivative(self, objective,surface={},mapping={},meshwarping={}, *args, **kwargs):
## 		'''
## 		Compute the derivative of the objective function wrt the
## 		surface.
## 		'''
## 		#Setup the partial derivative of the objective in sumb
## 		self.interface.setupGradientRHSVolume(objective)

## 		#setup the partial derivative of the volume coords. in sumb
## 		if not self.volumeMatrixInitialized:
## 			self.interface.setupGradientMatrixVolume()
## 			self.volumeMatrixInitialized = True
## 		#endif

## 		#compute and store the volume derivatives
## 		self.interface.computeTotalVolumeDerivative(objective)

## 		#Retrieve a vector of the volume derivatives
## 		volumeDerivative=self.interface.getTotalVolumeDerivatives(objective)
## 		#print volumeDerivative
## 		#stop

## 		#Get the global node ordering
## 		self.getGlobalNodeOrder(meshwarping=meshwarping)

## 		[xyz,conn,elemtype] = surface.getSurface()
## 		#now determine the surface derivatives using CS
## 		xyzref = copy.deepcopy(xyz)

## 		#initialize vector for the surface derivatives
## 		self.dIdxyz = numpy.zeros([len(xyzref[:,0]),len(xyzref[0,:])],'d')

## 		try: self.meshDerivatives
## 		except:
## 			self.meshDerivatives = []
## 		else:
## 			print 'meshDerivatives exists'
## 		#endif

## 		if (self.meshDerivatives == []):
## 			#it is necessary to compute the mesh derivatives

## 			#Setup Complex dummy array for the surface
## 			xyz_comp = numpy.zeros([len(xyz[:,0]),len(xyz[0,:])],'D')
		
			
## 			for i in xrange(len(xyzref[:,0])):
## 				if self.myid ==0:
## 					print "Coordinate %d of %d...."%(i,len(xyzref[:,0]))
## 				#endif

## 				#setup an empty list for this row
## 				rowDerivatives = []
## 				for j in xrange(len(xyzref[0,:])):
## 					if self.myid ==0:
## 						print 'j',j,'of',len(xyzref[0,:])
## 					#endif
## 					#set stepsize
## 					deltax = 1e-20j
					
##  				        #store reference design variables
## 					xref = xyz[i,j]
				
## 					#Copy design variables over to complex array
## 					xyz_comp[:,:] = xyz[:,:]
				
##    			                #perturb design variables
## 					xyz_comp[i,j] = xyz[i,j]+deltax
## 					#print 'xyz',xyz,'xyz_comp',xyz_comp

##                                         #Warp the mesh
## 					self.updateMesh(xyz_comp,mapping,meshwarping)
## 					#get the new Coordinates
## 					newMesh = meshwarping.getMeshCoordinates()
##  	                                #compute derivative
## 					newMeshDerivative = (newMesh.imag)/deltax.imag
## 					print 'newMeshDerivative',newMeshDerivative
## 					#append to mesh derivative list
## 					rowDerivatives.append(newMeshDerivative)


## 					#restore reference design variables
## 					xyz[i,j] = xref

## 				#endfor
## 				#append to mesh derivative list
## 				self.meshDerivatives.append(rowDerivatives)

## 			#endfor
## 		#endif
## 		timeA = time.time()
## 		for i in xrange(len(xyzref[:,0])):
## 			for j in xrange(len(xyzref[0,:])):
## 				#for k in xrange(len(volumeDerivative)):#self.meshDerivatives[i][j][:])):
## 				#print 'mesh derivatives',self.meshDerivatives[i][j][k],volumeDerivative[k]
## 				self.dIdxyz[i,j]= numpy.dot(self.meshDerivatives[i][j][:],volumeDerivative[:])
## 		        #endfor
			        
## 			#endfor
## 		#endfor				
## 		print 'Numpy dot products take:',time.time()-timeA,' seconds'
## 		return self.dIdxyz

# 	def computeAeroImplicitCoupling(self, objective,surface={},mapping={},meshwarping={}, *args, **kwargs):

# 		#setup the partial derivative of the volume coords. in sumb
# 		if not self.volumeMatrixInitialized:
# 			self.interface.setupGradientMatrixVolume()
# 			self.volumeMatrixInitialized = True
# 		#endif

# 		#compute and store the volume derivatives
# 		self.interface.computeAeroCouplingDerivative(objective)

# 		#Retrieve a vector of the volume derivatives
# 		couplingDerivative=self.interface.getAeroCouplingDerivatives(objective)
# 		#print volumeDerivative
# 		#stop

# 		#Get the global node ordering
# 		self.getGlobalNodeOrder(meshwarping=meshwarping)

# 		[xyz,conn,elemtype] = surface.getSurface()
# 		#now determine the surface derivatives using CS
# 		xyzref = copy.deepcopy(xyz)

# 		#initialize vector for the surface derivatives
# 		self.dJcdxyz = numpy.zeros([len(xyzref[:,0]),len(xyzref[0,:])],'d')
# 		if(self.interface.myid==0):print 'computeAeroImplicitCoupling'
# 		try: self.meshDerivatives
# 		except:
# 			self.meshDerivatives = []
# 		else:
# 			if(self.interface.myid==0):print 'meshDerivatives exists'
# 		#endif
# 		if self.meshDerivatives == []:
# 			filename = 'meshDerivatives+%04d'%self.myid+'.dat'
# 			try: 
# 				print 'trying to load mesh file:',filename
# 				pkl_file = open(filename, 'r')
# 				data = pickle.load(pkl_file)
# 				self.meshDerivatives= data['meshDerivatives']
# 			except:
# 				pass
# 			#end
# 		#end 
# 		#
	
	
# 		if (self.meshDerivatives == []):
# 			#it is necessary to compute the mesh derivatives

# 			#Setup Complex dummy array for the surface
# 			xyz_comp = numpy.zeros([len(xyz[:,0]),len(xyz[0,:])],'D')
		
			
# 			for i in xrange(xyzref.shape[0]):
# 				if self.myid ==0:
# 					print "Coordinate %d of %d...."%(i,xyzref.shape[0])
# 				#endif

# 				#setup an empty list for this row
# 				rowDerivatives = []
# 				for j in xrange(len(xyzref[0,:])):

# 					#set stepsize
# 					deltax = 1e-20j
					
#  				        #store reference design variables
# 					xref = xyz[i,j]
				
# 					#Copy design variables over to complex array
# 					xyz_comp[:,:] = xyz[:,:]
				
#    			                #perturb design variables
# 					xyz_comp[i,j] = xyz[i,j]+deltax
# 					#print 'xyz',xyz,'xyz_comp',xyz_comp

#                                         #Warp the mesh
# 					self.updateMesh(xyz_comp,mapping,meshwarping)
# 					#get the new Coordinates
# 					newMesh = meshwarping.getMeshCoordinates()
#  	                                #compute derivative
# 					newMeshDerivative = (newMesh.imag)/deltax.imag
# 					#print 'newMeshDerivative',newMeshDerivative
# 					#append to mesh derivative list
# 					rowDerivatives.append(newMeshDerivative)

# 					#restore reference design variables
# 					xyz[i,j] = xref
					
# 				#endfor
# 				#append to mesh derivative list
# 				self.meshDerivatives.append(rowDerivatives)

# 			#endfor
# 			filename = 'meshDerivatives+%04d'%self.myid+'.dat'
# 			data_save = {'meshDerivatives':self.meshDerivatives}
# 			pck_file = open(filename,'w')
# 			pickle.dump(data_save, pck_file)
# 		#end if 
# 		#endif
		
		
# 		for i in xrange(len(xyzref[:,0])):
# 			for j in xrange(len(xyzref[0,:])):
# 				self.dJcdxyz[i,j] = numpy.dot(self.meshDerivatives[i][j][:],couplingDerivative[:])
# 				#print 'mesh derivatives',self.meshDerivatives[i][j][k],volumeDerivative[k]
# 				#self.dJcdxyz[i,j]=self.dJcdxyz[i,j]+self.meshDerivatives[i][j][k]*\
# 				# couplingDerivative[k]
# 			#end for
# 		#endfor				

# 		return self.dJcdxyz

# 	def computeAeroExplicitCoupling(self, objective,surface={},mapping={},meshwarping={}, *args, **kwargs):

# 		#compute and store the volume derivatives
# 		self.interface.computeAeroExplicitCouplingDerivative(objective)

# 		#Retrieve a vector of the volume derivatives
# 		couplingDerivative=self.interface.getAeroExplicitCouplingDerivatives(objective)
# 		#print volumeDerivative
# 		#stop

# 		#Get the global node ordering
# 		self.getGlobalNodeOrder(meshwarping=meshwarping)

# 		[xyz,conn,elemtype] = surface.getSurface()
# 		#now determine the surface derivatives using CS
# 		xyzref = copy.deepcopy(xyz)

# 		#initialize vector for the surface derivatives
# 		self.dIcdxyz = numpy.zeros([len(xyzref[:,0]),len(xyzref[0,:])],'d')

# 		try: self.meshDerivatives
# 		except:
# 			self.meshDerivatives = []
# 		else:
# 			print 'meshDerivatives exists'
# 		#endif

# 		if (self.meshDerivatives == []):
# 			#it is necessary to compute the mesh derivatives

# 			#Setup Complex dummy array for the surface
# 			xyz_comp = numpy.zeros([len(xyz[:,0]),len(xyz[0,:])],'D')
		
			
# 			for i in xrange(len(xyzref[:,0])):
# 				if self.myid ==0:
# 					print "Coordinate %d of %d...."%(i,len(xyzref[:,0]))
# 				#endif

# 				#setup an empty list for this row
# 				rowDerivatives = []
# 				for j in xrange(len(xyzref[0,:])):

# 					#set stepsize
# 					deltax = 1e-20j
					
#  				        #store reference design variables
# 					xref = xyz[i,j]
				
# 					#Copy design variables over to complex array
# 					xyz_comp[:,:] = xyz[:,:]
				
#    			                #perturb design variables
# 					xyz_comp[i,j] = xyz[i,j]+deltax
# 					#print 'xyz',xyz,'xyz_comp',xyz_comp

#                                         #Warp the mesh
# 					self.updateMesh(xyz_comp,mapping,meshwarping)
# 					#get the new Coordinates
# 					newMesh = meshwarping.getMeshCoordinates()
#  	                                #compute derivative
# 					newMeshDerivative = (newMesh.imag)/deltax.imag
# 					#print 'newMeshDerivative',newMeshDerivative
# 					#append to mesh derivative list
# 					rowDerivatives.append(newMeshDerivative)

# 					#restore reference design variables
# 					xyz[i,j] = xref
					
# 				#endfor
# 				#append to mesh derivative list
# 				self.meshDerivatives.append(rowDerivatives)

# 			#endfor
# 		#endif
		
# 		for i in xrange(len(xyzref[:,0])):
# 			for j in xrange(len(xyzref[0,:])):
# 				self.dIcdxyz[i,j]=numpy.dot(self.meshDerivatives[i][j][:],couplingDerivative[:])
# 			#endfor
#     		#endfor				

# 		return self.dIcdxyz

# 	def aeroComputeTotalDerivatveStruct(self,objective,structAdjoint={}):
# 		'''
# 		compute the force portion of the total structural derivative
# 		based on the coupled structural adjoint
# 		'''

# 		self.interface.aeroComputeTotalDerivatveStruct(objective,structAdjoint=structAdjoint)
		
# 		return

# 	def getTotalDerivativeStruct(self,objective,meshwarping={},mapping={},surface={}):
# 		#Retrieve a vector of the volume derivatives
# 		structDerivative=self.interface.getTotalStructDerivatives(objective)
# 		#print volumeDerivative
# 		#stop

# 		#Get the global node ordering
# 		self.getGlobalNodeOrder(meshwarping=meshwarping)

# 		[xyz,conn,elemtype] = surface.getSurface()
# 		#now determine the surface derivatives using CS
# 		xyzref = copy.deepcopy(xyz)

# 		#initialize vector for the surface derivatives
# 		self.dStdxyz = numpy.zeros([len(xyzref[:,0]),len(xyzref[0,:])],'d')

# 		try: self.meshDerivatives
# 		except:
# 			self.meshDerivatives = []
# 		else:
# 			print 'meshDerivatives exists'
# 		#endif

# 		if (self.meshDerivatives == []):
# 			#it is necessary to compute the mesh derivatives

# 			#Setup Complex dummy array for the surface
# 			xyz_comp = numpy.zeros([len(xyz[:,0]),len(xyz[0,:])],'D')
		
			
# 			for i in xrange(len(xyzref[:,0])):
# 				if self.myid ==0:
# 					print "Coordinate %d of %d...."%(i,len(xyzref[:,0]))
# 				#endif

# 				#setup an empty list for this row
# 				rowDerivatives = []
# 				for j in xrange(len(xyzref[0,:])):
					
# 					#set stepsize
# 					deltax = 1e-20j
					
#  				        #store reference design variables
# 					xref = xyz[i,j]
				
# 					#Copy design variables over to complex array
# 					xyz_comp[:,:] = xyz[:,:]
				
#    			                #perturb design variables
# 					xyz_comp[i,j] = xyz[i,j]+deltax
# 					#print 'xyz',xyz,'xyz_comp',xyz_comp

#                                         #Warp the mesh
# 					self.updateMesh(xyz_comp,mapping,meshwarping)
# 					#get the new Coordinates
# 					newMesh = meshwarping.getMeshCoordinates()
#  	                                #compute derivative
# 					newMeshDerivative = (newMesh.imag)/deltax.imag
# 					#print 'newMeshDerivative',newMeshDerivative
# 					#append to mesh derivative list
# 					rowDerivatives.append(newMeshDerivative)

# 					#restore reference design variables
# 					xyz[i,j] = xref

# 				#endfor
# 				#append to mesh derivative list
# 				self.meshDerivatives.append(rowDerivatives)

# 			#endfor
# 		#endif
# 		for i in xrange(len(xyzref[:,0])):
# 			for j in xrange(len(xyzref[0,:])):
# 				self.dStdxyz[i,j] = numpy.dot(self.meshDerivatives[i][j][:],structDerivative[:])
# 			#endfor
# 		#endfor				

# 		return self.dStdxyz     


#	def solveCoupledAdjoint(self,objective,structsolver={},*args,**kwargs):
#
#
#		return
	
	def finalizeAdjoint(self):
		'''
		destroy the PESTcKSP context
		'''
		self.meshDerivatives =[]
		self.volumeMatrixInitialized = False
		self.flowMatrixInitialized=False	
		self.interface.releaseAdjointMemeory()
		
		return



# 	def getGlobalNodeOrder(self,meshwarping={}, *args, **kwargs):

# 		[_parallel,mpi]=self.interface.CheckIfParallel()
		
# 		if _parallel:
# 			#loop over the local blocks
# 			n = meshwarping.nBlocksLocal
# 			for i in xrange(n):
#                                 #retrieve the global node numbers from SUmb and store them in the meshwarping blocks
# 				ijk = meshwarping.blockList[i].ijk
# 				meshwarping.blockList[i].globalNode = self.interface.getGlobalNodesLocal(i+1,ijk[0],ijk[1],ijk[2])
			
                               
#                         #endfor
			
# 		else:
# 			if MPI.rank == 0:
# 				#assumes single processor is being used

#                                 #initialize size array and retrieve all block sizes
        
# 				ijkmax = numpy.zeros([3,self.interface.Mesh.GetNumberBlocks()+1],'i')
				
# 				for n in xrange(1, self.interface.Mesh.GetNumberBlocks()+1):
# 					ijkmax[:,n] = meshwarping.GetBlockDimensions(n)
# 				#endfor
            
#                                 #initialize the global node array
# 				ijkmax[0,0] = max(ijkmax[0,:])
# 				ijkmax[1,0] = max(ijkmax[1,:])
# 				ijkmax[2,0] = max(ijkmax[2,:])
                
#                                 #initialize the global node array
        
# 				self.globalNodes = numpy.zeros([ijkmax[0,0]+1,ijkmax[1,0]+1,ijkmax[2,0]+1,\
# 								self.interface.Mesh.GetNumberBlocks()+1],'i')
            
# 				for n in xrange(1, self.interface.Mesh.GetNumberBlocks()+1):
# 					q = int(ijkmax[0,n])
                                        
# 					self.interface.getGlobalNodes(n,ijkmax[0,n],ijkmax[1,n],ijkmax[2,n],\
# 								 self.globalNodes)
                                        
#                                 #endfor
#                         #endif
# 		#endif
# 		return

	def getSolution(self):
		'''
		retrieve the solution variables from the solver.
		'''
		solution = self.interface.getFunctionValues()
		
		return solution

# 	#====================
# 	#MD Coupling routines
# 	#=====================

# 	def getOMLForces(self,mapping={}):
# 		'''
# 		Compute the forces on the nodes and transfer them to the OML
# 		'''
# 		cfd_loads = self.interface.GetSurfaceLoadsLocal()
# 		cfdloads2 = self.interface.GetSurfaceLoads()

# #		print 'cfd_loads:'
# #		for i in xrange(cfd_loads.shape[0]):
# #			for j in xrange(cfd_loads.shape[1]):
# #				print cfd_loads[i][j]


# #               self.cfdloads2 = cfdloads2
# 		self.cfdloads2 = cfd_loads
# 		if self.myid == 0: 'OML Surface Shape',mapping.oml_surf_orig.shape
# 		oml_loads_local = numpy.zeros((3, len(mapping.oml_surf_orig[1])),'d')

		
# 		for icfd in range(cfd_loads.shape[1]):	
# 			ielem = mapping.nearest_ele[icfd] - 1
# 			for i in range(mapping.conn_oml.shape[0]):
# 				inode = int(mapping.conn_oml[i, ielem]) - 1
# 				weightm = mapping.weightt[:,:,i,icfd] + mapping.weightr[:,:,i,icfd]
# 				oml_loads_local[:,inode] = oml_loads_local[:,inode] + numpy.dot(numpy.transpose(weightm[:,:]), cfd_loads[:,icfd])
# 				#for k in xrange(weightm.shape[0]):
# 				#	for l in xrange(weightm.shape[1]):
# 				#		print real(weightm[k,l])
						
# 			#endfor
# 		#endfor
# 		#print 'local',oml_loads_local
# 		oml_loads = self.interface.AccumulateLoads(oml_loads_local)
# 		#print 'global',oml_loads
		
#                 mapping._WriteVectorsTecplot("oml_forces.dat", mapping.oml_surf_orig, oml_loads, 1.E-4)
# 		mapping._WriteVectorsTecplot("cfd_forces.dat", mapping.cfd_surf_orig, cfdloads2, 1.E-4)
# 		return oml_loads

# 	def setOMLDisplacements(self,oml_disp,mapping={},meshwarping={},surface={}):
# 		'''
# 		Compute the forces on the nodes and transfer them to the OML
# 		'''
# 		cfd_dispts = numpy.zeros((3, mapping.cfd_surf_orig.shape[1]), 'd')

# 		for icfd in range(cfd_dispts.shape[1]):	
# 			ielem = mapping.nearest_ele[icfd] - 1
# 			for i in range(mapping.conn_oml.shape[0]):
# 				inode = int(mapping.conn_oml[i, ielem]) - 1
# 				weightm = mapping.weightt[:,:,i,icfd] + mapping.weightr[:,:,i,icfd]
# 				cfd_dispts[:,icfd] = cfd_dispts[:,icfd] + numpy.dot(weightm[:,:], oml_disp[:,inode])
# 			#endfor
# 		#endfor

# 		self.cfd_dispts= cfd_dispts
# 		#get the original surface mesh
# 		[oml_surf,oml_conn,oml_elemtype]= surface.getSurface()
# 		if mapping.cfd_surf_orig.shape[1]!=0:
# 			new_cfd_surf = mapping.getMappedSurface(oml_surf)
# 		        #new_cfd_surf = copy.deepcopy(mapping.cfd_surf_orig)
		
#                 	#add the displacements
# 			new_cfd_surf = new_cfd_surf+cfd_dispts
# 		#endif

# 		#get the indices of the surface nodes on this block
# 		indices = self.Mesh.GetSurfaceIndicesLocal()
		
# 		[_parallel,mpi]=self.interface.CheckIfParallel()
			
# 		if _parallel :
# 			# Update the Local blocks based on the perturbed points
# 			# on the local processor
# 			if mapping.cfd_surf_orig.shape[1]!=0:
# 				meshwarping.updateLocalBlockFaces(new_cfd_surf,indices)
# 		        #endif
# 			mpi.WORLD.Barrier()
			
# 			# Propogate perturbations to all blocks based on
# 			# flow solver information
# 			meshwarping.synchronizeBlockFaces()
			
#                         #print 'Block faces synchronized...'
# 			#Update the local blocks with the fortran routine WARPBLK
# 			meshwarping.updateLocalBlockCoords()

# ##	            #repeat
# ##             self.meshwarping.synchronizeBlockFaces()

# ##             print 'Block faces synchronized...'
# ##             #Update the local blocks with the fortran routine WARPBLK
# ##             self.meshwarping.updateLocalBlockCoords()
	

# 		else:
# 			ijk_blnum = self.Mesh.GetSurfaceIndices()
# 			print mpi.rank, "Done with self.Mesh.GetSurfaceIndices('')"
# 			if mpi.rank == 0:
# 				print "calling self.meshwarping.Warp ...."
# 				meshwarping.Warp(new_cfd_surf, ijk_blnum[3,:], ijk_blnum[0:3,:])
# 				print "done with self.meshwarping.Warp"
				
# 				for n in range(1, self.Mesh.GetNumberBlocks()+1):
# 					#print "In n loop:", mpi.rank, n
				 
# 					if mpi.rank == 0:
# 						# Get the new coordinates from warp
# 						[blocknum,ijkrange,xyz_new] = self.meshwarping.GetCoordinates(n)
# 					# end if
				
# 					if mpi.rank != 0:
# 						blocknum = None
# 						ijkrange = None
# 						xyz_new = None
# 					# end if
				
# 					blocknum = self.comm_world.bcast(blocknum) 
# 					ijkrange = self.comm_world.bcast(ijkrange) 
# 					xyz_new = self.comm_world.bcast(xyz_new) 
                    
# 					# Set new mesh coordinates in SUmb
# 					self.Mesh.SetCoordinates(blocknum, ijkrange, xyz_new.real)	
# 				# end for
# 			#endif
# 		#endif
# 		if self.myid==0: 
# 			print " Displacements set in SUmb ... \n"
# 		#endif
# 		return 
		
	def _on_setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
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

	def _checkOptions(self,solver_options):
		'''Check the solver options against the default ones
		and add option iff it is NOT in solver_options
		'''
		for key in self.solver_options_default.keys():
			if not(key in solver_options.keys()):
				solver_options[key] = self.solver_options_default[key]
			else:
				self.solver_options_default[key]=solver_options[key]	
			# end if
		# end for
		return solver_options
			
		
	


#==============================================================================
# SUmb Analysis Test
#==============================================================================
if __name__ == '__main__':
	
	# Test SUmb
	print 'Testing ...'
	sumb = SUMB()
	print sumb
	
