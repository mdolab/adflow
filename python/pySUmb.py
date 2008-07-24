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

# =============================================================================
# External Python modules
# =============================================================================
import numpy

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
			
		self.sumb = SUmbInterface()
		self.Mesh = self.sumb.Mesh
		self.myid = self.sumb.myid
		self.callCounter = 0

	def __solve__(self, aero_problem, sol_type,grid_file='default', *args, **kwargs):
		
		'''
		Run Analyzer (Analyzer Specific Routine)
		
		Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
		'''
		
		# Pre-Processing
		try:  kwargs['reinitialize']
		except KeyError:
			test=1
		else:
			if kwargs['reinitialize']:
				self.sumb.initializeFlow(aero_problem,sol_type,grid_file)
				self.filename=grid_file
			#endif
		#endtry
		try:  kwargs['niterations']
		except KeyError:
			test=1
		else:
			niterations= kwargs['niterations']
		#endtry

		#if updategeometry:
		#	self.sumb.updateGeometry(geometry)
		##endif
		
		# Run Solver
		
		# get flow and ref from aero_problem
		#print 'niterations',niterations
		t0 = time.time()
		self.sumb.RunIterations(0,niterations)
		sol_time = time.time() - t0
		
		
		# Post-Processing
		#Write solutions
		volname=self.filename+'%d'%(self.callCounter)+'vol.cgns'
		surfname=self.filename+'%d'%(self.callCounter)+'surf.cgns'
		print volname,surfname
		self.sumb.WriteVolumeSolutionFile(volname)
		self.sumb.WriteSurfaceSolutionFile(surfname)
		
		# get forces? from SUmb attributes
		
		
		# Store Results
		#aero_problem.addSol(self.__class__.__name__, sol_name, sol_time, sol_inform, 
		#	sol_geom, sol_flows, sol_options, display_opts=disp_opts,
		#	#Lift,Drag,CL,CD,CDi,CDp,CDw,Distribs,etc...
		#	arguments=args, **kwargs)
		self.callCounter+=1
		return

	def updateMesh(self,xyz,mapping={},meshwarping={}, *args, **kwargs):
		'''
		Take in a new surface and update the mesh
		'''
		# Compute the new CFD coordinates for each set of local CFD
		# Surface points
		if mapping.cfd_surf_orig.shape[1]!=0:
			new_cfd_surf = mapping.getMappedSurface(xyz)
		#endif

		indices = self.Mesh.GetSurfaceIndicesLocal()
                #print 'surface indices',indices.shape,MPI.rank

		## # Translate new coordinates into perturbations
## 		if self.xyz_cfd_orig.shape[1]!=0:
## 			cfd_dispts = xyz_cfd_new - self.xyz_cfd_orig
## 		else:
## 			cfd_dispts = 0.0*self.xyz_cfd_orig
##                 #endif
##                 # what about case where initial block has no surface points but
##                 # association with another block causes perturbations? Handled
##                 # in warp!

##                 #return from perturbations to actual coordinates
## 	        cfd_xyz = self.xyz_cfd_orig + cfd_dispts
## 	        #print 'cfd',cfd_xyz
## 		# do we need this, can't we just pass the coordinates?
		[_parallel,mpi]=self.sumb.CheckIfParallel()
		
		if _parallel :
			# Update the Local blocks based on the perturbed points
			# on the local processor
			if mapping.cfd_surf_orig.shape[1]!=0:
			    meshwarping.updateLocalBlockFaces(new_cfd_surf,indices)
		        #endif
			mpi.WORLD.Barrier()
	            
			# Propogate perturbations to all blocks based on
			# flow solver information
			meshwarping.synchronizeBlockFaces()
			
                        #print 'Block faces synchronized...'
			#Update the local blocks with the fortran routine WARPBLK
			meshwarping.updateLocalBlockCoords()

##	            #repeat
##             self.meshwarping.synchronizeBlockFaces()

##             print 'Block faces synchronized...'
##             #Update the local blocks with the fortran routine WARPBLK
##             self.meshwarping.updateLocalBlockCoords()
	

		else:
			ijk_blnum = self.Mesh.GetSurfaceIndices()
			print mpi.rank, "Done with self.Mesh.GetSurfaceIndices('')"
			if mpi.rank == 0:
				print "calling self.meshwarping.Warp ...."
				meshwarping.Warp(new_cfd_surf, ijk_blnum[3,:], ijk_blnum[0:3,:])
				print "done with self.meshwarping.Warp"
				
				for n in range(1, self.Mesh.GetNumberBlocks()+1):
					#print "In n loop:", mpi.rank, n
				 
					if mpi.rank == 0:
						# Get the new coordinates from warp
						[blocknum,ijkrange,xyz_new] = self.meshwarping.GetCoordinates(n)
					# end if
				
					if mpi.rank != 0:
						blocknum = None
						ijkrange = None
						xyz_new = None
					# end if
				
					blocknum = self.comm_world.bcast(blocknum) 
					ijkrange = self.comm_world.bcast(ijkrange) 
					xyz_new = self.comm_world.bcast(xyz_new) 
                    
					# Set new mesh coordinates in SUmb
					self.Mesh.SetCoordinates(blocknum, ijkrange, xyz_new.real)	
				# end for
			#endif
		#endif
		if self.myid==0: 
			print " New coordinates set in SUmb ... \n"
		#endif
		return

	def initAdjoint(self, *args, **kwargs):
		'''
		Initialize the Ajoint problem for this test case
		in SUMB
		'''
		self.sumb.initializeADjoint()

		return

	def setupAdjointMatrix(self, *args, **kwargs):
		'''
		Setup the adjoint matrix for the current solution
		'''
		self.sumb.setupADjointMatrix()

		return
		

	def _on_adjoint(self,objective,*args,**kwargs):

		#print 'running adjoint',objective

		self.sumb.setupADjointRHS(objective)

		self.sumb.solveADjointPETSc()

		return

	def computeFlowDerivatives(self,objective,flowMatrixInitialized, *args, **kwargs):
		''' compute derivatives with respect to flow variables like alpha,Mach, beta....
		'''
		#Setup the partial derivative of the objective in sumb
		self.sumb.setupGradientRHSFlow(objective)

		#setup the partial derivative of the volume coords. in sumb
		if not flowMatrixInitialized:
			self.sumb.setupGradientMatrixFlow()
		#endif

		#compute and store the volume derivatives
		self.sumb.computeTotalFlowDerivative(objective)

		#Retrieve a vector of the volume derivatives
		flowDerivative=self.sumb.getTotalFlowDerivatives()

		return flowDerivative

	def computeSurfaceDerivative(self, objective,volumeMatrixInitialized,surface={},mapping={},meshwarping={}, *args, **kwargs):
		'''
		Compute the derivative of the objective function wrt the
		surface.
		'''
		#Setup the partial derivative of the objective in sumb
		self.sumb.setupGradientRHSVolume(objective)

		#setup the partial derivative of the volume coords. in sumb
		if not volumeMatrixInitialized:
			self.sumb.setupGradientMatrixVolume()
		#endif

		#compute and store the volume derivatives
		self.sumb.computeTotalVolumeDerivative(objective)

		#Retrieve a vector of the volume derivatives
		volumeDerivative=self.sumb.getTotalVolumeDerivatives(objective)

		[xyz,conn,elemtype] = surface.getSurface()
		#now determine the surface derivatives using CS
		xyzref = copy.deepcopy(xyz)

		#initialize vector for the surface derivatives
		dIdxyz = numpy.zeros([len(xyzref[:,0]),len(xyzref[0,:])],'d')

		#Setup Complex dummy array for the surface
		xyz_comp = numpy.zeros([len(xyz[:,0]),len(xyz[0,:])],'D')

		#Get the global node ordering
		self.getGlobalNodeOrder(meshwarping)

		for i in xrange(len(xyzref[:,0])):
			for j in xrange(len(xyzref[0,:])):
				#set stepsize
				deltax = 1e-20j
				
				#store reference design variables
				xref = xyz[i,j]

                                #Copy design variables over to complex array
				xyz_comp[:,:] = xyz[:,:]

			        #perturb design variables
				xyz_comp[i,j] = xyz[i,j]+deltax

                                #Warp the mesh
				self.updateMesh(xyz_comp,mapping,meshwarping)

  			        #get the new Coordinates
				newMesh = meshwarping.getMeshCoordinates()

 	                        #compute derivative
				newMeshDerivative = (newMesh.imag)/deltax.imag

				for k in xrange(len(newMeshDerivative)):
					dIdxyz[i,j]=dIdxyz[i,j]+newMeshDerivative[j]*\
					   volumeDerivative[j]
			        #endfor

			        #restore reference design variables
				xyz[i,j] = xref
			#endfor
		#endfor				

		return dIdxyz

	def getGlobalNodeOrder(self,meshwarping={}, *args, **kwargs):

		[_parallel,mpi]=self.sumb.CheckIfParallel()
		
		if _parallel:
			#loop over the local blocks
			n = meshwarping.nBlocksLocal
			for i in xrange(n):
                                #retrieve the global node numbers from SUmb and store them in the meshwarping blocks
				ijk = meshwarping.blockList[i].ijk
				meshwarping.blockList[i].globalNode = self.sumb.getGlobalNodesLocal(i+1,ijk[0],ijk[1],ijk[2])
			
                               
                        #endfor
			
		else:
			if MPI.rank == 0:
				#assumes single processor is being used

                                #initialize size array and retrieve all block sizes
        
				ijkmax = numpy.zeros([3,self.sumb.Mesh.GetNumberBlocks()+1],'i')
				
				for n in xrange(1, self.sumb.Mesh.GetNumberBlocks()+1):
					ijkmax[:,n] = meshwarping.GetBlockDimensions(n)
				#endfor
            
                                #initialize the global node array
				ijkmax[0,0] = max(ijkmax[0,:])
				ijkmax[1,0] = max(ijkmax[1,:])
				ijkmax[2,0] = max(ijkmax[2,:])
                
                                #initialize the global node array
        
				self.globalNodes = numpy.zeros([ijkmax[0,0]+1,ijkmax[1,0]+1,ijkmax[2,0]+1,\
								self.sumb.Mesh.GetNumberBlocks()+1],'i')
            
				for n in xrange(1, self.sumb.Mesh.GetNumberBlocks()+1):
					q = int(ijkmax[0,n])
                                        
					self.sumb.getGlobalNodes(n,ijkmax[0,n],ijkmax[1,n],ijkmax[2,n],\
								 self.globalNodes)
                                        
                                #endfor
                        #endif
		#endif
		return

	def getSolution(self):
		'''
		retrieve the solution variables from the solver.
		'''
		#print 'getting solution'
		solution = self.sumb.getFunctionValues()
		#print solution
		
		return solution

		
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
	


#==============================================================================
# SUmb Analysis Test
#==============================================================================
if __name__ == '__main__':
	
	# Test SUmb
	print 'Testing ...'
	sumb = SUMB()
	print sumb
	
