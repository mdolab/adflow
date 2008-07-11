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

	def __solve__(self, aero_problem, sol_type,reinitialize=False, updategeometry=False,niterations = 1,grid_file='default', *args, **kwargs):
		
		'''
		Run Analyzer (Analyzer Specific Routine)
		
		Documentation last updated:  July. 3, 2008 - C.A.(Sandy) Mader
		'''
		
		# Pre-Processing
		if reinitialize:
			self.sumb.initializeFlow(aero_problem,sol_type,grid_file)
		#endif

		if updategeometry:
			self.sumb.updateGeometry(geometry)
		#endif
		
		# Run Solver
		
		# get flow and ref from aero_problem
		#print 'niterations',niterations
		t0 = time.time()
		self.sumb.RunIterations(0,niterations)
		sol_time = time.time() - t0
		
		
		# Post-Processing
		
		# get forces? from SUmb attributes
		
		
		# Store Results
		#aero_problem.addSol(self.__class__.__name__, sol_name, sol_time, sol_inform, 
		#	sol_geom, sol_flows, sol_options, display_opts=disp_opts,
		#	#Lift,Drag,CL,CD,CDi,CDp,CDw,Distribs,etc...
		#	arguments=args, **kwargs)
		
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
# VORLIN Analysis Test
#==============================================================================
if __name__ == '__main__':
	
	# Test VORLIN
	print 'Testing ...'
	sumb = SUMB()
	print sumb
	
