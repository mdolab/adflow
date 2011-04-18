#!/usr/local/bin/python
'''
pyAero_solver

Holds the Python Aerodynamic Analysis Classes (base and inherited).

Copyright (c) 2008 by Dr. Ruben E. Perez
All rights reserved. Not to be used for commercial purposes.
Revision: 1.1   $Date: 21/05/2008 21:00$


Developers:
-----------
- Dr. Ruben E. Perez (RP)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2008)
'''

__version__ = '$Revision: $'

'''
To Do:goto-line
	- 
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import pdb

# =============================================================================
# External Python modules
# =============================================================================
#import external

# =============================================================================
# Extension modules
# =============================================================================
from pyAero_problem import AeroProblem

# =============================================================================
# Misc Definitions
# =============================================================================



# =============================================================================
# AeroSolver Class
# =============================================================================
class AeroSolver(object):
	
	'''
	Abstract Class for Aerodynamic Solver Object
	'''
	
	def __init__(self, name, category={}, def_options={}, informs={}, *args, **kwargs):
		
		'''
		AeroSolver Class Initialization
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.name = name
		self.category = category
		self.options = {}
		self.options['defaults'] = def_options
		self.informs = informs
		
		# Initialize Options
		def_keys = def_options.keys()
		for key in def_keys:
			self.options[key] = def_options[key]
		#end
		koptions = kwargs.pop('options',{})
		kopt_keys = koptions.keys()
		for key in kopt_keys:
			self.setOption(key,koptions[key])
		#end
		
		
	def __solve__(self, aero_problem, sol_type, *args, **kwargs):
		
		'''
		Run Analyzer (Analyzer Specific Routine)
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def __call__(self, aero_problem={}, sol_type='steady', *args, **kwargs):
		
		'''
		Run Analyzer (Calling Routine)
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Check Optimization Problem
		if not isinstance(aero_problem,AeroProblem):
			raise ValueError(" ")
		#end
		
		# Checks
		
		# Solve Analysis Problem
		self.__solve__(aero_problem, sol_type, *args, **kwargs)
		
		return 
		
		
	def _mesh(self, flow):
		
		'''
		Geometry Mesh
		'''
		
		# 
		
		
		
	def _solve(self, flow):
		
		'''
		Solve System of Equations
		
		Documentation last updated:  May. 26, 2008 - Ruben E. Perez
		'''
		
		# 
		
		
		
	def _forces(self, flow):
		
		'''
		Calculate Forces and Moments
		
		Documentation last updated:  May. 26, 2008 - Ruben E. Perez
		'''
		
		# 
		
		
		
	def _get_loads(self, flow):
		
		'''
		Get Aerodynamic Loads
		
		Documentation last updated:  May. 26, 2008 - Ruben E. Perez
		'''	
		
		#

	def initAdjoint(self, *args, **kwargs):
		'''
		Initialize the Ajoint problem for this test case
		'''
		pass

	def setupAdjointMatrix(self, *args, **kwargs):
		'''
		Setup the adjoint matrix for the current solution
		'''
		pass
		
	def solveAdjoint(self,objective, *args, **kwargs):
		'''
		Solve the adjoint problem for the desired objective functions.

		objectives - List of objective functions
		possibleObjectives - Generic objectives dictinoary
		'''
		possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
				       'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
				       'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
				       'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
				       'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
				       'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
				       'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
				       'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
				       'cMzAlpha':'cMzAlpha',\
				       'cM0':'cM0','cm0':'cM0',\
				       'clAlpha':'clAlpha',\
				       'cl0':'cl0',\
				       'cdAlpha':'cdAlpha',\
				       'cd0':'cd0'
				       }
		#for items in objectives:

		self._on_adjoint(objective,*args,**kwargs)
		#endfor
		
		
		
	def _on_adjoint(self, objective, *args, **kwargs):
		
		'''
		Adjoint
		
		Documentation last updated:  May. 26, 2008 - Ruben E. Perez
		'''	
		
		# 
		pass
		
	def computeSurfaceDerivative(self, objective, *args, **kwargs):
		'''
		Compute the derivative of the objective function wrt the surface.
		'''
		
		pass
		
	def computeDerivitiveCoupling(self, objective, *args, **kwargs):
		'''
		Compute the Aerodynmics portion of the augmented RHS in the 
			structural adjoint.
		This Boils down to dIdx_oml
		
		objective - objective function of interest.
		
		'''
		
		pass
		
	def getDerivativeCoupling(self, objective, *args, **kwargs):
		'''
		Retrieve the derivative vector computed and stored in 
			computeDerivative Coupling.
		
		objective - objective function of interest.
		'''
		
		pass
		
	def augmentAdjointRHS(self, aug_vector,objective, *args, **kwargs):
		'''
		Augment the RHS of the adjoint computation.
		
		aug_vector - the sensitivity vector from the other discipline with which to augment the RHS objective - objective function of interest.
		'''
		
		pass
		
	def _on_setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		raise NotImplementedError()
		
		
	def setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Calling Routine)
		
		Keyword arguments:
		-----------------
		name -> STRING: Option Name
		value -> SCALAR or BOOLEAN: Option Value
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		def_options = self.options['defaults']
		if def_options.has_key(name):
			if (type(value) == def_options[name][0]):
				self.options[name] = [type(value),value]
			else:
				raise IOError('Incorrect ' + repr(name) + ' value type')
			#end
		else:
			print '%s is not a valid option name'%(name)
			raise InputError('Not a valid option name')
		#end
		
		# 
	
		self._on_setOption(name, value)
		
		
	def _on_getOption(self, name):
		
		'''
		Get Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		raise NotImplementedError()
		
		
	def getOption(self, name):
		
		'''
		Get Optimizer Option Value (Calling Routine)
		
		Keyword arguments:
		-----------------
		name -> STRING: Option Name
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		def_options = self.options['defaults']
		if def_options.has_key(name):
			return self.options[name][1]
		else:	
			raise InputError(repr(name) + ' is not a valid option name')
		#end
		
		# 
		self._on_getOption(name)
		
		
	def _on_getInform(self, info):
		
		'''
		Get Optimizer Result Information (Optimizer Specific Routine)
		
		Keyword arguments:
		-----------------
		id -> STRING: Option Name
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		raise NotImplementedError()
		
		
	def getInform(self, infocode={}):
		
		'''
		Get Optimizer Result Information (Calling Routine)
		
		Keyword arguments:
		-----------------
		name -> STRING: Option Name
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		if (infocode == {}):
			return self.informs
		else:
			return self._on_getInform(infocode)
		#end
		
		
	def ListAttributes(self):
		
		'''
		Print Structured Attributes List
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		ListAttributes(self)
	


#==============================================================================
# 
#==============================================================================
def ListAttributes(self):
		
		'''
		Print Structured Attributes List
		
		Documentation last updated:  March. 24, 2008 - Ruben E. Perez
		'''
		
		print '\n'
		print 'Attributes List of: ' + repr(self.__dict__['name']) + ' - ' + self.__class__.__name__ + ' Instance\n'
		self_keys = self.__dict__.keys()
		self_keys.sort()
		for key in self_keys:
			if key != 'name':
				print str(key) + ' : ' + repr(self.__dict__[key])
			#end
		#end
		print '\n'
	


#==============================================================================
# Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	print 'Testing ...'
	
	# Test Optimizer
	azr = AeroSolver('Test')
	azr.ListAttributes()
	
