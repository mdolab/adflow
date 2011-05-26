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

		
	def _solve(self, flow):
		
		'''
		Solve System of Equations
		
		Documentation last updated:  May. 26, 2008 - Ruben E. Perez
		'''
		
		pass


	def resetFlow(self):
		'''
		Reset the flow to a uniform state
		'''
		
		pass


	def getSurfaceCoordinates(self,group_name):
		'''
		Return the set of surface coordinates cooresponding to a
		Particular group name
		'''
		
		pass


	def setSurfaceCoordinates(self,group_name,coordinates):
		'''
		Set the set of surface coordinates cooresponding to a
		Particular group name
		'''
		
		pass

	def getForces(self,group_name):
		'''
		Return the set of forces at the locations defined by 
		getSurfaceCoordinates
		'''
		
		pass


	def globalNKPreCon(self,in_vec):
		'''
		Precondition the residual in in_vec for a coupled 
		Newton-Krylov Method
		'''

		pass

	def totalSurfaceDerivative(self,objective):
		'''
		Return the total derivative of the objective at surface
		coordinates
		'''

		pass


	def totalAeroDerivative(self,objective):
		'''
		Return the total derivative of the objective with respect 
		to aerodynamic-only variables
		'''

		pass


	def getResNorms(self):
		'''
		Return the inital,starting and final residual norms for 
		the solver
		'''
		
		pass


	def getStateSize(self):
		'''
		Return the number of degrees of freedom (states) that are
		on this processor
		'''

		pass


	def getStates(self):
		'''
		Return the states on this processor.
		'''

		pass


	def setStates(self,states):
		''' Set the states on this processor.'''

		pass

	def getResidual(self):
		'''
		Return the reisudals on this processor.
		'''

		pass


	def getSolution(self):
		'''
		Retrieve the solution dictionary from the solver
		'''
		
		pass


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
	
		'''
	
		self._on_adjoint(objective,*args,**kwargs)
			
		
	def _on_adjoint(self, objective, *args, **kwargs):
		
		'''
		Adjoint
		
		Documentation last updated:  May. 26, 2008 - Ruben E. Perez
		'''	
		
		# 
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
	
