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

# =============================================================================
# Misc Definitions
# =============================================================================



# =============================================================================
# AeroSolver Class
# =============================================================================
class Base(object):
	
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
		
		
	def __solve__(self, *args, **kwargs):
		
		'''
		Run Analyzer (Analyzer Specific Routine)
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def __call__(self,  *args, **kwargs):
		
		'''
		Run Analyzer (Calling Routine)
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Checks
		
		# Solve Analysis Problem
		self.__solve__(*args, **kwargs)
		
		return 

		
	def _solve(self, flow):
		
		'''
		Solve System of Equations
		
		Documentation last updated:  May. 26, 2008 - Ruben E. Perez
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
	bse = Base('Test')
	bse.ListAttributes()
	
