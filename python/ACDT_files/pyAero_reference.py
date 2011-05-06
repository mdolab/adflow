#!/usr/local/bin/python
'''
pyAero_reference

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
To Do:
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
import numpy

# =============================================================================
# Extension modules
# =============================================================================


# =============================================================================
# Misc Definitions
# =============================================================================



# =============================================================================
# Reference Class
# =============================================================================
class Reference(object):
	
	'''
	Abstract Class for Aerodynamic Reference Object
	'''
	
	def __init__(self, name={}, sref=None, bref=None, cref=None, xref=0.0, yref=0.0, zref=0.0, *args, **kwargs):
		
		'''
		Reference Class Initialization
		
		Keyword Arguments:
		------------------
		sref -> SCALAR: Reference Area
		bref -> SCALAR: Reference Span
		cref -> SCALAR: Reference Chord
		xref -> SCALAR: Reference Position x-coordinate
		yref -> SCALAR: Reference Position y-coordinate
		zref -> SCALAR: Reference Position z-coordinate
		
		Attributes:
		-----------
		
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.name = name
		self.sref = sref
		self.bref = bref
		self.cref = cref
		self.xref = xref
		self.yref = yref
		self.zref = zref
	
		
	def ListAttributes(self):
		
		'''
		Print Structured Attributes List
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		ListAttributes(self)
		
		
	def __str__(self):
		
		'''
		Print Structured List of Variable
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		return ('name    mach	alpha	beta	phat	rhat	qhat\n'+'	 '+str(self.name).center(9) +'%8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n' %(self.mach, self.alpha, self.beta, self.phat, self.rhat, self.qhat))
	


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
# Reference Test
#==============================================================================
if __name__ == '__main__':
	
	print 'Testing ...'
	
	# Test Reference
	ref = Reference()
	ref.ListAttributes()
	
