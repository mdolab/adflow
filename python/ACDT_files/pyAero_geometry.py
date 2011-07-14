#!/usr/local/bin/python
'''
pyAero_geometry

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
# Geometry Class
# =============================================================================
class Geometry(object):
	
	'''
	Abstract Class for Geometry Object
	'''
	
	def __init__(self, name={},CGPercent = 0.25,ForeSparPercent = 0.25,
		     RearSparPercent = 0.75,StaticMarginPercent=0.05,
		     ForeThickCon = 0.01, RearThickCon = 0.99,
		     rootOffset = 0.01, tipOffset=0.01,
		     xRootec=0.0, yRootec=0.0, zRootec=0.0,
		     *args, **kwargs):
		
		'''
		Flow Class Initialization
		
		Keyword Arguments:
		------------------
		name -> STRING: Geometry Instance Name
		
		Attributes:
		-----------
		
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.name = name
		self.CGPercent = CGPercent
		self.ForeSparPercent = ForeSparPercent 
		self.RearSparPercent = RearSparPercent
		self.StaticMarginPercent = StaticMarginPercent
		self.ForeThickCon = ForeThickCon
		self.RearThickCon = RearThickCon
		self.tipOffset = tipOffset
		self.rootOffset = rootOffset
		self.xRootec = xRootec
		self.yRootec = yRootec
		self.zRootec = zRootec
		
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
		
		return ('name    \n'+'	 '+str(self.name).center(9) )
	


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
# Flow Test
#==============================================================================
if __name__ == '__main__':
	
	print 'Testing ...'
	
	# Test Variable
	geo = Geometry(name = 'test')
	geo.ListAttributes()
	print geo
