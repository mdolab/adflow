#!/usr/local/bin/python
'''
pyAero_flow

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
#import external

# =============================================================================
# Extension modules
# =============================================================================
#import extension

# =============================================================================
# Misc Definitions
# =============================================================================



# =============================================================================
# Flow Class
# =============================================================================
class Flow(object):
	
	'''
	Abstract Class for Flow Object
	'''
	def __init__(self, name={}, mach=0.0, alpha=0.0, beta=0.0,liftIndex=2,
		     phat=0.0, rhat=0.0, qhat=0.0, degreePol = 0,
		     coefPol = [0.0], degreeFourier = 0,
		     omegaFourier = 0.0, cosCoefFourier = [0.0], 
		     sinCoefFourier = [0.0],rho=1.225,P=101325.0,T=273.15,gamma=1.4,
		     mu=1.80e-5,nu=1.46e-5, *args, **kwargs):
		
		'''
		Flow Class Initialization
		
		Keyword Arguments:
		------------------
		name -> STRING: Flow Instance Name
		mach -> SCALAR: Mach Number
		alpha -> SCALAR: Angle of Attack
		beta -> SCALAR: Sideslip Angle
		phat -> SCALAR: Non-Dimensional Roll Rate
		rhat -> SCALAR: Non-Dimensional Yaw Rate
		qhat -> SCALAR: Non-Dimensional Pitch Rate
		
		Attributes:
		-----------
		
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.name = name
		if mach < 0:
			raise IOError('Mach number should be equal or greater than zero!')
		else:
			self.mach = mach
		#end

		self.alpha = alpha
		self.beta = beta
		self.liftIndex = liftIndex
		self.phat = phat
		self.rhat = rhat
		self.qhat = qhat
		self.degreePol = degreePol
		self.coefPol = coefPol
		self.degreeFourier = degreeFourier
		self.omegaFourier = omegaFourier
		
		self.cosCoefFourier = cosCoefFourier
		self.sinCoefFourier = sinCoefFourier
		
		self.rho = rho
		self.P   = P
		self.T = T
		self.gamma = gamma
		self.mu = mu
		self.nu = nu

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
# Flow Test
#==============================================================================
if __name__ == '__main__':
	
	print 'Testing ...'
	
	# Test Variable
	flow = Flow()
	flow.ListAttributes()
	
