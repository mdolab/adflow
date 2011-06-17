#!/usr/local/bin/python
'''
pyAero_problem

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
from pyAero_geometry import Geometry
from pyAero_flow import Flow

# =============================================================================
# Misc Definitions
# =============================================================================



# =============================================================================
# AeroProblem Class
# =============================================================================
class AeroProblem(object):
	
	'''
	Aerodynamic Problem Class
	'''
	
	def __init__(self, name, geom={}, flow_set={}, ref_set={}, *args, **kwargs):
		
		'''
		AeroProblem Class Initialization
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.name = name
		
		# Initialize Geometry
		self._geometry = geom
		
		# Initialize Flow Set
		self._flows = flow_set
		
		# Initialize Reference Set
		self._refs = ref_set
		
		
	def getFlow(self, i):
		
		'''
		Get Flow Instance from Flow Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Check Index
		if not (isinstance(i,int) and i >= 0):
			raise ValueError("Flow index must be an integer >= 0.")
		#end
		
		# 
		return self._flows[i]
		
		
	def addFlow(self, *args, **kwargs):
		
		'''
		Add Flow Instance into Flow Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.setFlow(self.firstavailableindex(self._flows),*args,**kwargs)
		
		
	def setFlow(self, i, *args, **kwargs):
		
		'''
		Set Flow Instance into Flow Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		if (len(args) > 0) and isinstance(args[0], Flow):
			self._flows[i] = args[0]
		else:
			try:
				self._flows[i] = Flow(*args,**kwargs)
			except:
				raise ValueError("Input is not a Valid for a Flow Object instance\n")
			#end
		#end
		
		
	def delFlow(self):
		
		'''
		Delete Flow Instance from Flow Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Check Index
		if not (isinstance(i,int) and i >= 0):
			raise ValueError("Flow index must be an integer >= 0.")
		#end
		
		# 
		del self._flows[i]
		
		
	def getFlowSet(self):
		
		'''
		Get Flows Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		return self._flows
		
		
	def getRef(self, i):
		
		'''
		Get Reference Point Instance from Reference Points Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Check Index
		if not (isinstance(i,int) and i >= 0):
			raise ValueError("Flow index must be an integer >= 0.")
		#end
		
		# 
		return self._refs[i]
		
		
	def addRef(self, *args, **kwargs):
		
		'''
		Add Reference Point Instance into Reference Points Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.setRef(self.firstavailableindex(self._refs),*args,**kwargs)
		
		
	def setRef(self, i, *args, **kwargs):
		
		'''
		Set Reference Point Instance into Reference Points Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		if (len(args) > 0) and isinstance(args[0], Reference):
			self._refs[i] = args[0]
		else:
			try:
				self._refs[i] = Reference(*args,**kwargs)
			except:
				raise ValueError("Input is not a Valid for a Flow Object instance\n")
			#end
		#end
		
		
	def delRef(self):
		
		'''
		Delete Flow Instance from Reference Points Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Check Index
		if not (isinstance(i,int) and i >= 0):
			raise ValueError("Reference index must be an integer >= 0.")
		#end
		
		# 
		del self._refs[i]
		
		
	def getRefSet(self):
		
		'''
		Get Reference Points Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		return self._refs
		
		
	def getSol(self, i):
		
		'''
		Get Solution from Solution Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Check Index
		if not (isinstance(i,int) and i >= 0):
			raise ValueError("Solution index must be an integer >= 0.")
		#end
		
		# 
		return self._solutions[i]
		
		
	def addSol(self, *args, **kwargs):
		
		'''
		Add Solution into Solution Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		self.setSol(self.firstavailableindex(self._solutions),*args,**kwargs)
		
		
	def setSol(self,i, *args, **kwargs):
		
		'''
		Set Solution into Solution Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		if (len(args) > 0) and isinstance(args[0], Solution):
			self._solutions[i] = args[0]
		else:
			#try:
			self._solutions[i] = Solution(*args,**kwargs)
			#except:
			#	print args
			#	print kwargs
			#	raise ValueError("Input is not a Valid for a Solution Object instance\n")
			#end
		#end
		
		
	def delSol(self):
		
		'''
		Delete Solution from Solutions Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# Check Index
		if not (isinstance(i,int) and i >= 0):
			raise ValueError("Solution index must be an integer >= 0.")
		#end
		
		# 
		del self._solutions[i]
		
		
	def getSolSet(self):
		
		'''
		Get Solutions Set
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		return self._solutions
		
		
	def firstavailableindex(self, set):
		
		'''
		List First Unused Index from Variable Objects List
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		i = 0
		while i in set: 
			i += 1
		#end
		
		return i
		
		
	def ListAttributes(self):
		
		'''
		Print Structured Attributes List
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		ListAttributes(self)
		
		
	def __str__(self):
		
		'''
		Print Structured Optimization Problem
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		text = """\Aerodynamic Problem -- %s\n%s\n"""
		
		return (text)
		
		
	def write2file(self, outfile='', disp_sols=False, **kwargs):
		
		'''
		Write Structured Optimization Problem to file
		
		Keyword arguments:
		-----------------
		outfile -> STRING/INST: File name or file instance.
		disp_sols -> BOOL: Display solutions falg  Default=False.
		solutions -> LIST: List of solution indixes.
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		if isinstance(outfile,str):
			if (outfile == ''):
				findir = os.listdir(os.curdir)
				tmpname = self.name.lower()
				tmpname = tmpname.split(' ')
				tmpname = tmpname[0]
				i = 0
				while (tmpname+'.txt') in findir:
					tmpname = tmpname.rstrip('_%d' %(i-1))
					tmpname = tmpname + '_' +str(i)
					i += 1
				#end
				tmpname += '.txt'
				outfile = open(tmpname,'w')
			else:
				outfile = open(outfile,'w')
		elif (not isinstance(outfile,str)) and (not isinstance(outfile,file)):
			raise IOError(repr(outfile) + 'is not a file or filename')
		#end
		ftext = self.__str__()
		outfile.write(ftext)
		if disp_sols or kwargs.has_key('solutions'):
			if kwargs.has_key('solutions'):
				sol_indices = kwargs['solutions']
			else:
				sol_indices = self._solutions.keys()
			#end
			for key in sol_indices:
				soltext = '\n' + self._solutions[key].__str__() 
				outfile.write(soltext)
			#end
		#end
		print 'Data written to file ', outfile.name
		outfile.close()
	


# =============================================================================
# Solution Class
# =============================================================================
class Solution(AeroProblem):
	
	'''
	Aerodynamic Analysis Solution Class
	'''
	
	def __init__(self, analyzer, name, sol_time, sol_inform, geom={}, flow_set={}, ref_set={}, options_set={}, *args, **kwargs):
		
		'''
		Solution Class Initialization
		
		Keyword arguments:
		-----------------
		
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		# 
		AeroProblem.__init__(self, name, geom, flow_set, ref_set, *args, **kwargs)
		self.solver = solver
		self.sol_time = sol_time
		self.sol_inform = sol_inform
		self.options_set = options_set
		
		if kwargs.has_key('display_opts'):
			self.display_opt = kwargs['display_opts']
			del kwargs['display_opts']
		else:
			self.display_opt = False
		#end
		self.parameters = kwargs
		
		
	def __str__(self):
		
		'''
		Print Structured Solution
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		text0 = Analysis.__str__(self)
		text1 = ''
		lines = text0.split('\n')
		lines[1] = lines[1].lstrip('Aerodynamic Analysis --')
		for i in xrange(5):
			text1 += lines[i] + '\n'
		#end
		if self.display_opt:
			text1 += '\n	Options:\n '
			opt_keys = self.options_set.keys()
			opt_keys.sort()
			for key in opt_keys:
				ns = 25-len(key)
				text1 += '		'+ key +':' + str(self.options_set[key][1]).rjust(ns,'.') + '\n'
			#end
		#end
		text1 += '\n    Solution: \n'
		text1 += ('-'*80) + '\n'
		text1 += '    Total Time: %25.4f\n' %(self.sol_time)
		for key in self.parameters.keys():
			if (isinstance(self.parameters[key],(dict,list,tuple))) and (len(self.parameters[key]) == 0):
				continue
			elif (isinstance(self.parameters[key],numpy.ndarray)) and (0 in (self.parameters[key]).shape):
				continue
			else:
				text1 += '    '+ key +': ' + str(self.parameters[key]).rjust(9) + '\n'
			#end
		#end
		for i in xrange(5,len(lines)):
			text1 += lines[i] + '\n'
		#end
		text1 += ('-'*80) + '\n'
		
		return text1
		
		
	def write2file(self, outfile):
		
		'''
		Write Structured Solution to file
		
		Documentation last updated:  May. 21, 2008 - Ruben E. Perez
		'''
		
		Analysis.write2file(self,outfile,False)
	


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
# Analysis Test
#==============================================================================
if __name__ == '__main__':
	
	print 'Testing ...'
	
	# Test Analysis
	
