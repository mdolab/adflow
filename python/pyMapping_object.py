#!/usr/local/bin/python
'''
pyMapping_object - A Generic Python class for surface to surface mapping.

Copyright (c) 2008 by Mr.C.A (Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 10/07/2008 11:00$


Developers:
-----------
- Mr. C.A.(Sandy) Mader (SM)

History
-------
	v. 1.0	- Original pyAero Framework Implementation (SM 2008)
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
import time

# =============================================================================
# Extension modules
# =============================================================================


# =============================================================================
# Misc Definitions
# =============================================================================

class Mapping(object):

    '''
    Abstract Class for Surface to Surface mapping
    '''

    def __init__(self, name, category={}, def_options={}, informs={}, *args, **kwargs):
		
		'''
		Surface Class Initialization
		
		Documentation last updated:  July. 10, 2008 - C.A.(Sandy) Mader
                
		'''

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

                return


    def createMapping(self, *args, **kwargs):
        '''
        create a mapping between two surfaces
        '''

        pass


    def getMappedSurface(self, *args, **kwargs):
        '''
        Get the current surface
        '''

        pass

 
