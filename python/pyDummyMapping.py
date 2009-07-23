#!/usr/local/bin/python
'''
pyDummyMapping - A dummy mapping that simply passes coords through

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
from pyMapping_object import Mapping

# =============================================================================
# Misc Definitions
# =============================================================================

class DUMMYMAPPING(Mapping):

    '''
    Abstract Class for Surface to Surface mapping
    '''

    def __init__(self, *args, **kwargs):
		
		'''
		DUMMYMAPPING Class Initialization
		
		Documentation last updated:  July. 21, 2008 - C.A.(Sandy) Mader
                
		'''
                name = 'DUMMYMAPPING'
		category = 'Surface to surface mapping tool'
		def_opts = {
			}
		informs = {
			}

                #run the Generic surface defintion
		Mapping.__init__(self, name, category, def_opts, informs, *args, **kwargs)

                return


    def createMapping(self,xyzmap,xyzref,conn,elemtype, *args, **kwargs):
        '''
        create a mapping between two surfaces
        '''
        #No need to do anything...
        

        return


    def getMappedSurface(self,xyznew, *args, **kwargs):
        '''
        Get the current surface
        '''
        xyz_mapped = xyznew

        return xyz_mapped

 

