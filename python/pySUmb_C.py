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
    v. 1.0  - Original pyAero Framework Implementation (RP,SM 2008)
'''

__version__ = '$Revision$'

'''
To Do:
    - 
'''


# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, copy

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================

from mdo_import_helper import import_modules
exec(import_modules('pyAero_problem'))
exec(import_modules('pySUMB'))
import sumb_cs

# =============================================================================
# Misc Definitions
# =============================================================================

# =============================================================================
# SUMB Class
# =============================================================================
class SUMB_C(SUMB):
    
    '''
    SUmb Aerodynamic Analysis Class - Inherited from the SUMB Class
    '''
    
    def __init__(self, *args, **kwargs):
        
        '''
        SUMB Class Initialization
        
        Documentation last updated:  July. 03, 2008 - C.A.(Sandy) Mader
        '''

        SUMB.__init__(self,sumb=sumb_cs,*args,**kwargs)        
        self.dtype = 'D'
        return


    def _on_setOption(self, name, value):
        
        '''
        Set Optimizer Option Value (Optimizer Specific Routine)
        
        Documentation last updated:  May. 21, 2008 - Ruben E. Perez
        '''

        # Ignored options do NOT get set in solver

        if name in self.ignore_options:
            return

        # Do special Options individually
        if name in self.special_options:
            if name in ['monitorvariables','surfacevariables','volumevariables']:
                varStr = ''
                for i in xrange(len(value)):
                    varStr = varStr + value[i] + '_'
                # end if
                varStr = varStr[0:-1] # Get rid of last '_'
                if name == 'monitorvariables':
                    self.sumb.monitorvariables(varStr)
                if name == 'surfacevariables':
                    self.sumb.surfacevariables(varStr)
                if name == 'volumevariables':
                    self.sumb.volumevariables(varStr)
            # end if
            if name == 'metricconversion':
                self.sumb.flowvarrefstate.lref = value
                self.sumb.flowvarrefstate.lrefspecified = True
                self.metricConversion = value
            return
        # end if

        # All other options do genericaly by setting value in module:
        # Check if there is an additional mapping to what actually
        # has to be set in the solver

        temp = copy.copy(self.optionMap[name]) # This is the dictionary
        temp.pop('location')
        try:
            temp.pop('len')
        except:
            pass
        # end if

        # If temp has anything left in it, we MUST be able to match to
        # one of them.

        if len(temp) == 0:
            pass
        else:
            #Convert the value to lower case:
            value = self.optionMap[name][value.lower()]
        # end if

        # If value is a string, put quotes around it and make it
        # the correct length, otherwise convert to string
        if isinstance(value,str): 
            spacesToAdd = self.optionMap[name]['len'] - len(value)
            value = '\'' + value + ' '*spacesToAdd + '\''
        elif isinstance(value,float):
            # Force floats to be complex
            value = str(complex(value))
        else:
            value = str(value)
        # end if
      
        # Exec str is what is actually executed:
        exec_str = 'self.sumb.'+self.optionMap[name]['location'] + '=' + value
        exec(exec_str)

        return
        
    def writeMeshFile(self, filename=None):
        mpiPrint('Output not supported with complex version.',comm=self.comm)
        return

    def writeVolumeSolutionFile(self, filename=None, writeGrid=True):
        mpiPrint('Output not supported with complex version.',comm=self.comm)
        return

    def writeSurfaceSolutionFile(self, *filename):
        mpiPrint('Output not supported with complex version.',comm=self.comm)
        return



#==============================================================================
# SUmb Analysis Test
#==============================================================================
if __name__ == '__main__':
    
    # Test SUmb
    print 'Testing ...'
    sumb = SUMB_C()
    print sumb
    
