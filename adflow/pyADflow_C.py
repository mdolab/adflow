#!/usr/bin/python
"""
pyADflow - A Python interface to ADflow.

Copyright (c) 2008 by Mr.C.A (Sandy) Mader
All rights reserved. Not to be used for commercial purposes.
Revision: 1.0   $Date: 03/07/2008 11:00$

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)
- Mr. C.A.(Sandy) Mader (SM)
- Dr. Ruben E. Perez (RP)

History
-------
v. 1.0  - Original pyAero Framework Implementation (RP,SM 2008)
"""

# =============================================================================
# Imports
# =============================================================================
import os
import time
import copy
import numpy
from mpi4py import MPI
from baseclasses import AeroSolver, AeroProblem
from . import MExt
from .pyADflow import ADFLOW
# =============================================================================
# ADFLOW_C Class
# =============================================================================
class ADFLOW_C(ADFLOW):

    def __init__(self, *args, **kwargs):
      
        # Load the compiled module using MExt, allowing multiple
        # imports
        debug = False
        if 'debug' in kwargs:
            debug=True
            
        curDir = os.path.dirname(os.path.realpath(__file__))
        self.adflow = MExt.MExt('libadflow_cs',[curDir], debug=debug)._module
        ADFLOW.__init__(self, dtype='D', *args, **kwargs)        

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
                for i in range(len(value)):
                    varStr = varStr + value[i] + '_'
                # end if
                varStr = varStr[0:-1] # Get rid of last '_'
                if name == 'monitorvariables':
                    self.adflow.monitorvariables(varStr)
                if name == 'surfacevariables':
                    self.adflow.surfacevariables(varStr)
                if name == 'volumevariables':
                    self.adflow.volumevariables(varStr)
            # end if
            if name == 'metricconversion':
                self.adflow.flowvarrefstate.lref = value
                self.adflow.flowvarrefstate.lrefspecified = True
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
      
        # Exec str is what is actually executed:
        exec_str = 'self.adflow.'+self.optionMap[name]['location'] + '=' + value
        exec(exec_str)
    
    def writeMeshFile(self, filename=None):
        if self.comm.rank == 0:
            print('Output not supported with complex version.')

    def writeVolumeSolutionFile(self, filename=None, writeGrid=True):
        if self.comm.rank == 0:
            print('Output not supported with complex version.')

    def writeSurfaceSolutionFile(self, *filename):
        if self.comm.rank == 0:
            print('Output not supported with complex version.')
