from __future__ import print_function

from .python.pyADflow import ADFLOW
from .python.pyADflow_C import ADFLOW_C
from .python.oversetCheck import OversetCheck
from .python.checkZipper import checkZipper


try: 
    import openmdao
except ImportError as err:  
    openmdao = None
    print('Warning: OpenMDAO dependency is not installed. OM_ADFLOW wrapper will not be active')

if openmdao is not None: 
    from .python.om_adflow import OM_ADFLOW




