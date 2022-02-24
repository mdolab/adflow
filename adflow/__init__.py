__version__ = "2.7.1"

from mpi4py import MPI

from .pyADflow import ADFLOW
from .pyADflow_C import ADFLOW_C
from .oversetCheck import OversetCheck
from .checkZipper import checkZipper


try:
    import openmdao
except ImportError as err:
    openmdao = None

    if MPI.COMM_WORLD.rank == 0:
        print("Warning: OpenMDAO dependency is not installed. OM_ADFLOW wrapper will not be active")

if openmdao is not None:
    from .om_adflow import OM_ADFLOW
