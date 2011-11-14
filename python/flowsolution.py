#! /usr/bin/env python
#
# File: flowsolution.py
# By: Patrick LeGresley
# Last Modified: 3/22/2004

__version__ = "0.1"

# Extension modules
from mesh import Mesh

class FlowSolution:
    """Represent a generic flow solution.

    This class does not actually implement any functionality but provides a
    template so that FlowSolution type classes can be defined for particular
    solvers.
    
    """

    def __init__(self):
        """Initialize the object.

        Each class deriving from FlowSolution will need to define its own set
        of appropriate keyword arguments.
        
        """
        self.Mesh = Mesh()
        return

    def GetSurfaceCoordinates(self):
        """Return the surface coordinates."""
        pass

    def GetSurfaceLoads(self,BoundaryCondition='Solid_Wall_Euler'):
        """Return the surface loads.

        Keyword arguments:
        BoundaryCondition -- the type of boundary condition
                             (default 'Solid_Wall_Euler')

        """
        pass

    def RunIterations(self,ncycles):
        """Run the flow solver for a given number of iterations.

        Keyword arguments:
        ncycles -- the number of iterations
        
        """
        pass

    def WriteCGNSFile(self,filename):
        """Write the current state of the flow solution to a CGNS file.

        Keyword arguments:
        filename -- the name of the file
        
        """
        pass
    
