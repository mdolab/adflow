#! /usr/bin/env python
#
# File: mesh.py
# By: Patrick LeGresley
# Last Modified: 3/22/2004

__version__ = "0.1"

class Mesh:
    """Represent a generic mesh.

    This class does not actually implement any functionality but provides a
    template so that more specific mesh type classes can be defined.
    
    """

    def __init__(self):
        """Initialize the object.

        Each class deriving from Mesh will need to define its own set
        of appropriate keyword arguments.
        
        """
        pass

    def GetSurfaceCoordinates(self):
        """Return the surface coordinates."""
        pass

    def WriteCGNSFile(self,filename):
        """Write the current state of the flow solution to a CGNS file.

        Keyword arguments:
        filename -- the name of the file
        
        """
        pass

    def SetCoordinates(self):
        """Set the coordinates for one or more portions of a mesh."""
        pass

# Import the multiblock meshwarping module
###import some_module
    
class MultiblockMesh(Mesh):
    """Represent a Multiblock Mesh."""

    def __init__(self,filename):
	"""Initialize the object.

        Keyword arguments:

        filename -- the name of the CGNS file to read.

        """
        if (os.path.isfile(filename)):	
            ###some_module.read_grid(filename)
            self.filename = filename
        else:
            print "Error: Could not find file %s" % filename
            self.filename = None
        self.master_lists = False
        return

    def GetCoordinates(self):
        """Return coordinates ...

        """
        print "Error: This method is not implemented yet."
        return

    def Warp(n,blnum,ijk,xyz):
        """Warp the mesh ...

        """
        if not(self.master_lists):
            ###some_module.create_master_corner_list()
            ###some_module.create_master_edge_list()
            pass ###

        ###some_module.warp(self.filename,n,blnum,ijk,xyz)

        return

    def ReadCGNSFile(self,filename):
         """Replace the existing mesh by reading a new CGNS file.

         Keyword arguments:

         filename -- the name of the CGNS file to read.

         """
         if (os.path.isfile(filename)):
            ###some_module.read_grid(filename)
            self.filename = filename
         else:
            print "Error: Could not find file %s" % filename
         self.master_lists = False
         return

    def WriteCGNSFile(self,*filename):
         """Write the current mesh to a CGNS file.

         Keyword arguments:

         filename -- the name of the file to write, if not provided the name
                     of the file originally read in is used.

         """
         print "Error: This method is not yet implemented."
         return

