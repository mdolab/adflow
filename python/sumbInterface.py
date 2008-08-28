#! /usr/bin/env python
#
# File: sumbInterface.py
# Authors: Patrick LeGresley, Andre Marta, Steve Repsher
# Modified by: C.A.(Sandy) Mader
# Last Modified: 03/07/2008

__version__ = "1.0"

# =============================================================================
# Standard Python modules
# =============================================================================
import sys
import os
import copy
import string
import types
from math import pi

# =============================================================================
# Extension modules
# =============================================================================

import numpy

# Try to import the mpi module
try:
    from mpi4py import MPI as mpi
#    import mpi
    _parallel = True
    print 'importingMPI',_parallel
except ImportError:
    try:
        import dummy_mpi as mpi
        _parallel = False
    except ImportError:
        print "Error: Failed to import mpi or dummy_mpi."

# Import the sumb module
if _parallel:
    try:
        import sumb_parallel as sumb
    except ImportError:
        try:
            import sumb
            print "Warning: Running in an MPI environment, but failed"
            print "         to import parallel version of SUmb.  Proceeding"
            print "         with a sequential version."
        except ImportError:
            print "Error: Failed to import parallel or sequential version"
            print "       of SUmb."
else:
    try:
        import sumb
    except ImportError:
        print "Error: Failed to import sequential version of SUmb."


# =============================================================================

class SUmbMesh(object):
    """Represents a SUmb mesh.
 
    """
 
    def __init__(self):
        """Initialize the object."""
        self._update_geom_info = False
        self._suggar_interface_initialized = False

    def GetSurfaceCoordinates(self, family=None, sps=1):
        """Return surface coordinates.
        
        Keyword arguments:
        
        family -- optional string specifying the return of surface coordinates
                  only for the specified family.
        sps    -- spectral time step (optional, default is set to 1).
                  (sps=1 for usual steady or unsteady models)

        """
        
        if(family):
            try:
                index = self.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdcreatesurfcoorlist(sps,index)
            return sumb.mddata.mdsurfxx[:,start_ind-1:end_ind]
        else:
            nfamilies = len(self.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdcreatesurfcoorlist(sps,0) 
            else:
                for n in range(nfamilies): 
                    [start_ind,end_ind] = sumb.mdcreatesurfcoorlist(sps,n+1)
            return sumb.mddata.mdsurfxx

    def GetSurfaceCoordinatesLocal(self, family=None, sps=1):
        """Return surface coordinates for this processor only.
        
        Keyword arguments:
        
        family -- optional string specifying the return of surface coordinates
                  only for the specified family.
        sps    -- spectral time step (optional, default is set to 1).
                  (sps=1 for usual steady or unsteady models)

        """
        
        if(family):
            try:
                index = self.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdcreatesurfcoorlistlocal(sps,index)
            return sumb.mddatalocal.mdsurfxxlocal[:,start_ind-1:end_ind]
        else:
            nfamilies = len(self.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdcreatesurfcoorlistlocal(sps,0) 
            else:
                for n in range(nfamilies): 
                    [start_ind,end_ind] = sumb.mdcreatesurfcoorlistlocal(sps,n+1)
            return sumb.mddatalocal.mdsurfxxlocal

    def DeallocateSurfaceCoordinates(self):
        """Deallocate memory used for the surface coordinates."""
        sumb.mddeletesurfcoorlist()

          
    def DeallocateSurfaceCoordinatesLocal(self):
        """Deallocate memory used for the surface coordinates."""
        sumb.mddeletesurfcoorlistlocal()

    def GetSurfacePatchDimensions(self, family=None):
        """Return surface patch dimensions.

        Keyword arguments:

        family -- optional string specifying the return of surface patch
                  dimensions only for the specified family.

        """
        if(family):
            try:
                index = self.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdsurfacepatchdim(index)
            return sumb.mddata.mdpatchdimensions[:,start_ind-1:end_ind]
        else:
            nfamilies = len(self.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdsurfacepatchdim(0)
            else:
                for n in range(nfamilies):
                    [start_ind,end_ind] = sumb.mdsurfacepatchdim(n+1)
            return sumb.mddata.mdpatchdimensions

    def DeallocateSurfacePatchDimensions(self):
        """Deallocate memory used for the surface patch dimensions."""
        sumb.mddeletesurfacepatchdim()
                                             
    def GetSurfaceIndices(self, family=None):
        """Return CGNS block IDs and indices.  This is a 2-D array, dimensions
        (4,npoints).  The first 3 indices are the i,j,k and the fourth is the
        CGNS block ID.
                                             
        Keyword arguments:
                                             
        family -- optional string specifying the return of indices only for the
                  specified family.
                                             
        """
        if(family):
            try:
                index = self.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdcreatesurfindlist(index)
            return sumb.mddata.mdsurfind[:,start_ind-1:end_ind]
        else:
            nfamilies = len(self.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdcreatesurfindlist(0)
            else:
                for n in range(nfamilies):
                    [start_ind,end_ind] = sumb.mdcreatesurfindlist(n+1)
            return sumb.mddata.mdsurfind

    def GetSurfaceIndicesLocal(self, family=None):
        """Return CGNS block IDs and indices.  This is a 2-D array, dimensions
        (4,npoints).  The first 3 indices are the i,j,k and the fourth is the
        CGNS block ID.
                                             
        Keyword arguments:
                                             
        family -- optional string specifying the return of indices only for the
                  specified family.
                                             
        """
        if(family):
            try:
                index = self.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdcreatesurfindlistlocal(index)
            return sumb.mddata.mdsurfindlocal[:,start_ind-1:end_ind]
        else:
            nfamilies = len(self.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdcreatesurfindlistlocal(0)
            else:
                for n in range(nfamilies):
                    [start_ind,end_ind] = sumb.mdcreatesurfindlistlocal(n+1)
            return sumb.mddatalocal.mdsurfindlocal
                                             
    def DeallocateSurfaceIndices(self):
        """Deallocate memory used for the surface indices."""
        sumb.mddeletesurfindlist()

    def DeallocateSurfaceIndicesLocal(self):
        """Deallocate memory used for the surface indices."""
        sumb.mddeletesurfindlistlocal()

    def WriteMeshFile(self,*filename):
        """Write the current state of the mesh to a CGNS file.
 
        Keyword arguments:
        
        filename -- the name of the file (optional)
         
        """
        if (filename):
            sumb.inputio.newgridfile[:] = ''
            sumb.inputio.newgridfile[0:len(filename[0])] = filename[0]
        sumb.monitor.writegrid=True
        sumb.monitor.writevolume=False
        sumb.monitor.writesurface=False
        sumb.writesol()

    def DummySetCoordinates(self, sps =1):
        """Dummy SetCoordinates routine to be called on processors where no
        new coordinates reside.
 
        Keyword arguments:

        sps    -- spectral time step (optional, default set to 1).
                        (sps= 1 for usual steady or unsteady models)

        """
        sumb.iteration.groundlevel = 1
        blocknums = numpy.zeros((0))
        ranges = numpy.zeros((3,2,0))
        xyz = numpy.zeros((3,0),'d')
        sumb.mdsetcoor(sps,blocknums,ranges,xyz)
        self._update_geom_info = True
 
    def SetCoordinates(self, blocknums, ranges, xyz, sps=1):
        """Set the coordinates for a given block or sub-blocks.
                                                                                
        Keyword arguments:
                                                                                
        blocknums -- the block number for each sub-block, dimension (nsubblocks)
        ranges -- the ijk ranges for each sub-block, dimension (3,2,nsubblocks)
        xyz -- the new xyz coordinates, dimension (3,npts)
        sps    -- spectral time step (optional, default is set to 1).
                  (sps=1 for usual steady or unsteady models)
                                                                                
        """
        sumb.iteration.groundlevel = 1
        sumb.mdsetcoor(sps,blocknums,ranges,xyz)
        self._update_geom_info = True

    def SetCoordinatesLocal(self, blocknum,il,jl,kl,xyz, sps=1):
        """Set the coordinates for a given block or sub-blocks.
                                                                                
        Keyword arguments:
                                                                                
        blocknums -- the block number for each sub-block, dimension (nsubblocks)
        ranges -- the ijk ranges for each sub-block, dimension (3,2,nsubblocks)
        xyz -- the new xyz coordinates, dimension (3,npts)
        sps    -- spectral time step (optional, default is set to 1).
                  (sps=1 for usual steady or unsteady models)
                                                                                
        """
        sumb.iteration.groundlevel = 1
        #sumb.mdsetcoor(sps,blocknums,ranges,xyz.real)#only the real part needs to be set in SUmb
        sumb.setblockcoords(blocknum,il,jl,kl,xyz.real)
        #sumb.mdsetcoor(sps,blocknums,ranges,xyz)
        self._update_geom_info = True

    def _UpdateGeometryInfo(self):
        """Update the SUmb internal geometry info, if necessary."""
        if (self._update_geom_info):
            sumb.updatecoordinatesalllevels()        
            sumb.updatewalldistancealllevels()
            sumb.updateslidingalllevels()
            sumb.updatemetricsalllevels()
            self._update_geom_info = False

    def GetNumberBlocks(self):
        """Get the number of blocks in the mesh."""
        return self.nmeshblocks

    def GetFamilyNames(self):
        """Return a list of the family names contained in the mesh."""
        return self.families.keys()
 
    def GetFilename(self):
        """Return the name of the CGNS file from which this mesh was orginally
        read.

        """
        return string.strip(sumb.inputio.gridfile.tostring())

    def _InitializeSuggarInterface(self):
        """Sets up ability to call SUGGAR interface methods."""

        self._suggar_interface_initialized = True
        sumb.initsuggarinterface()

        # Create a dictionary of zone names.
        self.zones = {}
        for i in range(sumb.suggardata.nzones[0]):
            self.zones[string.strip(sumb.suggardata.zonenames[i]
                       .tostring())] = sumb.suggardata.unsortedzone[i]

    def LoadSuggarDCIFiles(self, dcifilelist, sps=1):
        """Loads the SUGGAR ouput DCI files in the list to SUmb. The previous
           overset data for all zones affected by the files is cleared.

        Keyword arguments:

        dcifilelist -- list of SUGGAR output DCI filenames to load.
        sps         -- spectral mode (optional, default is set to 1).

        """
        # Initialize the interface with SUGGAR if not already done.
        if not self._suggar_interface_initialized:
            self._InitializeSuggarInterface()

        # Check to make sure each file exists and load it.
        for dcifile in dcifilelist:
            if os.path.isfile(dcifile):
                sumb.loadsuggardcifile(dcifile, sps)
            else:
                print 'Error: Could not find DCI file %s' % dcifile
                return None

    def WritePlot3DZoneFiles(self, zonelist=None, dirpath=None, sps=1, byteswap=False):
        """Writes a Plot3D grid file for each zone in zonelist to given path.

        Keyword arguments:

        zonelist -- list of zone names for which to write files.
                    (optional, default is all zones in mesh)
        dirpath  -- path to the directory in which files are to be written
                    (optional, default is current working directory)
        sps      -- spectral mode (optional, default is set to 1).
        byteswap -- Whether or not to apply byte-swapping to the output
                    files (optional, default is False so that the default 
                    byte-ordering or endian of the machine is used).

        """
        # Initialize the interface with SUGGAR if not already done.
        if not self._suggar_interface_initialized:
            self._InitializeSuggarInterface()

        # If no zone list was given then set it to the list of all zones.
        if zonelist is None:
            zonelist = self.zones.keys()

        # If no directory was given then set it to the current one.
        if dirpath is None:
            dirpath = "."

        # Write the file to dirpath/zonename.x for each zone in the list.
        for zonename in zonelist:
            try:
                index = self.zones[zonename]
            except KeyError:
                print "Error: No such zone '%s'" % zonename
                return None

            filename = dirpath + "/" + zonename + ".x"
            sumb.writeplot3dzonefile(index, filename, sps, byteswap)

        # Wait for all zone files to be written before returning since the
        # next step will be to run SUGGAR.
        self.sumb_comm_world.barrier

    def getBlockDimensions(self,blocknum):
        """ Get the i,j,k dimensions of block blocknum"""
        return  sumb.getblockdims(blocknum)

    def getBlockCoordinates(self,blocknum,il,jl,kl):
        """get the xyz coordinates from blockblocknum"""
        return sumb.getblockcoords(blocknum,il,jl,kl)
    def getBlockCGNSID(self,blocknum):
        """ Get original CGNS blockID for block blocknum"""
        return sumb.getblockcgnsid(blocknum)
    def getNSubfacesBlock(self,blocknum):
        """ get subface infor for block: blocknum"""
        return sumb.getnsubfacesblock(blocknum)
    def getBlockCommunicationInfo(self,blocknum,nSubface,n1to1,nNonMatch):
        """Get all fo the relevant MPI communication ifo for this block"""
        return sumb.getblockcommunicationinfo(blocknum,nSubface,n1to1,nNonMatch)
    def getNBlocksLocal(self):
        """ get the number of blocks present on the local processor"""
        return sumb.getnblockslocal()


# =============================================================================

class SUmbInterface(object):
    """Represents a SUmb flow solution."""

    def __init__(self,communicator=None,deforming_mesh=False):
        """Initialize the object.

        Keyword arguments:

        communicator   -- an MPI communicator for specifying which processors
                          to run on (optional)
        deforming_mesh -- set to True or False whether or not this is a 
                          deforming mesh case (optional, default is False)
                          

        """
        #Initialize the MPI in the flow solver
        print 'running initialization'
        sumb.sumb_init()

        # Setup a mesh object
        self.Mesh = SUmbMesh()

        # Create a Python communicator to mirror the Fortran
        # SUmb_COMM_WORLD and set the Fortran SUmb communicator to the
        # Python one.
        if (communicator is not None):
            self.sumb_comm_world = communicator
        else:
#            self.sumb_comm_world = mpi.COMM_WORLD.comm_create(
#        				     mpi.COMM_WORLD[:])
            self.sumb_comm_world = mpi.COMM_WORLD
        #sumb.communication.sumb_comm_world = int(self.sumb_comm_world)
        self.Mesh.sumb_comm_world = self.sumb_comm_world

##         # Store the name of the input file
##         self.startfile = startfile
##         sumb.inputio.paramfile[0:len(startfile)] = startfile

        # Determine the rank and number of processors inside the group
        # defined by sumb_comm_world.
        self.myid = sumb.communication.myid = self.sumb_comm_world.rank
        self.nproc = sumb.communication.nproc = self.sumb_comm_world.size

        # Allocate the memory for SENDREQUESTS and RECVREQUESTS.
        try:
            sumb.communication.sendrequests = numpy.zeros(
        			  (self.sumb_comm_world.size))
            sumb.communication.recvrequests = numpy.zeros(
        			  (self.sumb_comm_world.size))
        except:
            print "Memory allocation failure for SENDREQUESTS " \
        	  "and RECVREQUESTS."
            return

        # Set the SUmb module value of standalonemode to false and
        # the value of deforming_grid to the input value.
        sumb.iteration.standalonemode = False
        sumb.iteration.deforming_grid = deforming_mesh

        # Write the intro message
        sumb.writeintromessage()


        return

    def initializeFlow(self,aero_problem,sol_type,grid_file,*args,**kwargs):
        ''' Take in flow object,reference object, solutionType and
            desired meshfile and reinitialize the flow solution in SUmb

            Keyword arguments

            flow - flow condition to be solved
            ref  - reference condition
            sol_type - solution type - Steady,Unsteady,Time Spectral
            grid_file - name of 3-d Mesh file
            '''
        #Generate Input File from Options
        self.generateInputFile(aero_problem,sol_type,grid_file,*args,**kwargs)
        
        startfile = 'autogen.input'
        #startfile      -- the name of the SUmb input parameter file

        # Make sure the parameter file exists
        if not os.path.isfile(startfile):
            print 'Error: Could not find file %s' % startfile
            return None

        # Store the name of the input file
        self.startfile = startfile
        
        sumb.inputio.paramfile[0:len(startfile)] = startfile
        
        # Read the parameter file
        sumb.readparamfile()

        # Partition the blocks and read the grid
        sumb.partitionandreadgrid()
        
        # Perform the preprocessing task
        sumb.preprocessing()
        
        # Initialize the flow variables
        sumb.initflow()

        # Create dictionary of variables we are monitoring
        nmon = sumb.monitor.nmon
        self.monnames = {}
        for i in range(nmon):
            self.monnames[string.strip(
        		   sumb.monitor.monnames[i].tostring())] = i
        
        # Create dictionary of the family names
        sumb.mdgetfamilynames()
        nfamilies = sumb.mddata.mdnfamilies
        self.Mesh.families = {}
        for i in range(nfamilies):
            self.Mesh.families[string.strip(sumb.mddata.mdfamilynames[i]
        		       .tostring())] = i + 1
        sumb.mdcreatensurfnodes()

        # Determine the total number of blocks in the mesh and store it
        self.Mesh.nmeshblocks = self.sumb_comm_world.Allreduce(
        			     sumb.block.ndom,mpi.SUM)

        #Set flags for ADjoint initialization
        self.adjointInitialized = False
        
        return

    def setInflowAngle(self,aero_problem):
        '''
        Set the alpha and beta fromthe desiggn variables
        '''
        [velDir,liftDir,dragDir]= sumb.adjustinflowangleadj((aero_problem._flows.alpha*(pi/180.0)),(aero_problem._flows.beta*(pi/180.0)),aero_problem._flows.liftIndex)
        sumb.inputphysics.veldirfreestream = velDir
        sumb.inputphysics.liftdirection = liftDir
        sumb.inputphysics.dragdirection = dragDir

        return
        
    def generateInputFile(self,aero_problem,sol_type,grid_file,file_type='cgns',eqn_type='Euler',*args,**kwargs):
        ''' Code to generate an SUmb Input File on the fly'''
        print 'generating input file'
        print 'flow',aero_problem._flows

        #Convert alpha and beta to a freestream vector
        [velDir,liftDir,dragDir]= sumb.adjustinflowangleadj((aero_problem._flows.alpha*(pi/180.0)),(aero_problem._flows.beta*(pi/180.0)),aero_problem._flows.liftIndex)
        
        autofile = open("autogen.input",'w')

        # Write the header of the file

        autofile.write("===============================================================================\n")
        autofile.write("This is an automatically created parameter file for SUmb.\n")
        autofile.write("It contains all the parameters as specified in aero_problem.\n\n")
        autofile.write("Comments are indicated by #. All info following a # is ignored.\n")
        autofile.write("The sequence of the parameters is arbitrary and keywords are case insensitive.\n")
        autofile.write("When a keyword occurs multiple times,the last value is taken.\n")
        autofile.write("===============================================================================\n\n")

        # Write the keywords and default values for the IO parameters.

        autofile.write("-------------------------------------------------------------------------------\n")
        autofile.write("     IO Parameters\n")
        autofile.write("-------------------------------------------------------------------------------\n")
        if file_type =='cgns':
            autofile.write("                  File format Read: CGNS\n")
            autofile.write("                 File format write: CGNS\n\n")
        elif file_type == 'plot3d':
            autofile.write("                  File format Read: PLOT3D\n")
            autofile.write("                 File format write: PLOT3D\n\n")
        else:
            print 'invalid file type'
            #raiseError
        #endif

        if file_type =='cgns':
            autofile.write("                         Grid file: %s.cgns \n"%(grid_file))
            autofile.write("\n")
            autofile.write("                      Restart file: restart.cgns\n")
            autofile.write("                           Restart: no\n")
            autofile.write("       Check nondimensionalization: yes\n\n")
            autofile.write("                     New grid file: %s_NewGrid.cgns\n\n"%(grid_file))
            autofile.write("                     Solution file: %s_SolSUmb.cgns\n"%(grid_file))
            autofile.write("             Surface solution file: %s_SolSUmb_surface.cgns\n"%(grid_file))
            
        elif file_type == 'plot3d':
            autofile.write("                         Grid file: %s.mesh \n"%(grid_file))
            autofile.write("          PLOT3D Connectivity file: %s.conn\n"%(grid_file))
            
        else:
            print 'invalid file type'
            #raiseError
        #endif
        autofile.write("      Rind layer in solution files: yes\n")
        autofile.write("        Write coordinates in meter: yes\n")
        autofile.write("\n")
        autofile.write("        Automatic parameter update: no\n")
        autofile.write("\n")
        autofile.write("#                Cp curve fit file: \n")
        autofile.write("\n")
        
        autofile.write("       Write precision grid: double\n")
        #autofile.write("# Default is executable precision\n")
        #autofile.write(writeUnit,"(a)") "            # Possibilities: single"
        #autofile.write(writeUnit,"(a)") "            #              : double"
        autofile.write("   Write precision solution: double\n")
        #autofile.write(writeUnit,"(a)") "    # Default is executable precision"
        #autofile.write(writeUnit,"(a)") "            # Possibilities: single"
        #autofile.write(writeUnit,"(a)") "            #              : double"
        autofile.write("\n")#
        #
        autofile.write("Store convergence inner iterations: no\n")
        autofile.write("\n")
        autofile.write("\n")

        #! Write the keywords and default values for the physics parameters.

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Physics Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "                  Equations: %s\n"%(eqn_type))
        #autofile.write(  "            # Possibilities: Euler")
        #autofile.write(  "            #              : Laminar NS")
        #autofile.write(  "            #              : RANS")
        autofile.write("\n")

        autofile.write(  "                       Mode: %s\n"%(sol_type))
        #autofile.write(  "            # Possibilities: Steady")
        #autofile.write(  "            #              : Unsteady")
        #autofile.write(  "            #              : Time spectral")
        autofile.write("\n" )

        autofile.write(  "                  Flow type: External flow\n")
        #autofile.write(  "            # Possibilities: Internal flow")
        #autofile.write(  "            #              : External flow")
        autofile.write( "\n")

        autofile.write(  "                   Cp model: Constant\n")
        #autofile.write(  "        # Other possibility: Temperature curve fits")
        autofile.write( "\n")

        if eqn_type=='RANS':
            autofile.write(  "           Turbulence model: Baldwin Lomax\n")
            #autofile.write(  "            # Possibilities: Baldwin Lomax")
            autofile.write(  "            #              : Spalart Allmaras")
            autofile.write(  "            #              : Spalart Allmaras Edwards")
            autofile.write(  "            #              : KOmega Wilcox")
            autofile.write(  "            #              : KOmega Modified")
            autofile.write(  "            #              : KTau")
            autofile.write(  "            #              : Menter SST")
            autofile.write(  "            #              : V2F")
            autofile.write( "\n")
            
            autofile.write(  "     V2F version (n1 or n6): 1")
            autofile.write(  "        # Other possibility: 6")
            autofile.write("\n" )
            autofile.write(  "       V2F with upper bound: yes")
            autofile.write(  "        # Other possibility: no")
            autofile.write( "\n")

            autofile.write(  " Turbulence production term: Strain")
            autofile.write(  "      # Other possibilities: Vorticity")
            autofile.write(  "            #              : Kato-Launder")
            autofile.write( "\n")

            autofile.write(  "                Use wall functions: no")
            autofile.write(  "Offset from wall in wall functions: 0.0")
            autofile.write( "\n")
        #endif
        
        autofile.write(  "   Constant specific heat ratio: 1.4    #Air\n")
        autofile.write(  "   Gas constant (J/(kg K))     : 287.87 #Air\n")
        autofile.write(  "   Prandtl number              : 0.72   #Air\n")
        autofile.write(  "   Turbulent Prandtl number    : 0.90\n")
        autofile.write(  "   Max ratio k-prod/dest       : 20.0\n")
        autofile.write( "\n")

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Free Stream Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "                             Mach: %3.3e\n"%(aero_problem._flows.mach))
        autofile.write(  "            Mach for coefficients: %3.3e\n"%(aero_problem._flows.mach))
        autofile.write(  "                          # Default is Mach\n")
        autofile.write(  "#                         Reynolds: 100000\n")
        autofile.write(  "       Reynolds length (in meter): 1.0\n")
        #autofile.write(  "   Free stream velocity direction: 1.0 0.05 0.0\n")
        #autofile.write(  "                   Lift direction: -0.05 1.0 0.0\n")
        autofile.write(  "   Free stream velocity direction: %12.12e %12.12e %12.12e\n"%(velDir[0],velDir[1],velDir[2]))
        autofile.write(  "                   Lift direction: %12.12e %12.12e %12.12e\n"%(liftDir[0],liftDir[1],liftDir[2]))
        autofile.write(  "     # Default is normal to free stream without y-component\n")
        autofile.write(  "   Free stream temperature (in K): 288.15\n")
        autofile.write(  " Free stream eddy viscosity ratio: 0.01\n")
        autofile.write(  "  Free stream turbulent intensity: 0.001\n")
        autofile.write(  "\n")
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Reference State\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "            Reference pressure (in Pa): 101325.0\n")
        autofile.write(  "         Reference density (in kg/m^3): 1.25\n")
        autofile.write(  "          Reference temperature (in K): 273.15\n")
        autofile.write(  " Conversion factor grid units to meter: 1.0\n")
        autofile.write( "\n")
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Geometrical Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "           Reference surface: %2.1f\n"%(aero_problem._refs.sref))
        #autofile.write(  "           Reference surface: %2.1f\n"%(1.0))
        autofile.write(  "            Reference length: %2.1f\n"%(aero_problem._refs.cref))
        autofile.write(  "    Moment reference point x:  %2.1f\n"%(aero_problem._refs.xref))
        autofile.write(  "    Moment reference point y:  %2.1f\n"%(aero_problem._refs.yref))
        autofile.write(  "    Moment reference point z:  %2.1f\n"%(aero_problem._refs.zref))
        autofile.write( "\n")
        
        #! Write the keywords and default values for the discretization
        #! parameters.

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Fine Grid Discretization Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "       Discretization scheme: Upwind\n")
        autofile.write(  "             # Possibilities: Central plus scalar dissipation\n")
        autofile.write(  "             #              : Central plus matrix dissipation\n")
        autofile.write(  "             #              : Central plus CUSP dissipation\n")
        autofile.write(  "             #              : Upwind\n")
        autofile.write( "\n")
        
        if eqn_type=='RANS':
            autofile.write(  "   Order turbulent equations: First order\n")
            autofile.write(  "         # Other possibility: Second order\n")
            autofile.write( "\n")
        #endif
       
        autofile.write(  "              Riemann solver: Roe\n")
        autofile.write(  "       # Other possibilities: Van Leer\n")
        autofile.write(  "       #                    : AUSMDV\n\n")
        
        autofile.write(  "   Limiter                  : No limiter\n")
        autofile.write(  "       # Other possibilities: First order\n")
        autofile.write(  "       #                    : Van Albeda\n")
        autofile.write(  "       #                    : Minmod\n")
        autofile.write( "\n")
        
        autofile.write(  "              Preconditioner: No preconditioner\n")
        autofile.write(  "       # Other possibilities: Turkel\n")
        autofile.write(  "       #                    : Choi Merkle\n")
        autofile.write( "\n")
        
        autofile.write(  "     Wall boundary treatment: Linear extrapolation pressure\n")
        autofile.write(  "       # Other possibilities: Constant pressure\n")
        autofile.write(  "       #                    : Linear extrapolation pressure\n")
        autofile.write(  "       #                    : Quadratic extrapolation pressure\n")
        autofile.write( "\n")
        
        autofile.write(  "  Outflow boundary treatment: Constant extrapolation\n")
        autofile.write(  "        # Other possibility: Linear extrapolation\n")
        autofile.write( "\n")
        
        autofile.write(  "non-matching block to block treatment: NonConservative\n")
        autofile.write(  "                  # Other possibility: Conservative\n")
        autofile.write( "\n")
        
        autofile.write(  "                           Vis2: 0.5\n")
        autofile.write(  "                           Vis4: 0.015625  # 1/64\n")
        autofile.write(  "Directional dissipation scaling: yes\n")
        autofile.write(  "   Exponent dissipation scaling: 0.0\n")
        autofile.write( "\n")
        autofile.write(  "   Total enthalpy scaling inlet: no\n")
        autofile.write( "\n")
        autofile.write(  "      Kappa interpolation value: 0.33333\n")
        autofile.write( "\n")
        
        autofile.write(  "              Vortex correction: no\n")
        autofile.write( "\n")
        
        #! Write the keywords and default values for the unsteady
        #! parameters.
        if sol_type=='Unsteady':
            autofile.write(  "-------------------------------------------------------------------------------")
            autofile.write(  "     Unsteady Parameters")
            autofile.write(  "-----------------------------------------\
            &--------------------------------------")
            
            autofile.write(  "     Time integration scheme: BDF")
            autofile.write(  "       # Other possibilities: explicit \
            &Runge-Kutta")
            autofile.write(  "       # Other possibilities: implicit \
            &Runge-Kutta")
            
            autofile.write(  "      Time accuracy unsteady: Second")
            autofile.write(  "       # Other possibilities: \
            &First to Fifth")
            autofile.write( "\n")
            
            autofile.write(  "Number of unsteady time steps coarse grid: \
            &-1  # Means same as on fine grid")
            autofile.write(  "  Number of unsteady time steps fine grid: \
            &MISSING PARAMETER")
            autofile.write( "\n")
            
            autofile.write(  "              Unsteady time step (in sec): \
            &MISSING PARAMETER")
            autofile.write( "\n")
            
            autofile.write(  "      Update wall distance unsteady mode: \
            &yes")
            autofile.write( "\n")
        #endif

        if eqn_type=='Time_Spectral':
            #! Write the keywords and default values for the time spectral
            #! parameters.

            autofile.write("-------------------------------------------------------------------------------\n")
##             write(writeUnit,"(a)") "     Time Spectral Parameters"
##             write(writeUnit,"(a)") "-----------------------------------------&
##             &--------------------------------------"
##             write(writeUnit,"(a)") "             Number time intervals spectral: &
##             &MISSING PARAMETER"
##             write(writeUnit,"(a)")
##             write(writeUnit,"(a)") "            Write file for unsteady restart: &
##             &no"
##             write(writeUnit,"(a)") "    Time step (in sec) for unsteady restart: &
##             &MISSING PARAMETER"
##             write(writeUnit,"(a)")
##             write(writeUnit,"(a)") "       Write unsteady volume solution files: &
##             &no"
##             write(writeUnit,"(a)") "      Write unsteady surface solution files: &
##             &no"
##             write(writeUnit,"(a)") "          Number of unsteady solution files: &
##             &MISSING PARAMETER"
##             write(writeUnit,"(a)")
        #endif

        #! Write the keywords and default values for the iteration
        #! parameters.

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Iteration Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "                                   Smoother: Runge Kutta\n")
        autofile.write(  "                            # Possibilities: Runge Kutta\n")
        autofile.write(  "                            #              : Nonlinear LUSGS\n")
        autofile.write(  "                            #              : Nonlinear LUSGS Line\n")
        autofile.write( "\n")

        autofile.write(  "               Number of Runge Kutta stages: 5\n")
        autofile.write( "\n")

        if eqn_type =='RANS':
            autofile.write(  "              Treatment turbulent equations: Segregated\n")
            autofile.write(  "                        # Other possibility: Coupled\n")
            autofile.write(  "    Number additional turbulence iterations: 0\n")
            autofile.write( "\n")

            autofile.write(  "                         Turbulent smoother: ADI\n")
            autofile.write(  "                        # Other possibility: GMRES\n")
            autofile.write( "\n")
        #endif
       
        autofile.write(  "                        Update bleeds every: 5\n")
        autofile.write(  "Relaxation factor bleed boundary conditions: 0.1\n")

        autofile.write(  "                         Residual averaging: all stages\n")
        autofile.write(  "                      # Other possibilities: no\n")
        autofile.write(  "                      #                    : alternate stages\n")
        autofile.write(  "     Residual averaging smoothing parameter: 1.5\n")

        autofile.write(  "                 Number of multigrid cycles: 100\n")
        autofile.write(  "   Number of single grid startup iterations: 0\n")
        autofile.write(  "                                 Save every: 0\n")
        autofile.write(  "                         Save surface every: 0\n")
        autofile.write(  "                                 CFL number: 1.5\n")
        autofile.write( "\n")
        if eqn_type=='RANS':
            autofile.write(  "                       Turbulent relaxation: Explixit\n")
            autofile.write(  "                            # Possibilities: Explicit\n")
            autofile.write(  "                            #              : Implicit\n")
            autofile.write(  "                     Alpha turbulent DD-ADI: 0.8\n")
            autofile.write(  "                      Beta turbulent DD-ADI: -1  # Same as alpha\n")
        #endif
        autofile.write(  "           Relative L2 norm for convergence: 1.e-12\n")
        autofile.write( "\n")

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Multigrid Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "      Number of multigrid cycles coarse grid:  -1  # Means same as on fine grid\n")
        autofile.write(  "                      CFL number coarse grid: -1.0  # Means same as on fine grid\n")
        autofile.write(  "Relative L2 norm for convergence coarse grid: 1.e-2\n")
        autofile.write( "\n")
        
        autofile.write(  "#        Discretization scheme coarse grid:  # Default fine grid scheme\n")
        autofile.write(  "#               Riemann solver coarse grid:  # Default fine grid solver\n")
        autofile.write(  "                         Vis2 coarse grid: 0.5\n")
        autofile.write( "\n")
        
        autofile.write(  "        Freeze turbulent source terms in MG: yes\n")
        autofile.write(  "                        # Other possibility: no\n")
        autofile.write( "\n")

        autofile.write(  " Treatment boundary multigrid corrections: Zero Dirichlet\n")
        autofile.write(  "            Restriction relaxation factor: 1.0\n")
        autofile.write(  "#                    Multigrid start level:  # Default is coarsest MG level\n")
        autofile.write(  "                 Multigrid cycle strategy: 2v\n")
        autofile.write( "\n")

        #! Write the keywords and default values for the parallel, i.e.
        #! load balance parameters.

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Load balancing Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "        Allowable load imbalance: 0.1\n")
        autofile.write(  "   Split blocks for load balance: no#yes\n")
        autofile.write( "\n")

##        ! Write the visualization parameters.

##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( ) "     Visualization Parameters"
##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( ) "          pV3 visualization only: no"
##        autofile.write( )

##        ! Write the keywords and example values for the grid motion.

##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( ) "     Grid motion Parameters"
##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"

##        autofile.write( ) "     Rotation point body (x,y,z): 0.0 0.0 0.0"
##        autofile.write( )

##        autofile.write( ) "    Degree polynomial x-rotation: 0"
##        autofile.write( ) "    Degree polynomial y-rotation: 0"
##        autofile.write( ) "    Degree polynomial z-rotation: 1"
##        autofile.write( )
       
##        autofile.write( ) "Polynomial coefficients x-rotation: 0.0"
##        autofile.write( ) "Polynomial coefficients y-rotation: 0.0"
##        autofile.write( ) "Polynomial coefficients z-rotation: 0.0 1.e-3"
##        autofile.write( )

##        autofile.write( ) "       Degree fourier x-rotation: 0"
##        autofile.write( ) "       Degree fourier y-rotation: 0"
##        autofile.write( ) "       Degree fourier z-rotation: 1"
##        autofile.write( )

##        autofile.write( ) "        Omega fourier x-rotation: 0.25"
##        autofile.write( ) "        Omega fourier y-rotation: 0.32"
##        autofile.write( ) "        Omega fourier z-rotation: 0.41"
##        autofile.write( )

##        autofile.write( ) "Fourier cosine coefficients x-rotation: 0.0"
##        autofile.write( ) "Fourier cosine coefficients y-rotation: 0.0"
##        autofile.write( ) "Fourier cosine coefficients z-rotation: 0.0 0.0"
##        autofile.write( )

##        autofile.write( ) "Fourier sine coefficients z-rotation: 0.1"
##        autofile.write( )

        #! Write the monitor, surface output and volume output variables.

        autofile.write("-------------------------------------------------------------------------------\n")
        autofile.write( "     Monitoring and output variables\n")
        autofile.write( "-------------------------------------------------------------------------------\n")
        autofile.write( "                Monitoring variables: resrho_cl_cd_cmx\n")
        autofile.write( " Monitor massflow sliding interfaces: no\n")
        autofile.write( "            Surface output variables: rho_cp_vx_vy_vz_mach\n")
        autofile.write( "           Volume output variables: ptloss_resrho\n")
        autofile.write( "\n")

##        ! The section to overwrite the rotation info for the families.

##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( ) "     Family rotation info "
##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( )

##        autofile.write( ) "                               Rotation &
##                               &center  Rotation rate (rad/s)"
##        autofile.write( ) "Rotating family <family_name1> : 0.0 0.0 0.0 &
##                               &   1.e+4 0.e+0 0.e+0"
##        autofile.write( ) "Rotating family <family_name2> : 0.0 0.0 0.0 &
##                               &   1.e+3 0.e+0 0.e+0"
##        autofile.write( ) "Etc."
##        autofile.write( )

##        ! The section to overwrite the boundary condition data sets
##        ! for the families.

##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( ) "     Family boundary condition data sets "
##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( )
##        autofile.write( ) "      # Subsonic outflow boundary with &
##                                   &varying pressure in radial direction"
##        autofile.write( ) "Boundary family <family_name1> : CoordinateR &
##                               & Pressure"
##        autofile.write( ) "            Npoints = 3"
##        autofile.write( ) "                                     0.2     &
##                               &   85360  #Pa"
##        autofile.write( ) "                                     0.6     &
##                               &   88900"
##        autofile.write( ) "                                     1.0     &
##                               &   92730"
##        autofile.write( )

##        autofile.write( ) "      # Subsonic inflow boundary with a &
##                                   &constant state"
##        autofile.write( ) "Boundary family <family_name2> :  &
##                               &PressureStagnation TemperatureStagnation &
##                               &VelocityAngleX VelocityAngleY VelocityAngleZ"
##        autofile.write( ) "            Npoints = 1"
##        autofile.write( ) "                                  &
##                               &    112850                  300.0        &
##                               &      0.0       1.570796327    1.570796327"
##        autofile.write( )

##        autofile.write( ) "      # Isothermal wall boundary with a &
##                               &variable temperature"
##        autofile.write( ) "Boundary family <family_name3> : CoordinateX &
##                               &CoordinateY CoordinateZ Temperature"
##        autofile.write( ) "            Npoints = 2  2  # 2D structured"
##        autofile.write( ) "                                     0.0     &
##                               &    0.0         0.0         312.5  # (1,1)"
##        autofile.write( ) "                                     0.5     &
##                               &    0.0         0.1         315.0  # (2,1)"
##        autofile.write( ) "                                     0.0     &
##                               &    0.5         0.0         310.5  # (1,2)"
##        autofile.write( ) "                                     0.5     &
##                               &    0.5         0.1         313.5  # (2,2)"
##        autofile.write( )

##        autofile.write( ) " # Whether or not to monitor the mass flow&
##                               & for certain families."
##        autofile.write( ) "Boundary family <family_name4> : &
##                               & monitor mass flow: yes   #no"
        autofile.write( "\n")

 
        autofile.write("\n" )
        
        autofile.flush()
        
        autofile.close()
        
        return
    
    def GetMesh(self):
        """Return a reference to the mesh object for this solution."""
        return self.Mesh

    def RunIterations(self,*ncycles):
        """Run the flow solver for a given number of time steps and
           multigrid cycles.

        Keyword arguments:

        ntimesteps -- passed as ncycles[0] is the number of unsteady
                      physical time steps which defaults to the value
                      read from the input file

        ncycles -- passed as ncycles[1], which implies that the number
                   of time steps has to be passed as ncycles[0], is the
                   number of multigrid cycles to run on the finest mesh,
                   which defaults to the previous value used when calling this
                   method OR the value read from the input file if this method
                   has not been called before

        """
        # print ncycles,sumb.monitor.niterold, sumb.monitor.nitercur, sumb.iteration.itertot ,'storeconv',sumb.inputio.storeconvinneriter,True,False
        #sumb.inputio.storeconvinneriter=True#False
        
        if (sumb.monitor.niterold == 0 and
            sumb.monitor.nitercur == 0 and
            sumb.iteration.itertot == 0):
            # No iterations have been done
            if (ncycles):
                # Set new value of unsteady physical time steps
                if (ncycles[0] != sumb.inputunsteady.ntimestepsfine):
                    sumb.inputunsteady.ntimestepsfine = ncycles[0]
                # Set new value of MG cycles
                if (len(ncycles) > 1):
                  if (ncycles[1] != sumb.inputiteration.ncycles):
                      sumb.inputiteration.ncycles = ncycles[1]
                # Reallocate convergence history array and
                # time array with new size
                if (self.myid == 0):
                    sumb.monitor.timearray = None
                    sumb.monitor.timedataarray = None
                    sumb.alloctimearrays(sumb.inputunsteady.ntimestepsfine)
                    if (sumb.inputio.storeconvinneriter):
                        nn = sumb.inputiteration.nsgstartup+sumb.inputiteration.ncycles
                        if (sumb.inputphysics.equationmode==2):#2 is unsteady
                            nn=sumb.inputunsteady.ntimestepsfine*sumb.inputiteration.ncycles
                        #endif
                        sumb.monitor.convarray = None
                        #print 'nn',nn
                        sumb.allocconvarrays(nn)
                        #print 'convarray',sumb.monitor.convarray
                        

        elif (sumb.monitor.nitercur == 0 and
              sumb.iteration.itertot == 0):

            # Read in a restart file but no new iterations
            if (ncycles):
                # Set new value of unsteady physical time steps
                if (ncycles[0] != sumb.inputunsteady.ntimestepsfine):
                    sumb.inputunsteady.ntimestepsfine = ncycles[0]
                # Set new value of MG cycles
                if (len(ncycles) > 1):
                    if (ncycles[1] != sumb.inputiteration.ncycles):
                        sumb.inputiteration.ncycles = ncycles[1]
                # Reallocate convergence history array and
                # time array with new size, storing old values from restart
                if (self.myid == 0):
                    # number of time steps from restart
                    ntimestepsrestart = sumb.monitor.ntimestepsrestart[0]
                    # store restart time history and deallocate arrays
                    temp_t = copy.deepcopy(sumb.monitor.timearray[:ntimestepsrestart])
                    sumb.monitor.timearray = None
                    temp_td = copy.deepcopy(sumb.monitor.timedataarray[:ntimestepsrestart,:])
                    sumb.monitor.timedataarray = None
                    # allocate time history arrays with new extended size
                    sumb.alloctimearrays(ntimestepsrestart+sumb.inputunsteady.ntimestepsfine)
                    # recover values from restart and deallocate temporary arrays
                    sumb.monitor.timearray[:temp_td.shape[0]] = temp_t
                    sumb.monitor.timedataarray[:temp_td.shape[0],:temp_td.shape[1]] = temp_td
                    temp_t = None
                    temp_td = None

                    if (sumb.inputio.storeconvinneriter):
                        # number of iterations from restart
                        niterold = sumb.monitor.niterold[0]
                        # store restart convergence history and deallocate array
                        temp = copy.deepcopy(sumb.monitor.convarray[:niterold+1,:])
                        sumb.monitor.convarray = None
                        # allocate convergence history array with new extended size
                        sumb.allocconvarrays(temp.shape[0]
                                             +sumb.inputunsteady.ntimestepsfine
                                             *sumb.inputiteration.ncycles-1)
                        # recover values from restart and deallocate temporary array
                        sumb.monitor.convarray[:temp.shape[0],:temp.shape[1]] = temp
                        temp = None
                  
        else:
            # More Time Steps / Iterations in the same session
            #print 'more cycles',ncycles[0],ncycles[1],sumb.monitor.convarray,self.myid
            
            if (ncycles):
                # Set new value of unsteady physical time steps to run
                if (ncycles[0] != sumb.inputunsteady.ntimestepsfine):
                    sumb.inputunsteady.ntimestepsfine = ncycles[0]
                # Set new value of MG cycles ro run
                if (len(ncycles) > 1):
                    if (ncycles[1] != sumb.inputiteration.ncycles):
                        # print 'sumb cycles',sumb.inputiteration.ncycles
                        sumb.inputiteration.ncycles = ncycles[1]
                        #print sumb.inputiteration.ncycles

            # Reallocate convergence history array and
            # time array with new size, storing old values from previous runs
            if (self.myid == 0):
                #print 'more cycles',ncycles[0],ncycles[1],sumb.monitor.convarray
                
                #print 'equation mode',sumb.inputphysics.equationmode
                if (sumb.inputphysics.equationmode==2):#2 is unsteady
                    #print 'unsteady'
                    #print 'checking time array'
                    # store previous time history and deallocate arrays
                    temp_t = copy.deepcopy(sumb.monitor.timearray)
                    #print 'deallocating time array',sumb.monitor.timearray
                    sumb.monitor.timearray = None
                    #print 'checking time data array',sumb.monitor.timedataarray
                    temp_td = copy.deepcopy(sumb.monitor.timedataarray)
                    #print 'deallocating timedata array'
                    sumb.monitor.timedataarray = None
                    #print 'allocating new time array'
                    # allocate time history arrays with new extended size
                    sumb.alloctimearrays(temp_td.shape[0]+sumb.inputunsteady.ntimestepsfine)
                    # recover values from previous runs and deallocate temporary arrays
                    sumb.monitor.timearray[:temp_td.shape[0]] = temp_t
                    sumb.monitor.timedataarray[:temp_td.shape[0],:temp_td.shape[1]] = temp_td
                    temp_t = None
                    temp_td = None

                if (sumb.inputio.storeconvinneriter):
                    #print 'store conv?',sumb.inputio.storeconvinneriter,sumb.inputiteration.ncycles-1
                    #print 'conv, arreay',sumb.monitor.convarray[:,0,1]
                    #sdfg
                    # store previous convergence history and deallocate array
                    temp = copy.deepcopy(sumb.monitor.convarray)
                    #temp = copy.deepcopy(sumb.monitor.convarray)
                    #print 'conv, arreay',sumb.monitor.convarray
                    sumb.monitor.convarray = None
                    #print 'allocating convergence arrays for new size',sumb.monitor.convarray,temp.shape[0],sumb.inputunsteady.ntimestepsfine,sumb.inputiteration.ncycles,sumb.inputiteration.ncycles-1
                    
                    #print 'testing',temp.shape[0]
                    # allocate convergence history array with new extended size
                    nn = sumb.inputiteration.nsgstartup+sumb.inputiteration.ncycles
                    if (sumb.inputphysics.equationmode==2):#2 is unsteady
                        nn=sumb.inputunsteady.ntimestepsfine*sumb.inputiteration.ncycles
                    #endif
                    sumb.allocconvarrays(temp.shape[0]+nn-1)
                    #print 'convergence shape',sumb.monitor.convarray.shape,temp.shape, sumb.monitor.convarray[:temp.shape[0],:].shape,temp[:,0,:].shape
                    # recover values from previous runs and deallocate temporary array
                    sumb.monitor.convarray[:temp.shape[0],:] = copy.deepcopy(temp)
##                     # allocate convergence history array with new extended size
##                     sumb.allocconvarrays(temp.shape[0]
##                                          +sumb.inputunsteady.ntimestepsfine
##                                          *sumb.inputiteration.ncycles-1)
##                     # recover values from previous runs and deallocate temporary array
##                     sumb.monitor.convarray[:temp.shape[0],:temp.shape[1]] = temp
                    temp = None

            # re-initialize iteration variables
	    sumb.inputiteration.mgstartlevel = 1
            sumb.monitor.niterold  = sumb.monitor.nitercur
            sumb.monitor.nitercur  = 0
            sumb.iteration.itertot = 0
            # update number of time steps from restart
            sumb.monitor.ntimestepsrestart = sumb.monitor.ntimestepsrestart \
                                           + sumb.monitor.timestepunsteady
            # re-initialize number of time steps previously run (excluding restart)             
            sumb.monitor.timestepunsteady = 0
            # update time previously run
            sumb.monitor.timeunsteadyrestart = sumb.monitor.timeunsteadyrestart \
                                             + sumb.monitor.timeunsteady
            # re-initialize time run
            sumb.monitor.timeunsteady = 0.0

        #endif
        if self.myid ==0: print 'calling solver'
        self.GetMesh()._UpdateGeometryInfo()
        sumb.solver()
## ################################################################
## # The code below reproduces solver.F90 in python
## # See file solver.F90
## ################################################################

##         if (sumb.inputdiscretization.wallbctreatment == sumb.inputdiscretization.normalmomentum):
##            sumb.iteration.exchangepressureearly = sumb.eulerwallspresent()
##         else:
##            sumb.iteration.exchangepressureearly = False

##         sumb.killsignals.localsignal == sumb.killsignals.nosignal
##         try:
##             sumb.connectsignals()
##         except AttributeError:
##             pass

##         sumb.iteration.t0solver = sumb.su_wtime()

##         sumb.monitor.timeunsteady = 0.0

##         if (not sumb.iteration.pv3initialized):
##             sumb.iteration.groundlevel = 1
##             try:
##                 sumb.initalizepv3()
##                 sumb.iteration.pv3initialized = True
##             except AttributeError:
##                 pass
        
##         for sumb.iteration.groundlevel in range(sumb.inputiteration.mgstartlevel,0,-1):
##             if (sumb.inputphysics.equationmode == sumb.inputphysics.steady or \
##                 sumb.inputphysics.equationmode == sumb.inputphysics.timespectral):
##                 sumb.solversteady()
##             elif (sumb.inputphysics.equationmode == sumb.inputphysics.unsteady):
## #                sumb.solverunsteady()
## ################################################################
## # The code below reproduces solverUnsteady.F90 in python
## # See file solverUnsteady.F90
## ################################################################
## #               sumb.monitor.writevolume  = False
## #               sumb.monitor.writesurface = False
## #               sumb.monitor.writegrid    = False

##                 sumb.monitor.timestepunsteady = 0

##                 ntimesteps = sumb.inputunsteady.ntimestepscoarse
##                 if (sumb.iteration.groundlevel == 1):
##                     ntimesteps =  sumb.inputunsteady.ntimestepsfine

##                 for iter in range(1,ntimesteps+1):
##                     sumb.inittimesteppart1()
##                     ##############################  
##                     #sumb.additional_routine_here#
##                     ##############################
##                     sumb.inittimesteppart2()

##                     sumb.solvestate()
## #                   sumb.iteration.noldsolavail = sumb.iteration.noldsolavail + 1

## #ifdef USE_PV3
## #                   CALL PV_UPDATE(REAL(TIMESTEPUNSTEADY,REALPV3TYPE))
## #endif

##                     sumb.checkwriteunsteadyinloop()

##                     if (sumb.killsignals.globalsignal == sumb.killsignals.signalwritequit):
##                         break

##                 sumb.checkwriteunsteadyendloop()
## ################################################################
## # end of solverUnsteady.F90
## ################################################################

##             if (sumb.iteration.groundlevel > 1):

##                 sumb.iteration.currentlevel = sumb.iteration.groundlevel - 1

##                 if (sumb.communication.myid == 0):
##                     print "#"
##                     print "# Going down to grid level %i" % sumb.iteration.currentlevel
##                     print "#"

##                 sumb.transfertofinegrid(False)

##                 if (sumb.inputphysics.equationmode == sumb.inputphysics.unsteady and \
##                     sumb.iteration.changing_grid):
##                     sumb.updatecoorfinemesh(sumb.monitor.timeunsteady, 1)

##                 sumb.iteration.noldsolavail = 1
## ################################################################
## # end of solver.F90
## ################################################################

    def CheckIfParallel(self):
        """Return true if we are running with MPI."""
        return _parallel,mpi
        
    def WriteVolumeSolutionFile(self,*filename):
        """Write the current state of the volume flow solution to a CGNS file.
 
        Keyword arguments:
        
        filename -- the name of the file (optional)

        """
        if (filename):
	   sumb.inputio.solfile[:] = ''
           sumb.inputio.solfile[0:len(filename[0])] = filename[0]
        if(sumb.iteration.changing_grid or sumb.inputmotion.gridmotionspecified):
            sumb.monitor.writegrid=True
	else:
            sumb.monitor.writegrid=False
        sumb.monitor.writegrid=True
        sumb.monitor.writevolume=True
        sumb.monitor.writesurface=False
        sumb.writesol()

    def WriteSurfaceSolutionFile(self,*filename):
        """Write the current state of the surface flow solution to a CGNS file.
 
        Keyword arguments:
        
        filename -- the name of the file (optional)

        """
        if (filename):
            sumb.inputio.surfacesolfile[:] = ''
            sumb.inputio.surfacesolfile[0:len(filename[0])] = filename[0]
        sumb.monitor.writegrid=False
        sumb.monitor.writevolume=False
        sumb.monitor.writesurface=True
        sumb.writesol()

    def GetSurfaceLoads(self, family=None, sps=1):
        """Return an array of the surface forces.
         
        Keyword arguments:
         
        family -- optional string specifying the return of forces
                  only for the specified family.
        sps    -- spectral time step (optional, default is set to 1).
                  (sps=1 for usual steady or unsteady models)
                                   
        """
#        if(sumb.flovarrefstate.viscous):
#            sumb.iteration.rfil = 1.
#            for i in range(sumb.block.ndom):
#                sumb.setpointers(i,1)
#                sumb.viscousflux() 
        if(family):
            try:
                index = self.Mesh.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdcreatesurfforcelist(sps,index)
            return sumb.mddata.mdsurfforce[:,start_ind-1:end_ind]
        else:
            nfamilies = len(self.Mesh.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdcreatesurfforcelist(sps,0)
            else:
                for n in range(nfamilies):
                    [start_ind,end_ind] = sumb.mdcreatesurfforcelist(sps,n+1)
            return sumb.mddata.mdsurfforce

    def DeallocateSurfaceLoads(self):
        """Deallocate memory used for the surface loads."""
        sumb.mddeletesurfforcelist()

    def GetSurfaceLoadsLocal(self, family=None, sps=1):
        """Return an array of the surface forces.
         
        Keyword arguments:
         
        family -- optional string specifying the return of forces
                  only for the specified family.
        sps    -- spectral time step (optional, default is set to 1).
                  (sps=1 for usual steady or unsteady models)
                                   
        """
#        if(sumb.flovarrefstate.viscous):
#            sumb.iteration.rfil = 1.
#            for i in range(sumb.block.ndom):
#                sumb.setpointers(i,1)
#                sumb.viscousflux() 
        if(family):
            try:
                index = self.Mesh.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdcreatesurfforcelistlocal(sps,index)
            return sumb.mddatalocal.mdsurfforcelocal[:,start_ind-1:end_ind]
        else:
            nfamilies = len(self.Mesh.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdcreatesurfforcelistlocal(sps,0)
            else:
                for n in range(nfamilies):
                    [start_ind,end_ind] = sumb.mdcreatesurfforcelistlocal(sps,n+1)
            return sumb.mddatalocal.mdsurfforcelocal

    def DeallocateSurfaceLoadsLocal(self):
        """Deallocate memory used for the surface loads."""
        sumb.mddeletesurfforcelistlocal()

    def AccumulateLoads(self,oml_loads_local):
        '''
        Sum up all of the local loads to get the total force on the oml
        '''
        oml_loads = self.sumb_comm_world.Allreduce(oml_loads_local,mpi.SUM)

        return oml_loads


    def GetSurfaceCp(self, family=None, sps=1):
        """Return an array of the surface pressure coefficients.

        Keyword arguments:

        family -- optional string specifying the return of surface Cp
                  only for the specified family.
        sps    -- spectral time step (optional, default is set to 1).
                  (sps=1 for usual steady or unsteady models)

        """
        if(family):
            try:
                index = self.Mesh.families[family]
            except KeyError:
                print "Error: No such family '%s'" % family
                return None
            [start_ind,end_ind] = sumb.mdcreatesurfvarlist(sps,index)
            return sumb.mddata.mdsurfval[start_ind-1:end_ind]
        else:
            nfamilies = len(self.Mesh.families)
            if (nfamilies == 0):
                [start_ind,end_ind] = sumb.mdcreatesurfvarlist(sps,0)
            else:
                for n in range(nfamilies):
                    [start_ind,end_ind] = sumb.mdcreatesurfvarlist(sps,n+1)
            return sumb.mddata.mdsurfval

    def DeallocateSurfaceCp(self):
        """Deallocate memory used for the surface pressure coefficients."""
        sumb.mddeletesurfvallist()

    def GetMach(self):
        """Get the current freestream Mach number."""
        return sumb.inputparam.mach

    def SetMach(self,mach):
        """Set the freestream Mach number."""
        sumb.inputparam.mach = mach

    def GetConvergenceHistory(self,name):
        """Return an array of the convergence history for a particular quantity.

        Keyword arguments:

        name -- the text string for a particular quantity

        """
        try:
            index = self.monnames[name]
        except KeyError:
            print "Error: No such quantity '%s'" % name
            return None
	if (self.myid == 0):
            if (sumb.monitor.niterold == 0 and
                sumb.monitor.nitercur == 0 and
                sumb.iteration.itertot == 0):
	        history = None
            elif (sumb.monitor.nitercur == 0 and
                  sumb.iteration.itertot == 0):
	        niterold = sumb.monitor.niterold[0]	    
                history = sumb.monitor.convarray[:niterold+1,index]
            else:
	        history = sumb.monitor.convarray[:,index]
	else:
            history = None
	history = self.sumb_comm_world.bcast(history)
	return history

    def GetMonitoringVariables(self):
        """Return a list of the text strings describing the variables being
        monitored.

        """
        return self.monnames.keys()

    ###########################
    #ADjoint routines
    ###########################

    def initializeADjoint(self):
        '''
        Initialize the Ajoint problem for this test case
        in SUMB
        '''
	#Check to see if initialization has already been performed
        if(self.adjointInitialized):
            return
        
        #Set the mesh level and timespectral instance for this
        #computation
        self.level = 1
        self.sps = 1
        
        sumb.iteration.currentlevel=1
        sumb.iteration.groundlevel=1
        
        #Run the preprocessing routine. Sets the node numbering and
        #allocates memory.
        print 'preprocessing adjoint'
        sumb.preprocessingadjoint(self.level)
        
        #Initialize the design variable and function storage
        print 'Before design init'
        sumb.designinit()

        #initalize PETSc
        print 'before petsc'
        sumb.initializepetsc()

        #create the neccesary PETSc objects
        print 'before createpetsecars'
        sumb.createpetscvars()

        #mark the ADjoint as initialized
        self.adjointInitialized=True
        if(self.myid==0):
            print 'ADjoint Initialized Succesfully...'
        #endif
        print 'before nspatial'
        self.nSpatial = sumb.adjointvars.ndesignspatial

        return

    def setupADjointMatrix(self):
        '''
        Setup the ADjoint dRdw matrix and create the PETSc
        Solution KSP object
        '''
        sumb.setupadjointmatrix(self.level)

        sumb.createpetscksp()

        return

    def releaseAdjointMemeory(self):
        '''
        release the KSP memory...
        '''
        sumb.destroypetscksp()

        return

    def setupADjointRHS(self,objective):
        '''
        setup the RHS vector for a given cost function
        '''
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        if self.myid==0:
            print 'Analyzing ADjoint for costfuntion: ',objective
            print 'SUmb index:',SUmbCostfunctions[objective]
    
        #print SUmbCostfunctions[objective]

        sumb.setupadjointrhs(self.level,self.sps,SUmbCostfunctions[objective])

        return

    def solveADjointPETSc(self):
        '''
        Solve the ADjoint system using PETSc
        '''
        sumb.solveadjointpetsc()

        return

    def setupGradientRHSVolume(self,objective):
        '''
        setup the rhs partial for the mesh derivatives
        '''
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
				       'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
				       'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
				       'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
				       'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
				       'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
				       'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
				       'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
				       }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing RHS Partial for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        sumb.setupgradientrhsspatial(self.level,SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return

    
    def setupGradientRHSFlow(self,objective):
        '''
        setup the rhs partial for the mesh derivatives
        '''
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
				       'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
				       'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
				       'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
				       'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
				       'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
				       'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
				       'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
				       }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing RHS Partial for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        sumb.setupgradientrhsextra(self.level,SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return
    
        
    def setupGradientMatrixVolume(self):

        """Set up the residual sensitivity w.r.t. spatial design variables."""

        sumb.setupgradientmatrixspatial(self.level)

	if (self.myid == 0):
            print "ADjoint: Spatial residual sensitivity set up successfully."

        return

    def setupGradientMatrixFlow(self):

        """Set up the residual sensitivity w.r.t. spatial design variables."""

        sumb.setupgradientmatrixextra(self.level)

	if (self.myid == 0):
            print "ADjoint: Extra Vars residual sensitivity set up successfully."

        return
    
    def computeTotalVolumeDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
				       'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
				       'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
				       'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
				       'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
				       'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
				       'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
				       'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
				       }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing total mesh derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        sumb.computeadjointgradientspatial(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return
        
    def getTotalVolumeDerivatives(self,objective):

        """
        Get the  sensitivities from SUmb.
                
	"""
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               }
        
        grad = numpy.zeros((sumb.adjointvars.ndesignspatial),float)
        grad[:] = sumb.adjointvars.functiongradspatial[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for item in objective:
        
        #for i in xrange(sumb.adjointvars.ndesignspatial):
            #grad[i] = sumb.adjointvars.functiongradspatial[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
           
        #endfor
        #endfor
        
        return grad

    def computeTotalFlowDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
				       'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
				       'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
				       'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
				       'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
				       'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
				       'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
				       'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
				       }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing total flow derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        sumb.computeadjointgradientextra(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return

    def getTotalFlowDerivatives(self,objective):

        """
        Get the  sensitivities from SUmb.
                
	"""
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               }
        
        grad = numpy.zeros((sumb.adjointvars.ndesignextra),float)
        grad[:] = sumb.adjointvars.functiongrad[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for item in objective:
        #for i in xrange(sumb.adjointvars.ndesignextra):
            #grad[i] = sumb.adjointvars.functiongrad[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
        #endfor
        #endfor
        
        return grad

    def computeAeroCouplingDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
				       'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
				       'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
				       'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
				       'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
				       'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
				       'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
				       'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
				       }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing Aero Coupling derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        sumb.computeaerocoupling(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return
        
    def getAeroCouplingDerivatives(self,objective):

        """
        Get the  sensitivities from SUmb.
                
	"""
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               }
        
        grad = numpy.zeros((sumb.adjointvars.ndesignspatial),float)
        #for item in objective:
        grad[:] = sumb.adjointvars.functiongradcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for i in xrange(sumb.adjointvars.ndesignspatial):
            #grad[i] = sumb.adjointvars.functiongradcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
            
        #endfor
        #endfor
        
        return grad

    def computeAeroExplicitCouplingDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
				       'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
				       'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
				       'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
				       'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
				       'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
				       'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
				       'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
				       }
        try:
            #for item in objective:
            if self.myid==0:
                print 'Computing Aero Coupling derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
                print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
            #endif

            sumb.computeaeroexpcoupling(SUmbCostfunctions[possibleObjectives[objective]])
        except:
            print 'not an aerodynamic cost function'
        #end

        return
        
    def getAeroExplicitCouplingDerivatives(self,objective):

        """
        Get the  sensitivities from SUmb.
                
	"""
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               }
        
        grad = numpy.zeros((sumb.adjointvars.ndesignspatial),float)
        #for item in objective:
        try:
            grad[:] = sumb.adjointvars.functiongradexpcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
            #for i in xrange(sumb.adjointvars.ndesignspatial):
                #grad[i] = sumb.adjointvars.functiongradexpcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
            #endfor
        except:
            print 'not an aerodynamic cost function'
        #endtry
        
        return grad
    
    def getGlobalNodesLocal(self,blocknum,il,jl,kl):
        #get global node ordering from sumb
        print 'in sumbInterface'
        globalNodes = sumb.getglobalnodes(blocknum,il,jl,kl)
        
        return globalNodes

    def getFunctionValues(self):
        '''
        retrieve the solution values from SUmb
        '''
        #print 'interface getting solution'
        # Map cost functions
        sumb.getsolution()

        #print 'solution mapped'
        SUmbsolutions = {'cl':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncliftcoef-1],\
                         'cd':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncdragcoef-1],\
                         'cFx':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncforcexcoef-1],\
                         'cFy':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncforceycoef-1],\
                         'cFz':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncforcezcoef-1],\
                         'cMx':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncmomxcoef-1],\
                         'cMy':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncmomycoef-1],\
                         'cMz':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncmomzcoef-1],\
                         }

        return SUmbsolutions

    def augmentADjointRHS(self,objective,structAdjoint):
        '''
        run the routines to augment the RHS of the ADjoint
        '''

        sumb.setupcouplingmatrixstruct(1)
        sumb.setupadjointrhsstruct(structAdjoint)
        

        return

    def getAdjoint(self,objective):

        """
        Get the  sensitivities from SUmb.
                
	"""
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               }

        sumb.getadjoint(SUmbCostfunctions[possibleObjectives[objective]])
        nstate = sumb.flowvarrefstate.nw*sumb.adjointvars.ncellsglobal
        adjoint = numpy.zeros((nstate),float)
        #for item in objective:
        #This should work
        adjoint[:] = sumb.adjointvars.adjoint[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for i in xrange(nstate):
            # adjoint[i] = sumb.adjointvars.adjoint[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
        #endfor
       
        return adjoint

    def getTotalStructDerivatives(self,objective):

        """
        Get the  sensitivities from SUmb.
                
	"""
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               }
        
        grad = numpy.zeros((sumb.adjointvars.ndesignspatial),float)
        #for item in objective:
        grad[:] = sumb.adjointvars.functiongradstruct[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for i in xrange(sumb.adjointvars.ndesignspatial):
            #grad[i] = sumb.adjointvars.functiongradstruct[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
        #endfor
        
        return grad
    
    def aeroComputeTotalDerivatveStruct(self,objective,structAdjoint={}):
        '''
        compute the force portion of the total structural derivative
        based on the coupled structural adjoint
        '''
        
        sumb.setupcouplingtotalstruct(1)
        
        SUmbCostfunctions = {'cl':sumb.adjointvars.costfuncliftcoef,\
                             'cd':sumb.adjointvars.costfuncdragcoef,\
                             'cFx':sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':sumb.adjointvars.costfuncforceycoef,\
                             'cFz':sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':sumb.adjointvars.costfuncmomycoef,\
                             'cMz':sumb.adjointvars.costfuncmomzcoef,\
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               }
        try:
            #for item in objective:
            if self.myid==0:
                print 'Computing Aero Coupling derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
                print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
            #endif
                
            sumb.setupadjointtotalstruct(structAdjoint,SUmbCostfunctions[possibleObjectives[objective]])
        except:
            print 'not an aerodynamic cost function'
        #end

        return
    
