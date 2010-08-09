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
import bisect
from math import pi

# =============================================================================
# Extension modules
# =============================================================================

import numpy
try:
    from mdo_import_helper import *
    
except:
    print 'Warning: mdo_import_helper not found'
# end try

# try to import the mpi module
try:
    from mpi4py import MPI as mpi
    _parallel = True
except ImportError:
    try:
        import dummy_mpi as mpi
        _parallel = False
    except ImportError:
        print "Error: Failed to import mpi or dummy_mpi."

# # Import the sumb module


# if _parallel:
#     try:
#         import sumb_parallel as sumb
#     except ImportError:
#         try:
#             import sumb
#             if mpi.COMM_WORLD.rank==0:
#                 print "Warning: Running in an MPI environment, but failed"
#                 print "         to import parallel version of SUmb.  Proceeding"
#                 print "         with a sequential version."
#             #endif

#         except ImportError:
#             if mpi.COMM_WORLD.rank==0:
#                 print "Error: Failed to import parallel or sequential version"
#                 print "       of SUmb."
#             #endif
# else:
#     try:
#         import sumb
#     except ImportError:
#         print "Error: Failed to import sequential version of SUmb."


# =============================================================================

class SUmbMesh(object):
    """Represents a SUmb mesh.
 
    """
 
    def __init__(self,comm,sumb):
        """Initialize the object."""
        self._update_geom_info = False
        self._suggar_interface_initialized = False
        self.comm = comm
        self.myid = comm.rank
        self.sumb = sumb
#     def GetSurfaceCoordinates(self, family=None, sps=1):
#         """Return surface coordinates.
        
#         Keyword arguments:
        
#         family -- optional string specifying the return of surface coordinates
#                   only for the specified family.
#         sps    -- spectral time step (optional, default is set to 1).
#                   (sps=1 for usual steady or unsteady models)

#         """
        
#         if(family):
#             try:
#                 index = self.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdcreatesurfcoorlist(sps,index)
#             return sumb.mddata.mdsurfxx[:,start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdcreatesurfcoorlist(sps,0) 
#             else:
#                 for n in range(nfamilies): 
#                     [start_ind,end_ind] = sumb.mdcreatesurfcoorlist(sps,n+1)
#             return sumb.mddata.mdsurfxx

#     def GetSurfaceCoordinatesLocal(self, family=None, sps=1):
#         """Return surface coordinates for this processor only.
        
#         Keyword arguments:
        
#         family -- optional string specifying the return of surface coordinates
#                   only for the specified family.
#         sps    -- spectral time step (optional, default is set to 1).
#                   (sps=1 for usual steady or unsteady models)

#         """
        
#         if(family):
#             try:
#                 index = self.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdcreatesurfcoorlistlocal(sps,index)
#             return sumb.mddatalocal.mdsurfxxlocal[:,start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdcreatesurfcoorlistlocal(sps,0) 
#             else:
#                 for n in range(nfamilies): 
#                     [start_ind,end_ind] = sumb.mdcreatesurfcoorlistlocal(sps,n+1)
#             return sumb.mddatalocal.mdsurfxxlocal

#     def DeallocateSurfaceCoordinates(self):
#         """Deallocate memory used for the surface coordinates."""
#         sumb.mddeletesurfcoorlist()

          
#     def DeallocateSurfaceCoordinatesLocal(self):
#         """Deallocate memory used for the surface coordinates."""
#         sumb.mddeletesurfcoorlistlocal()

#     def GetSurfacePatchDimensions(self, family=None):
#         """Return surface patch dimensions.

#         Keyword arguments:

#         family -- optional string specifying the return of surface patch
#                   dimensions only for the specified family.

#         """
#         if(family):
#             try:
#                 index = self.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdsurfacepatchdim(index)
#             return sumb.mddata.mdpatchdimensions[:,start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdsurfacepatchdim(0)
#             else:
#                 for n in range(nfamilies):
#                     [start_ind,end_ind] = sumb.mdsurfacepatchdim(n+1)
#             return sumb.mddata.mdpatchdimensions

#     def DeallocateSurfacePatchDimensions(self):
#         """Deallocate memory used for the surface patch dimensions."""
#         sumb.mddeletesurfacepatchdim()
                                             
#     def GetSurfaceIndices(self, family=None):
#         """Return CGNS block IDs and indices.  This is a 2-D array, dimensions
#         (4,npoints).  The first 3 indices are the i,j,k and the fourth is the
#         CGNS block ID.
                                             
#         Keyword arguments:
                                             
#         family -- optional string specifying the return of indices only for the
#                   specified family.
                                             
#         """
#         if(family):
#             try:
#                 index = self.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdcreatesurfindlist(index)
#             return sumb.mddata.mdsurfind[:,start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdcreatesurfindlist(0)
#             else:
#                 for n in range(nfamilies):
#                     [start_ind,end_ind] = sumb.mdcreatesurfindlist(n+1)
#             return sumb.mddata.mdsurfind

#     def GetSurfaceIndicesLocal(self, family=None):
#         """Return CGNS block IDs and indices.  This is a 2-D array, dimensions
#         (4,npoints).  The first 3 indices are the i,j,k and the fourth is the
#         CGNS block ID.
                                             
#         Keyword arguments:
                                             
#         family -- optional string specifying the return of indices only for the
#                   specified family.
                                             
#         """
#         if(family):
#             try:
#                 index = self.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdcreatesurfindlistlocal(index)
#             return sumb.mddata.mdsurfindlocal[:,start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdcreatesurfindlistlocal(0)
#             else:
#                 for n in range(nfamilies):
#                     [start_ind,end_ind] = sumb.mdcreatesurfindlistlocal(n+1)
#             return sumb.mddatalocal.mdsurfindlocal
                                             
#     def DeallocateSurfaceIndices(self):
#         """Deallocate memory used for the surface indices."""
#         sumb.mddeletesurfindlist()

#     def DeallocateSurfaceIndicesLocal(self):
#         """Deallocate memory used for the surface indices."""
#         sumb.mddeletesurfindlistlocal()

    def WriteMeshFile(self,*filename):
        """Write the current state of the mesh to a CGNS file.
 
        Keyword arguments:
        
        filename -- the name of the file (optional)
         
        """

        print 'Writing Mesh File'
        self.sumb.iteration.groundlevel = 1
        if (filename):
            self.sumb.inputio.newgridfile[:] = ''
            self.sumb.inputio.newgridfile[0:len(filename[0])] = filename[0]
        self.sumb.monitor.writegrid=True
        self.sumb.monitor.writevolume=True#False
        self.sumb.monitor.writesurface=True
        self.sumb.writesol()

#     def DummySetCoordinates(self, sps =1):
#         """Dummy SetCoordinates routine to be called on processors where no
#         new coordinates reside.
 
#         Keyword arguments:

#         sps    -- spectral time step (optional, default set to 1).
#                         (sps= 1 for usual steady or unsteady models)

#         """
#         sumb.iteration.groundlevel = 1
#         blocknums = numpy.zeros((0))
#         ranges = numpy.zeros((3,2,0))
#         xyz = numpy.zeros((3,0),'d')
#         sumb.mdsetcoor(sps,blocknums,ranges,xyz)
#         self._update_geom_info = True
 
#     def SetCoordinates(self, blocknums, ranges, xyz, sps=1):
#         """Set the coordinates for a given block or sub-blocks.
                                                                                
#         Keyword arguments:
                                                                                
#         blocknums -- the block number for each sub-block, dimension (nsubblocks)
#         ranges -- the ijk ranges for each sub-block, dimension (3,2,nsubblocks)
#         xyz -- the new xyz coordinates, dimension (3,npts)
#         sps    -- spectral time step (optional, default is set to 1).
#                   (sps=1 for usual steady or unsteady models)
                                                                                
#         """
#         sumb.iteration.groundlevel = 1
#         sumb.mdsetcoor(sps,blocknums,ranges,xyz)
#         self._update_geom_info = True

#     def SetCoordinatesLocal(self, blocknum,il,jl,kl,xyz, sps=1):
#         """Set the coordinates for a given block or sub-blocks.
                                                                                
#         Keyword arguments:
                                                                                
#         blocknums -- the block number for each sub-block, dimension (nsubblocks)
#         ranges -- the ijk ranges for each sub-block, dimension (3,2,nsubblocks)
#         xyz -- the new xyz coordinates, dimension (3,npts)
#         sps    -- spectral time step (optional, default is set to 1).
#                   (sps=1 for usual steady or unsteady models)
                                                                                
#         """
#         sumb.iteration.groundlevel = 1
#         #sumb.mdsetcoor(sps,blocknums,ranges,xyz.real)#only the real part needs to be set in SUmb
#         #print 'setting block coords',il,jl,kl,xyz.real.shape
#         sumb.setblockcoords(blocknum,il,jl,kl,xyz.real)
#         #print 'bloock coords set'
#         #sumb.mdsetcoor(sps,blocknums,ranges,xyz)
#         self._update_geom_info = True

    def _UpdateGeometryInfo(self):
        """Update the SUmb internal geometry info, if necessary."""
        if (self._update_geom_info):
            self.sumb.updatecoordinatesalllevels()
            self.sumb.updatewalldistancealllevels()
            self.sumb.updateslidingalllevels()
            self.sumb.updatemetricsalllevels()
            self.sumb.updategridvelocitiesalllevels()
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
        return string.strip(self.sumb.inputio.gridfile.tostring())

    def _InitializeSuggarInterface(self):
        """Sets up ability to call SUGGAR interface methods."""

        self._suggar_interface_initialized = True
        self.sumb.initsuggarinterface()

        # Create a dictionary of zone names.
        self.zones = {}
        for i in range(self.sumb.suggardata.nzones[0]):
            self.zones[string.strip(self.sumb.suggardata.zonenames[i]
                       .tostring())] = self.sumb.suggardata.unsortedzone[i]

    def LoadSuggarDCIFiles(self, dcifilelist, sps=1):
        """Loads the SUGGAR ouput DCI files in the list to Self.Sumb. The previous
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
                self.sumb.loadsuggardcifile(dcifile, sps)
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
            self.sumb.writeplot3dzonefile(index, filename, sps, byteswap)

        # Wait for all zone files to be written before returning since the
        # next step will be to run SUGGAR.
        self.sumb_comm_world.barrier

    def getBlockDimensions(self,blocknum):
        """ Get the i,j,k dimensions of block blocknum"""
        return  self.sumb.getblockdims(blocknum)

    def getBlockCoordinates(self,blocknum,il,jl,kl):
        """get the xyz coordinates from blockblocknum"""
        return self.sumb.getblockcoords(blocknum,il,jl,kl)

    def getBlockCGNSID(self,blocknum):
        """ Get original CGNS blockID for block blocknum"""
        return self.sumb.getblockcgnsid(blocknum)

    def getNSubfacesBlock(self,blocknum):
        """ get subface infor for block: blocknum"""
        return self.sumb.getnsubfacesblock(blocknum)

    def getBlockCommunicationInfo(self,blocknum,nSubface,n1to1,nNonMatch):
        """Get all fo the relevant MPI communication ifo for this block"""
        return self.sumb.getblockcommunicationinfo(blocknum,nSubface,n1to1,nNonMatch)

    def getNBlocksLocal(self):
        """ get the number of blocks present on the local processor"""
        return self.sumb.getnblockslocal()

    def getSingleState(self,blocknum,i,j,k,l):
        '''
        get the requested state value
        '''
        return self.sumb.getsinglestate(blocknum,i,j,k,l)

    def setSingleState(self,blocknum,i,j,k,l,state):
        '''
        set the requested state value
        '''
        self.sumb.setsinglestate(blocknum,i,j,k,l,state)
        
        return 
##     def internalWarping(self,new_cfd_surf,indices):
##         '''
##         Call the integrated warping routines...
##         '''
##         self.sumb.integratedwarp(new_cfd_surf,indices)

##         return

    def initializeInternalWarping(self):
        '''
        Initialize the required variables for the internal
        Meshwarping and derivatives
        '''
        if(self.sumb.cgnsgrid.cgnsnfamilies>0):
            famID = 1
        else:
            famID = 0
        #endif
        
        self.sumb.initializewarping(famID)

        self.nGlobalSurfNodes = self.sumb.mddata.mdnsurfnodescompact
        return

    def SetGlobalSurfaceCoordinates(self,xyz,reinitialize=True):
        """Set the surface coordinates for the mesh                                                                                
        Keyword arguments:
                                                                                
        xyz -- the new xyz coordinates, dimension (3,ncoords)
                                                                                        
        """
        self.sumb.iteration.groundlevel = 1
        xyz = self.metricConversion*xyz
        self.sumb.updatefacesglobal(xyz,reinitialize)
        self._update_geom_info = True

        return

    def GetGlobalSurfaceCoordinates(self):
        """
        Get the surface coordinates for the mesh    
        """
        if self.sol_type.lower() == 'steady':
            return self.sumb.mddata.mdglobalsurfxx/self.metricConversion
        elif self.sol_type.lower() == 'time spectral':
            return self.sumb.mddata.mdglobalsurfxx[:,:,0]/self.metricConversion
        else:
            print 'invalid solutions type for surface coords...Exiting'
            sys.exit(0)
        #endif
        return

    def warpMesh(self):
        '''
        run the internal meshwarping scheme
        '''
        self.sumb.warpmesh()

        return

    def initializeExternalWarping(self,ndoflocal):
        self.sumb.initializeexternalwarping(ndoflocal)

        return

    def setGrid(self,externaldof):

        self.sumb.setgrid(externaldof)
        self.sumb.setgrid(externaldof)
        self._update_geom_info = True
        return
    
    def getForces(self,cgnsdof):
        self.sumb.getforces1()
        if cgnsdof > 0:
            return self.sumb.getforces2(cgnsdof)
        else:
            return numpy.empty([0],dtype='d')
        # end if

# =============================================================================

class SUmbInterface(object):
    """Represents a SUmb flow solution."""

    def __init__(self,*args,**kwargs):#deforming_mesh=False):
        """Initialize the object.

        Keyword arguments:

        communicator   -- an MPI communicator for specifying which processors
                          to run on (optional)
        deforming_mesh -- set to True or False whether or not this is a 
                          deforming mesh case (optional, default is False)
                          

        """

        if 'sumb' in kwargs:
            self.sumb = kwargs['sumb']
        else:
            sumb_mod = MExt('sumb_parallel')
            self.sumb = sumb_mod._module
        # end if
        
        # The very first thing --> Set the MPI Communicators
        if 'comm' in kwargs:
            self.sumb.communication.sumb_comm_world = kwargs['comm'].py2f()
            self.sumb.communication.sumb_petsc_comm_world = kwargs['comm'].py2f()
            self.sumb.communication.sumb_comm_self  = mpi.COMM_SELF.py2f()
            self.sumb_comm_world = kwargs['comm']
        else:
            self.sumb.communication.sumb_comm_world = mpi.COMM_WORLD.py2f()
            self.sumb.communication.sumb_comm_self  = mpi.COMM_SELF.py2f()
            self.sumb_comm_world = mpi.COMM_WORLD
        # end if
        
        if 'init_petsc' in kwargs:
            if kwargs['init_petsc']:
                self.sumb.initializepetsc()
        # end if

        # Setup the mesh object with sumb_comm_world
        self.Mesh = SUmbMesh(self.sumb_comm_world,self.sumb)

        # Determine the rank sumb_comm_world size
        self.myid = self.sumb.communication.myid = self.sumb_comm_world.rank
        self.nproc = self.sumb.communication.nproc = self.sumb_comm_world.size

        # Allocate the memory for SENDREQUESTS and RECVREQUESTS.
        try:
            self.sumb.communication.sendrequests = numpy.zeros(
        			  (self.sumb_comm_world.size))
            self.sumb.communication.recvrequests = numpy.zeros(
        			  (self.sumb_comm_world.size))
        except:
            print "Memory allocation failure for SENDREQUESTS " \
        	  "and RECVREQUESTS."
            sys.exit(1)
        # end try

        # Set the SUmb module value of standalonemode to false and
        # the value of deforming_grid to the input value.
        self.sumb.iteration.standalonemode = False
        self.sumb.iteration.deforming_grid = deforming_mesh = False

        # Write the intro message
        self.sumb.writeintromessage()

        # Set the frompython flag to true
        self.sumb.killsignals.frompython=True

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

        if kwargs['options']:
            self.probName = kwargs['options']['probName'][1]
            self.OutputDir = kwargs['options']['OutputDir'][1]
        else:
            self.probName = ''
            self.OutputDir = './'
        # end if 

        self.Mesh.sol_type = sol_type

        startfile = self.OutputDir+self.probName+grid_file+'_autogen.input'
        
        self.generateInputFile(aero_problem,sol_type,grid_file,startfile,*args,**kwargs)
        
        # Make sure the parameter file exists
        if not os.path.isfile(startfile):
            print 'Error: Could not find file %s' % startfile
            return None

        # Store the name of the input file
        self.startfile = startfile
        
        self.sumb.inputio.paramfile[0:len(startfile)] = startfile

        # Read the parameter file
        self.sumb.readparamfile()

        # Set printIteration Flag
        self.sumb.inputiteration.printiterations = kwargs['options']['printIterations'][1]

        #This is just to flip the -1 to 1 possibly a memory issue?
        self.sumb.inputio.storeconvinneriter=abs(self.sumb.inputio.storeconvinneriter)

        if(self.myid ==0):print ' -> Partitioning and Reading Grid'
        self.sumb.partitionandreadgrid()

        if(self.myid==0):print ' -> Preprocessing'
        self.sumb.preprocessing()

        if(self.myid==0):print ' -> Initializing flow'
        self.sumb.initflow()

        # Create dictionary of variables we are monitoring
        nmon = self.sumb.monitor.nmon
        self.monnames = {}
        for i in range(nmon):
            self.monnames[string.strip(
        		   self.sumb.monitor.monnames[i].tostring())] = i
        
        # Create dictionary of the family names
        self.sumb.mdgetfamilynames()
        nfamilies = self.sumb.mddata.mdnfamilies
        self.Mesh.families = {}
        for i in range(nfamilies):
            self.Mesh.families[string.strip(self.sumb.mddata.mdfamilynames[i]
        		       .tostring())] = i + 1

        # Create Surface Node list
        if(self.myid==0): print ' -> Creating Surface Node List'
        self.sumb.mdcreatensurfnodes()

        # Reduce the total number of blocks
        self.Mesh.nmeshblocks = self.sumb_comm_world.allreduce(
        			     self.sumb.block.ndom,mpi.SUM)
        #Set flags for ADjoint initialization
        self.adjointInitialized = False
        
        return

    def setInflowAngle(self,aero_problem):
        '''
        Set the alpha and beta fromthe desiggn variables
        '''
        [velDir,liftDir,dragDir]= self.sumb.adjustinflowangleadj((aero_problem._flows.alpha*(pi/180.0)),(aero_problem._flows.beta*(pi/180.0)),aero_problem._flows.liftIndex)
        self.sumb.inputphysics.veldirfreestream = velDir
        self.sumb.inputphysics.liftdirection = liftDir
        self.sumb.inputphysics.dragdirection = dragDir

        #if (self.myid==0):print '-> Alpha...',aero_problem._flows.alpha*(pi/180.0),aero_problem._flows.alpha,velDir,liftDir,dragDir
        #update the flow vars
        self.sumb.updateflow()
        return

    def resetFlow(self):
        '''
        Reset the flow for the complex derivative calculation
        '''

        self.sumb.setuniformflow()

        return
    
    def generateInputFile(self,aero_problem,sol_type,grid_file,startfile,file_type='cgns',*args,**kwargs):
        ''' Code to generate an SUmb Input File on the fly'''
        if (self.myid==0): print ' -> Generating Input File'
        
        #Convert alpha and beta to a freestream vector
        [velDir,liftDir,dragDir]= self.sumb.adjustinflowangleadj((aero_problem._flows.alpha*(pi/180.0)),(aero_problem._flows.beta*(pi/180.0)),aero_problem._flows.liftIndex)
       
        autofile = open(startfile,'w')

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
            autofile.write("                      Restart file: %s_restart.cgns \n"%(grid_file))
            autofile.write("                           Restart: %s\n"%(kwargs['options']['sol_restart'][1]))
            autofile.write("       Check nondimensionalization: yes\n\n")
            autofile.write("                     New grid file: %s%s%s_NewGrid.cgns\n\n"%(self.OutputDir,self.probName,grid_file))
            autofile.write("                     Solution file: %s%s%s_SolSUmb.cgns\n"%(self.OutputDir,self.probName,grid_file))
            autofile.write("             Surface solution file: %s%s%s_SolSUmb_surface.cgns\n"%(self.OutputDir,self.probName,grid_file))
            
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
        #print 'solver options',kwargs['options']['Equation Type'][1]
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Physics Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "                  Equations: %s\n"%(kwargs['options']['Equation Type'][1]))
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

        if kwargs['options']['Equation Type'][1]=='RANS':
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
        try: kwargs['options']['FamilyRot'][1]
        except KeyError:
            Rotating = False
        else:
            Rotating = True
        #endif
        #print 'Rotating',Rotating,kwargs['options']['FamilyRot'][1]
        if Rotating or sol_type=='Time Spectral':
            autofile.write(  "                             Mach: %12.12e\n"%(0))
            autofile.write(  "            Mach for coefficients: %12.12e\n"%(aero_problem._flows.mach))
            autofile.write(  "            Mach for mesh velocity: %12.12e\n"%(aero_problem._flows.mach))
        else:
            autofile.write(  "                             Mach: %12.12e\n"%(aero_problem._flows.mach))
            autofile.write(  "            Mach for coefficients: %12.12e\n"%(aero_problem._flows.mach))
            autofile.write(  "            Mach for mesh velocity: %12.12e\n"%(0))
        #endif
        autofile.write(  "                          # Default is Mach\n")
        autofile.write(  "                         Reynolds: 100000\n")
        autofile.write(  "       Reynolds length (in meter): %12.12e\n"%(aero_problem._refs.cref*kwargs['options']['MetricConversion'][1]))
        #autofile.write(  "   Free stream velocity direction: 1.0 0.05 0.0\n")
        #autofile.write(  "                   Lift direction: -0.05 1.0 0.0\n")
        autofile.write(  "   Free stream velocity direction: %12.12e %12.12e %12.12e\n"%(velDir[0],velDir[1],velDir[2]))
        autofile.write(  "                   Lift direction: %12.12e %12.12e %12.12e\n"%(liftDir[0],liftDir[1],liftDir[2]))
        autofile.write(  "     # Default is normal to free stream without y-component\n")
        autofile.write(  "   Free stream temperature (in K): %12.12e\n"%(kwargs['options']['Reference Temp.'][1]))
        autofile.write(  " Free stream eddy viscosity ratio: 0.01\n")
        autofile.write(  "  Free stream turbulent intensity: 0.001\n")
        autofile.write(  "\n")
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Reference State\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "            Reference pressure (in Pa): %12.12e\n"%(kwargs['options']['Reference Pressure'][1]))
        autofile.write(  "         Reference density (in kg/m^3): %12.12e\n"%(kwargs['options']['Reference Density'][1]))
        autofile.write(  "          Reference temperature (in K): %12.12e\n"%(kwargs['options']['Reference Temp.'][1]))
        autofile.write(  " Conversion factor grid units to meter: %6.4f\n"%(kwargs['options']['MetricConversion'][1]))
        self.Mesh.metricConversion = kwargs['options']['MetricConversion'][1]
        autofile.write( "\n")
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Geometrical Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "           Reference surface: %12.12e\n"%(aero_problem._refs.sref*kwargs['options']['MetricConversion'][1]**2))

        #autofile.write(  "           Reference surface: %2.1f\n"%(1.0))

        autofile.write(  "            Reference length: %12.12e\n"%(aero_problem._refs.cref*kwargs['options']['MetricConversion'][1]))
        autofile.write(  "    Moment reference point x:  %12.12e\n"%(aero_problem._refs.xref*kwargs['options']['MetricConversion'][1]))
        autofile.write(  "    Moment reference point y:  %12.12e\n"%(aero_problem._refs.yref*kwargs['options']['MetricConversion'][1]))
        autofile.write(  "    Moment reference point z:  %12.12e\n"%(aero_problem._refs.zref*kwargs['options']['MetricConversion'][1]))
        autofile.write( "\n")
        
        #! Write the keywords and default values for the discretization
        #! parameters.

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Fine Grid Discretization Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "       Discretization scheme: %s\n"%(kwargs['options']['Discretization'][1]))
        autofile.write(  "             # Possibilities: Central plus scalar dissipation\n")
        autofile.write(  "             #              : Central plus matrix dissipation\n")
        autofile.write(  "             #              : Central plus CUSP dissipation\n")
        autofile.write(  "             #              : Upwind\n")
        autofile.write( "\n")
        
        if kwargs['options']['Equation Type'][1]=='RANS':
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
      
        autofile.write(  "                           Vis2: %10.8f\n"%(kwargs['options']['Dissipation Coefficients'][1][0]))
        autofile.write(  "                           Vis4: %10.8f  # 1/64\n"%(kwargs['options']['Dissipation Coefficients'][1][1]))
        autofile.write(  "Directional dissipation scaling: yes\n")
        autofile.write(  "   Exponent dissipation scaling: %4.2f\n"%(kwargs['options']['Dissipation Scaling Exponent'][1]))
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

        if sol_type=='Time Spectral':
            #! Write the keywords and default values for the time spectral
            #! parameters.

            autofile.write("-------------------------------------------------------------------------------\n")
            autofile.write("     Time Spectral Parameters\n")
            autofile.write("-------------------------------------------------------------------------------\n")
            autofile.write("             Number time intervals spectral: %d\n"%(kwargs['options']['Time Intervals'][1]))
            autofile.write("\n")
            autofile.write("            Write file for unsteady restart: no\n")
            autofile.write("    Time step (in sec) for unsteady restart: 1000\n")
            autofile.write("\n")
            autofile.write("       Write unsteady volume solution files: no\n")
            autofile.write("      Write unsteady surface solution files: no\n")
            autofile.write("          Number of unsteady solution files: 0\n")
            autofile.write("\n")
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

        if kwargs['options']['Equation Type'][1] =='RANS':
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
        autofile.write(  "     Residual averaging smoothing parameter: .5\n")

        autofile.write(  "                 Number of multigrid cycles: 200\n")
        autofile.write(  "   Number of single grid startup iterations: 0\n")
        autofile.write(  "                                 Save every: 10\n")
        autofile.write(  "                         Save surface every: 10\n")
        autofile.write(  "                                 CFL number: %2.1f\n"%(kwargs['options']['CFL'][1]))
        autofile.write( "\n")
        if kwargs['options']['Equation Type'][1] =='RANS':
            autofile.write(  "                       Turbulent relaxation: Explixit\n")
            autofile.write(  "                            # Possibilities: Explicit\n")
            autofile.write(  "                            #              : Implicit\n")
            autofile.write(  "                     Alpha turbulent DD-ADI: 0.8\n")
            autofile.write(  "                      Beta turbulent DD-ADI: -1  # Same as alpha\n")
        #endif
        autofile.write(  "           Relative L2 norm for convergence: %3.2e\n"%(kwargs['options']['L2Convergence'][1]))
        autofile.write( "\n")

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Multigrid Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "      Number of multigrid cycles coarse grid:  -1  # -1 Means same as on fine grid\n")
        autofile.write(  "                      CFL number coarse grid: -1  # -1 Means same as on fine grid\n")

        autofile.write(  "Relative L2 norm for convergence coarse grid: 1.e-1\n")
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
        autofile.write(  "                 Multigrid cycle strategy: %s\n"%(kwargs['options']['MGCycle'][1]))
        autofile.write( "\n")

#
# ADjoint Parameters
#
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "      ADjoint Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "                       solve ADjoint :  %s\n"%(kwargs['options']['solveADjoint'][1]))
        autofile.write(  "                  # Other possibility: no\n")
        autofile.write(  "	                 set monitor   :  %s\n"%(kwargs['options']['set Monitor'][1]))
        autofile.write(  "                  # Other possibility: no\n")
        autofile.write(  "       Use Approximate Preconditioner: %s\n"%(kwargs['options']['Approx PC'][1]))
        autofile.write(  "                  # Other possibility: yes/no\n")
        autofile.write(  " 	            Adjoint solver type: %s\n"%(kwargs['options']['Adjoint solver type'][1]))
        autofile.write(  "                # Other possibilities: BiCGStab\n")
        autofile.write(  "                #                      CG\n")
        autofile.write(  "                #                      GMRES\n")
        autofile.write(  "                #                      FGMRES\n")
        autofile.write(  "        adjoint relative tolerance   : %3.2e\n"%(kwargs['options']['adjoint relative tolerance'][1]))
        autofile.write(  "        adjoint absolute tolerance   : %3.2e\n"%(kwargs['options']['adjoint absolute tolerance'][1]))
        autofile.write(  "        adjoint divergence tolerance : 1e5\n")
        autofile.write(  "        adjoint max iterations       : %d\n"%kwargs['options']['adjoint max iterations'][1])
        autofile.write(  "        adjoint restart iteration    : %d\n"%kwargs['options']['adjoint restart iteration'][1])
        autofile.write(  "        adjoint monitor step         : %d\n"%kwargs['options']['adjoint monitor step'][1])
        autofile.write(  "        dissipation lumping parameter: %d\n"%kwargs['options']['dissipation lumping parameter'][1])
        
        autofile.write(  "                  Preconditioner Side: %s\n"%(kwargs['options']['Preconditioner Side'][1]))
        autofile.write(  "                # Other possibilities: Left\n")
        autofile.write(  "                #                      Right\n")
        autofile.write(  "	             Matrix Ordering   : %s\n"%(kwargs['options']['Matrix Ordering'][1]))
        autofile.write(  "	               #                 ReverseCuthillMckee\n")
        autofile.write(  "	               #                 Natural\n")
        autofile.write(  "	               #                 NestedDissection\n")
        autofile.write(  "                     #                 OnewayDissection\n")
        autofile.write(  "                     #                 QuotientMinimumDegree\n")
        autofile.write(  "          Global Preconditioner Type : %s\n"%(kwargs['options']['Global Preconditioner Type'][1]))
        autofile.write(  "          #                            Jacobi\n")
        autofile.write(  "          #                            Block Jacobi\n")
        autofile.write(  "          #                            Additive Schwartz\n")
        autofile.write(  "                         ASM Overlap : 5\n")
        autofile.write(  "            Local Preconditioner Type: %s\n"%(kwargs['options']['Local Preconditioner Type'][1]))
        autofile.write(  "            #                          ILU\n")
        autofile.write(  "            #                          ICC\n")
        autofile.write(  "            #                          LU\n")
        autofile.write(  "            #                          Cholesky\n")
        autofile.write(  "                    ILU Fill Levels  : %d\n"%kwargs['options']['ILU Fill Levels'][1])
        autofile.write(  "            Jacobi Scale Factor Type : RowMax\n")
        autofile.write(  "            #                          RowMax\n")
        autofile.write(  "            #                          RowSum\n")
        autofile.write(  "            #                          RowAbs\n")

        #Write the options for the Time Spectral stability derivatives to the input file
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "TS Stability Derivative Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
	autofile.write(  "compute TS stability derivatives : %s\n"%(kwargs['options']['TS Stability'][1]))
	autofile.write(  "#Other Possibilities:               no\n")
        if kwargs['options']['TS Stability'][1]== 'yes':
            autofile.write(  "TS Alpha mode: %s\n"%(kwargs['options']['Alpha Mode'][1]))
            autofile.write(  "TS Beta mode: %s\n"%(kwargs['options']['Beta Mode'][1]))
            autofile.write(  "TS p mode: %s\n"%(kwargs['options']['p Mode'][1]))
            autofile.write(  "TS q mode: %s\n"%(kwargs['options']['q Mode'][1]))
            autofile.write(  "TS r mode: %s\n"%(kwargs['options']['r Mode'][1]))
            autofile.write(  "TS Mach number mode: %s\n"%(kwargs['options']['Mach Mode'][1]))
            autofile.write(  "TS Altitude mode: %s\n"%(kwargs['options']['Altitude Mode'][1]))
        #endif

        #! Write the keywords and default values for the parallel, i.e.
        #! load balance parameters.

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Load balancing Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "        Allowable load imbalance: 0.1\n")
        autofile.write(  "   Split blocks for load balance: %s\n"%(kwargs['options']['Allow block splitting'][1]))
        autofile.write( "\n")

##        ! Write the visualization parameters.

##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( ) "     Visualization Parameters"
##        autofile.write( ) "-----------------------------------------&
##                               &--------------------------------------"
##        autofile.write( ) "          pV3 visualization only: no"
##        autofile.write( )

        if sol_type=='Time Spectral':
##        ! Write the keywords and example values for the grid motion.
            
            autofile.write("-------------------------------------------------------------------------------\n")
            autofile.write("     Grid motion Parameters\n")
            autofile.write( "-------------------------------------------------------------------------------\n")

            autofile.write( "     Rotation point body (x,y,z): %12.12e %12.12e %12.12e\n"%(kwargs['options']['rotCenter'][1][0],kwargs['options']['rotCenter'][1][1],kwargs['options']['rotCenter'][1][2]))
            autofile.write("\n" )

            # inputs for TS stability derivatives
            autofile.write( "Degree polynomial Alpha: 0\n")
            autofile.write( "Degree polynomial Beta: 0\n")
            autofile.write( "Degree polynomial Mach: 0\n")
            
            autofile.write( "Polynomial coefficients Alpha: 0.0\n")
            autofile.write( "Polynomial coefficients Beta: 0.0\n")
            autofile.write( "Polynomial coefficients Mach: 0.0\n")

            if kwargs['options']['Alpha Mode'][1] =='yes':
                autofile.write( "Degree fourier Alpha: 1\n")
                autofile.write( "Degree fourier Beta: 0\n")
                autofile.write( "Degree fourier Mach: 0#1\n")
            elif kwargs['options']['Beta Mode'][1] =='yes':
                autofile.write( "Degree fourier Alpha: 0\n")
                autofile.write( "Degree fourier Beta: 1\n")
                autofile.write( "Degree fourier Mach: 0\n")
            elif kwargs['options']['Mach Mode'][1] =='yes':
                autofile.write( "Degree fourier Alpha: 0\n")
                autofile.write( "Degree fourier Beta: 0\n")
                autofile.write( "Degree fourier Mach: 1\n")
            #endif

            if kwargs['options']['Alpha Mode'][1] =='yes':
                autofile.write( "Omega fourier Alpha: %f\n"%(kwargs['options']['Omega fourier'][1]))
                autofile.write( "Omega fourier Beta: 0.0\n")
                autofile.write( "Omega fourier Mach: 0.0\n")
            elif kwargs['options']['Beta Mode'][1] =='yes':
                autofile.write( "Omega fourier Alpha: 0.0\n")
                autofile.write( "Omega fourier Beta: %f\n"%(kwargs['options']['Omega fourier'][1]))
                autofile.write( "Omega fourier Mach: 0.0\n")
            elif kwargs['options']['Mach Mode'][1] =='yes':
                autofile.write( "Omega fourier Alpha: 0.0\n")
                autofile.write( "Omega fourier Beta: 0.0\n")
                autofile.write( "Omega fourier Mach: %f\n"%(kwargs['options']['Omega fourier'][1]))
            #endif

            if kwargs['options']['Alpha Mode'][1] =='yes':
                autofile.write( "Fourier cosine coefficients Alpha: 0.0 0.0\n")
                autofile.write( "Fourier sine coefficients Alpha: %12.12e\n"%(kwargs['options']['Fourier sine coefficient'][1]))
            elif kwargs['options']['Beta Mode'][1] =='yes':
                autofile.write( "Fourier cosine coefficients Beta: 0.0 0.0\n")
                autofile.write( "Fourier sine coefficients Beta: %12.12e\n"%(kwargs['options']['Fourier sine coefficient'][1]))
            elif kwargs['options']['Mach Mode'][1] =='yes':
                autofile.write( "Fourier cosine coefficients Mach: 0.0 0.0\n")
                autofile.write( "Fourier sine coefficients Mach: %12.12e\n"%(kwargs['options']['Fourier sine coefficient'][1]))
            #endif
        #endif
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
        #autofile.write( "                Monitoring variables: resrho_cl_cd_cmx_cmy_cmz\n")
        autofile.write( "                Monitoring variables: resrho_cl_cd_cfx_cfy_cfz\n")
        autofile.write( " Monitor massflow sliding interfaces: no\n")
        autofile.write( "            Surface output variables: rho_cp_vx_vy_vz_mach\n")
        autofile.write( "           Volume output variables: ptloss_resrho\n")
        autofile.write( "\n")

        # The section to overwrite the rotation info for the families.
        try:
            autofile.write( "------------------------------------------------------------------------------\n")
            autofile.write( "     Family rotation info \n")
            autofile.write( "------------------------------------------------------------------------------\n")
            autofile.write( "\n")
            
            autofile.write( "                               Rotation center  Rotation rate (rad/s)\n")
            autofile.write( "Rotating family %s : %6.6f %6.6f %6.6f    %6.6f %6.6f %6.6f    \n"%(kwargs['options']['FamilyRot'][1],kwargs['options']['rotCenter'][1][0],kwargs['options']['rotCenter'][1][1],kwargs['options']['rotCenter'][1][2],kwargs['options']['rotRate'][1][0],kwargs['options']['rotRate'][1][1],kwargs['options']['rotRate'][1][2]))
        except:
            if(self.myid ==0): print ' -> No rotating families Present'
            #endif
        #endtry
        
#        autofile.write( "Rotating family <family_name2> : 0.0 0.0 0.0    1.e+3 0.e+0 0.e+0    \n")
#        autofile.write( ) "Etc."
#        autofile.write( )

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

    def RunIterations(self, *args, **kwargs):
        """Run the flow solver for a given number of time steps and
           multigrid cycles.

        Keyword arguments:

        sol_type  -- Solution type: Steady, Time Spectral or Unsteady

        ncycles -- The number of multigrid cycles to run on the finest mesh,
                   which defaults to the previous value used when calling this
                   method OR the value read from the input file if this method
                   has not been called before

        ntimesteps -- The number of unsteady
                      physical time steps which defaults to the value
                      read from the input file

        """
        
        self.sumb.monitor.niterold = self.sumb_comm_world.bcast(self.sumb.monitor.niterold,root=0)

        try: kwargs['sol_type']
        except:
            print 'sol_type not defined...exiting...'
            sys.exit(0)
        else:
            sol_type = kwargs['sol_type']
        #endtry

        try: kwargs['ncycles']
        except:
            print 'ncycles not defined...exiting...'
            sys.exit(0)
        else:
            ncycles = kwargs['ncycles']
        #endtry

        try: kwargs['storeHistory']
        except:
            storeHistory=False
        else:
            storeHistory= kwargs['storeHistory']
        #endtry

        #reset python failute check to false
        self.sumb.killsignals.routinefailed=False

        if sol_type.lower() in ['steady', 'time spectral']:
           
            #set the number of cycles for this call
            self.sumb.inputiteration.ncycles = ncycles
            
            if (self.sumb.monitor.niterold == 0 and self.sumb.monitor.nitercur == 0 and self.sumb.iteration.itertot == 0):
                # No iterations have been done
                if (self.sumb.inputio.storeconvinneriter):
                    nn = self.sumb.inputiteration.nsgstartup+self.sumb.inputiteration.ncycles
                    if(self.myid==0):
                        # self.sumb.monitor.convarray = None
                       self.sumb.deallocconvarrays()
                       self.sumb.allocconvarrays(nn)
                    #endif
                #endif

            elif(self.sumb.monitor.nitercur == 0 and  self.sumb.iteration.itertot == 0):

                # Reallocate convergence history array and
                # time array with new size, storing old values from restart
                if (self.myid == 0):
                    # number of time steps from restart
                    ntimestepsrestart = self.sumb.monitor.ntimestepsrestart
                                        
                    if (self.sumb.inputio.storeconvinneriter):
                        # number of iterations from restart
                        niterold = self.sumb.monitor.niterold#[0]
                        if storeHistory:
                            # store restart convergence history and deallocate array
                            temp = copy.deepcopy(self.sumb.monitor.convarray[:niterold+1,:])
                            self.sumb.deallocconvarrays()
                            # allocate convergence history array with new extended size
                            self.sumb.allocconvarrays(temp.shape[0]
                                                 +self.sumb.inputiteration.ncycles-1)
                            # recover values from restart and deallocate temporary array
                            self.sumb.monitor.convarray[:temp.shape[0],:temp.shape[1]] = temp
                            temp = None
                        else:
                            temp=copy.deepcopy(self.sumb.monitor.convarray[0,:,:])
                            self.sumb.deallocconvarrays()
                            # allocate convergence history array with new extended size
                            self.sumb.allocconvarrays(self.sumb.inputiteration.ncycles+1+niterold)
                            self.sumb.monitor.convarray[0,:,:] = temp
                            self.sumb.monitor.convarray[1,:,:] = temp
                            temp = None
                        #endif
                    #endif
                #endif
            else:

                # More Time Steps / Iterations in the same session
            
                # Reallocate convergence history array and
                # time array with new size, storing old values from previous runs
                if (self.myid == 0):
                    if (self.sumb.inputio.storeconvinneriter):
                        if storeHistory:
                            # store previous convergence history and deallocate array
                            temp = copy.deepcopy(self.sumb.monitor.convarray)

                            self.sumb.deallocconvarrays()
                            # allocate convergence history array with new extended size
                            nn = self.sumb.inputiteration.nsgstartup+self.sumb.inputiteration.ncycles
                           
                            self.sumb.allocconvarrays(temp.shape[0]+nn-1)
                        
                            # recover values from previous runs and deallocate temporary array
                            self.sumb.monitor.convarray[:temp.shape[0],:] = copy.deepcopy(temp)
                            
                            temp = None
                        else:
                            temp=copy.deepcopy(self.sumb.monitor.convarray[0,:,:])
                            self.sumb.deallocconvarrays()
                            self.sumb.allocconvarrays(self.sumb.inputiteration.ncycles+1)
                            self.sumb.monitor.convarray[0,:,:] = temp
                            temp = None
                        #endif
                    #endif
                #endif
                # re-initialize iteration variables
                self.sumb.inputiteration.mgstartlevel = 1
                if storeHistory:
                    self.sumb.monitor.niterold  = self.sumb.monitor.nitercur
                else:
                    self.sumb.monitor.niterold  = 1#0#self.sumb.monitor.nitercur
                #endif
                self.sumb.monitor.nitercur  = 0#1
                self.sumb.iteration.itertot = 0#1

               
                # update number of time steps from restart
                self.sumb.monitor.ntimestepsrestart = self.sumb.monitor.ntimestepsrestart \
                                                 + self.sumb.monitor.timestepunsteady
                # re-initialize number of time steps previously run (excluding restart)             
                self.sumb.monitor.timestepunsteady = 0
                # update time previously run
                self.sumb.monitor.timeunsteadyrestart = self.sumb.monitor.timeunsteadyrestart \
                                                   + self.sumb.monitor.timeunsteady
                # re-initialize time run
                self.sumb.monitor.timeunsteady = 0.0
                
            #endif
        elif sol_type.lower()=='unsteady':
            print 'unsteady not implemented yet...'
            sys.exit(0)
        else:
            print 'unknown sol_type:',sol_type,'exiting...'
            sys.exit(0)
        #endif

    
        self.GetMesh()._UpdateGeometryInfo()
      
        self.routineFailed = self.sumb_comm_world.allreduce(self.sumb.killsignals.routinefailed,mpi.MIN)
       
        if (abs(self.routineFailed)==True):
            if self.myid ==0:print 'Error raise in updateGeometry'
            raise ValueError
        #endif
     
        self.sumb.inputiteration.l2convrel = kwargs['options']['L2ConvergenceRel'][1]
                
        # Now coll the solver

        self.sumb.solver()

        #Check to see whether we have a valid solution
        #in this case routineFailed will be triggered on all processors
        self.routineFailed = self.sumb_comm_world.allreduce(abs(self.sumb.killsignals.routinefailed),mpi.SUM)/mpi.COMM_WORLD.size
        #therefore sum up and devide by nProc
        if (abs(self.routineFailed)==True):
            if self.myid ==0: print 'Error Raised in solver'
            raise ValueError
        #endif

#         if kwargs['options']['printIterations'][1] == False:
#             # If we weren't printing iterations, output the 0th (initial) 1st (start of this set of iterations) and final
#             if self.myid == 0:
#                 print '  -> CFD Initial Rho Norm: %10.5e'%(self.sumb.monitor.convarray[0,0,0])
#                 print '  -> CFD Start   Rho Norm: %10.5e'%(self.sumb.monitor.convarray[1,0,0])
#                 print '  -> CFD Final   Rho Norm: %10.5e'%(self.sumb.monitor.convarray[self.sumb.monitor.nitercur,0,0])

        return

    def getResiduals(self):
        if self.myid == 0:
            return self.sumb.monitor.convarray[0,0,0],self.sumb.monitor.convarray[1,0,0],self.sumb.monitor.convarray[self.sumb.monitor.nitercur,0,0]
        else:
            return None,None,None
        # end if

        

 #        sys.exit(0)
        
#         if (sumb.monitor.niterold == 0 and
#             sumb.monitor.nitercur == 0 and
#             sumb.iteration.itertot == 0):
#             # No iterations have been done
#             if (ncycles):
#                 # Set new value of unsteady physical time steps
#                 if (ncycles[0] != sumb.inputunsteady.ntimestepsfine):
#                     sumb.inputunsteady.ntimestepsfine = ncycles[0]
#                 # Set new value of MG cycles
#                 if (len(ncycles) > 1):
#                   if (ncycles[1] != sumb.inputiteration.ncycles):
#                       if(self.myid==0):print 'setting ncycles'
#                       sumb.inputiteration.ncycles = ncycles[1]
#                 # Reallocate convergence history array and
#                 # time array with new size
#                 if (self.myid == 0):
#                     print 'sumb.monitor.timearray'#,sumb.monitor.timearray 
#                     sumb.monitor.timearray = None
#                     print 'sumb.monitor.timearray',sumb.monitor.timearray 
#                     sumb.monitor.timedataarray = None
#                     sumb.alloctimearrays(sumb.inputunsteady.ntimestepsfine)
#                     if (sumb.inputio.storeconvinneriter):
#                         nn = sumb.inputiteration.nsgstartup+sumb.inputiteration.ncycles
#                         if (sumb.inputphysics.equationmode==2):#2 is unsteady
#                             nn=sumb.inputunsteady.ntimestepsfine*sumb.inputiteration.ncycles
#                         #endif
#                         sumb.monitor.convarray = None
#                         #print 'nn',nn
#                         sumb.allocconvarrays(nn)
#                         #print 'convarray',sumb.monitor.convarray
                        

#         elif (sumb.monitor.nitercur == 0 and
#               sumb.iteration.itertot == 0):

#             # Read in a restart file but no new iterations
#             if (ncycles):
#                 # Set new value of unsteady physical time steps
#                 if (ncycles[0] != sumb.inputunsteady.ntimestepsfine):
#                     sumb.inputunsteady.ntimestepsfine = ncycles[0]
#                 # Set new value of MG cycles
#                 if (len(ncycles) > 1):
#                     if (ncycles[1] != sumb.inputiteration.ncycles):
#                         sumb.inputiteration.ncycles = ncycles[1]
#                 # Reallocate convergence history array and
#                 # time array with new size, storing old values from restart
#                 if (self.myid == 0):
#                     # number of time steps from restart
#                     ntimestepsrestart = sumb.monitor.ntimestepsrestart#[0]
#                     # store restart time history and deallocate arrays
#                     temp_t = copy.deepcopy(sumb.monitor.timearray[:ntimestepsrestart])
#                     sumb.monitor.timearray = None
#                     temp_td = copy.deepcopy(sumb.monitor.timedataarray[:ntimestepsrestart,:])
#                     sumb.monitor.timedataarray = None
#                     # allocate time history arrays with new extended size
#                     sumb.alloctimearrays(ntimestepsrestart+sumb.inputunsteady.ntimestepsfine)
#                     # recover values from restart and deallocate temporary arrays
#                     sumb.monitor.timearray[:temp_td.shape[0]] = temp_t
#                     sumb.monitor.timedataarray[:temp_td.shape[0],:temp_td.shape[1]] = temp_td
#                     temp_t = None
#                     temp_td = None

#                     if (sumb.inputio.storeconvinneriter):
#                         # number of iterations from restart
#                         niterold = sumb.monitor.niterold[0]
#                         # store restart convergence history and deallocate array
#                         temp = copy.deepcopy(sumb.monitor.convarray[:niterold+1,:])
#                         sumb.monitor.convarray = None
#                         # allocate convergence history array with new extended size
#                         sumb.allocconvarrays(temp.shape[0]
#                                              +sumb.inputunsteady.ntimestepsfine
#                                              *sumb.inputiteration.ncycles-1)
#                         # recover values from restart and deallocate temporary array
#                         sumb.monitor.convarray[:temp.shape[0],:temp.shape[1]] = temp
#                         temp = None
                  
#         else:
#             # More Time Steps / Iterations in the same session
#             #print 'more cycles',ncycles[0],ncycles[1],sumb.monitor.convarray,self.myid
            
#             if (ncycles):
#                 # Set new value of unsteady physical time steps to run
#                 if (ncycles[0] != sumb.inputunsteady.ntimestepsfine):
#                     sumb.inputunsteady.ntimestepsfine = ncycles[0]
#                 # Set new value of MG cycles ro run
#                 if (len(ncycles) > 1):
#                     if (ncycles[1] != sumb.inputiteration.ncycles):
#                         # print 'sumb cycles',sumb.inputiteration.ncycles
#                         sumb.inputiteration.ncycles = ncycles[1]
#                         #print sumb.inputiteration.ncycles

#             # Reallocate convergence history array and
#             # time array with new size, storing old values from previous runs
#             if (self.myid == 0):
#                 #print 'more cycles',ncycles[0],ncycles[1],sumb.monitor.convarray
                
#                 #print 'equation mode',sumb.inputphysics.equationmode
#                 if (sumb.inputphysics.equationmode==2):#2 is unsteady
#                     #print 'unsteady'
#                     #print 'checking time array'
#                     # store previous time history and deallocate arrays
#                     temp_t = copy.deepcopy(sumb.monitor.timearray)
#                     #print 'deallocating time array',sumb.monitor.timearray
#                     sumb.monitor.timearray = None
#                     #print 'checking time data array',sumb.monitor.timedataarray
#                     temp_td = copy.deepcopy(sumb.monitor.timedataarray)
#                     #print 'deallocating timedata array'
#                     sumb.monitor.timedataarray = None
#                     #print 'allocating new time array'
#                     # allocate time history arrays with new extended size
#                     sumb.alloctimearrays(temp_td.shape[0]+sumb.inputunsteady.ntimestepsfine)
#                     # recover values from previous runs and deallocate temporary arrays
#                     sumb.monitor.timearray[:temp_td.shape[0]] = temp_t
#                     sumb.monitor.timedataarray[:temp_td.shape[0],:temp_td.shape[1]] = temp_td
#                     temp_t = None
#                     temp_td = None

#                 if (sumb.inputio.storeconvinneriter):
#                     #print 'store conv?',sumb.inputio.storeconvinneriter,sumb.inputiteration.ncycles-1
#                     #print 'conv, arreay',sumb.monitor.convarray[:,0,1]
#                     #sdfg
#                     # store previous convergence history and deallocate array
#                     temp = copy.deepcopy(sumb.monitor.convarray)
#                     #temp = copy.deepcopy(sumb.monitor.convarray)
#                     #print 'conv, arreay',sumb.monitor.convarray
#                     sumb.monitor.convarray = None
#                     #print 'allocating convergence arrays for new size',sumb.monitor.convarray,temp.shape[0],sumb.inputunsteady.ntimestepsfine,sumb.inputiteration.ncycles,sumb.inputiteration.ncycles-1
                    
#                     #print 'testing',temp.shape[0]
#                     # allocate convergence history array with new extended size
#                     nn = sumb.inputiteration.nsgstartup+sumb.inputiteration.ncycles
#                     if (sumb.inputphysics.equationmode==2):#2 is unsteady
#                         nn=sumb.inputunsteady.ntimestepsfine*sumb.inputiteration.ncycles
#                     #endif
#                     sumb.allocconvarrays(temp.shape[0]+nn-1)
#                     #print 'convergence shape',sumb.monitor.convarray.shape,temp.shape, sumb.monitor.convarray[:temp.shape[0],:].shape,temp[:,0,:].shape
#                     # recover values from previous runs and deallocate temporary array
#                     sumb.monitor.convarray[:temp.shape[0],:] = copy.deepcopy(temp)
# ##                     # allocate convergence history array with new extended size
# ##                     sumb.allocconvarrays(temp.shape[0]
# ##                                          +sumb.inputunsteady.ntimestepsfine
# ##                                          *sumb.inputiteration.ncycles-1)
# ##                     # recover values from previous runs and deallocate temporary array
# ##                     sumb.monitor.convarray[:temp.shape[0],:temp.shape[1]] = temp
#                     temp = None

#             # re-initialize iteration variables
# 	    sumb.inputiteration.mgstartlevel = 1
#             sumb.monitor.niterold  = sumb.monitor.nitercur
#             sumb.monitor.nitercur  = 0
#             sumb.iteration.itertot = 0
#             # update number of time steps from restart
#             sumb.monitor.ntimestepsrestart = sumb.monitor.ntimestepsrestart \
#                                            + sumb.monitor.timestepunsteady
#             # re-initialize number of time steps previously run (excluding restart)             
#             sumb.monitor.timestepunsteady = 0
#             # update time previously run
#             sumb.monitor.timeunsteadyrestart = sumb.monitor.timeunsteadyrestart \
#                                              + sumb.monitor.timeunsteady
#             # re-initialize time run
#             sumb.monitor.timeunsteady = 0.0

#         #endif
#         if self.myid ==0:print 'setupfinished'
#         #self.Mesh.WriteMeshFile('newmesh.cgns')
#         if self.myid ==0: print 'calling solver'
#         self.GetMesh()._UpdateGeometryInfo()

#         if sumb.killsignals.routinefailed==True:
#             raise ValueError
#             return    
#         #endif
        
#         #print 'mesh info updated'
#         #self.Mesh.WriteMeshFile('newmesh0.cgns')
#         #print 'new file written'
#         sumb.solver()
#         if self.myid ==0:print 'solver called'
        
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
        
    def WriteVolumeSolutionFile(self,filename=None,gridname=None):
        """Write the current state of the volume flow solution to a CGNS file.
 
        Keyword arguments:
        
        filename -- the name of the file (optional)

        """

        self.sumb.monitor.writegrid = False
        if (filename):
            self.sumb.inputio.solfile[:] = ''
            self.sumb.inputio.solfile[0:len(filename)] = filename
        if (gridname):
            self.sumb.inputio.newgridfile[:] = ''
            self.sumb.inputio.newgridfile[0:len(gridname)] = gridname
            self.sumb.monitor.writegrid = True
#         if(self.sumb.iteration.changing_grid or self.sumb.inputmotion.gridmotionspecified):
#             self.sumb.monitor.writegrid=True
# 	else:
#             self.sumb.monitor.writegrid=False


        self.sumb.monitor.writevolume=True
        self.sumb.monitor.writesurface=False
        self.sumb.writesol()

    def WriteSurfaceSolutionFile(self,*filename):
        """Write the current state of the surface flow solution to a CGNS file.
 
        Keyword arguments:
        
        filename -- the name of the file (optional)

        """
        if (filename):
            self.sumb.inputio.surfacesolfile[:] = ''
            self.sumb.inputio.surfacesolfile[0:len(filename[0])] = filename[0]
        self.sumb.monitor.writegrid=False
        self.sumb.monitor.writevolume=False
        self.sumb.monitor.writesurface=True
        self.sumb.writesol()

 #    def GetSurfaceLoads(self, family=None, sps=1):
#         """Return an array of the surface forces.
         
#         Keyword arguments:
         
#         family -- optional string specifying the return of forces
#                   only for the specified family.
#         sps    -- spectral time step (optional, default is set to 1).
#                   (sps=1 for usual steady or unsteady models)
                                   
#         """
# #        if(sumb.flovarrefstate.viscous):
# #            sumb.iteration.rfil = 1.
# #            for i in range(sumb.block.ndom):
# #                sumb.setpointers(i,1)
# #                sumb.viscousflux() 
#         if(family):
#             try:
#                 index = self.Mesh.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdcreatesurfforcelist(sps,index)
#             return sumb.mddata.mdsurfforce[:,start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.Mesh.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdcreatesurfforcelist(sps,0)
#             else:
#                 for n in range(nfamilies):
#                     [start_ind,end_ind] = sumb.mdcreatesurfforcelist(sps,n+1)
#             return sumb.mddata.mdsurfforce

#     def DeallocateSurfaceLoads(self):
#         """Deallocate memory used for the surface loads."""
#         sumb.mddeletesurfforcelist()

#     def GetSurfaceLoadsLocal(self, family=None, sps=1):
#         """Return an array of the surface forces.
         
#         Keyword arguments:
         
#         family -- optional string specifying the return of forces
#                   only for the specified family.
#         sps    -- spectral time step (optional, default is set to 1).
#                   (sps=1 for usual steady or unsteady models)
                                   
#         """
# #        if(sumb.flovarrefstate.viscous):
# #            sumb.iteration.rfil = 1.
# #            for i in range(sumb.block.ndom):
# #                sumb.setpointers(i,1)
# #                sumb.viscousflux() 
#         if(family):
#             try:
#                 index = self.Mesh.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdcreatesurfforcelistlocal(sps,index)
#             return sumb.mddatalocal.mdsurfforcelocal[:,start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.Mesh.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdcreatesurfforcelistlocal(sps,0)
#             else:
#                 for n in range(nfamilies):
#                     [start_ind,end_ind] = sumb.mdcreatesurfforcelistlocal(sps,n+1)
#             return sumb.mddatalocal.mdsurfforcelocal

#     def DeallocateSurfaceLoadsLocal(self):
#         """Deallocate memory used for the surface loads."""
#         sumb.mddeletesurfforcelistlocal()

#     def AccumulateLoads(self,oml_loads_local):
#         '''
#         Sum up all of the local loads to get the total force on the oml
#         '''
#         oml_loads = self.sumb_comm_world.Allreduce(oml_loads_local,mpi.SUM)

#         return oml_loads


#     def GetSurfaceCp(self, family=None, sps=1):
#         """Return an array of the surface pressure coefficients.

#         Keyword arguments:

#         family -- optional string specifying the return of surface Cp
#                   only for the specified family.
#         sps    -- spectral time step (optional, default is set to 1).
#                   (sps=1 for usual steady or unsteady models)

#         """
#         if(family):
#             try:
#                 index = self.Mesh.families[family]
#             except KeyError:
#                 print "Error: No such family '%s'" % family
#                 return None
#             [start_ind,end_ind] = sumb.mdcreatesurfvarlist(sps,index)
#             return sumb.mddata.mdsurfval[start_ind-1:end_ind]
#         else:
#             nfamilies = len(self.Mesh.families)
#             if (nfamilies == 0):
#                 [start_ind,end_ind] = sumb.mdcreatesurfvarlist(sps,0)
#             else:
#                 for n in range(nfamilies):
#                     [start_ind,end_ind] = sumb.mdcreatesurfvarlist(sps,n+1)
#             return sumb.mddata.mdsurfval

#     def DeallocateSurfaceCp(self):
#         """Deallocate memory used for the surface pressure coefficients."""
#         sumb.mddeletesurfvallist()

    def GetMach(self):
        """Get the current freestream Mach number."""
        return self.sumb.inputparam.mach

    def SetMach(self,mach):
        """Set the freestream Mach number."""
        self.sumb.inputparam.mach = mach

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
            if (self.sumb.monitor.niterold == 0 and
                self.sumb.monitor.nitercur == 0 and
                self.sumb.iteration.itertot == 0):
	        history = None
            elif (self.sumb.monitor.nitercur == 0 and
                  self.sumb.iteration.itertot == 0):
	        niterold = self.sumb.monitor.niterold[0]	    
                history = self.sumb.monitor.convarray[:niterold+1,index]
            else:
	        history = self.sumb.monitor.convarray[:,index]
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
	        
        #Set the mesh level and timespectral instance for this
        #computation
        self.level = 1
        self.sps = 1
        
        self.sumb.iteration.currentlevel=1
        self.sumb.iteration.groundlevel=1

        #Check to see if initialization has already been performed
        if(self.adjointInitialized):
            return
        
        #Run the preprocessing routine. Sets the node numbering and
        #allocates memory.
        if(self.myid==0): print 'preprocessing adjoint'
        self.sumb.preprocessingadjoint(self.level)
        
        #Initialize the design variable and function storage
        if(self.myid==0):print 'Before design init'
        self.sumb.designinit()

        #initalize PETSc
        if(self.myid==0):print 'before petsc'
        self.sumb.initializepetsc()

        #create the neccesary PETSc objects
        if(self.myid==0):print 'before createpetsecars'
        self.sumb.createpetscvars()

        #mark the ADjoint as initialized
        self.adjointInitialized = True
        if(self.myid==0):
            print 'ADjoint Initialized Succesfully...'
        #endif
        if(self.myid==0):print 'before nspatial..   '
        self.nSpatial = self.sumb.adjointvars.ndesignspatial
        if(self.myid==0):print 'returning....'
        return

    def setupADjointMatrix(self):
        '''
        Setup the ADjoint dRdw matrix and create the PETSc
        Solution KSP object
        '''
        #create the neccesary PETSc objects
        if(self.myid==0):print 'before createpetscvars'
        #self.sumb.createpetscvars()
        #self.sumb.createpetscmat()
        #self.sumb.setupadjointmatrix(self.level)
        #self.sumb.setupadjointmatrixtranspose(self.level)
        self.sumb.setupallresidualmatrices(self.level)
        self.sumb.setuppetscksp(self.level)

        return

    def releaseAdjointMemeory(self):
        '''
        release the KSP memory...
        '''
        #self.sumb.destroypetscksp()
        #self.sumb.destroypetscvars()
        #self.sumb.destroypetscmat()
        return

    def setupADjointRHS(self,objective):
        '''
        setup the RHS vector for a given cost function
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        if self.myid==0:
            print 'Analyzing ADjoint for costfuntion: ',objective
            print 'SUmb index:',SUmbCostfunctions[objective]
    
        #print SUmbCostfunctions[objective]

        self.sumb.setupadjointrhs(self.level,SUmbCostfunctions[objective])

        return

    def solveADjointPETSc(self):
        '''
        Solve the ADjoint system using PETSc
        '''
        #self.sumb.solveadjointpetsc()
        self.sumb.solveadjointtransposepetsc()

        return

    def verifyPartials(self):
        '''
        Run solverADjoint to verify the partial derivatives in the ADjoint
        '''
        print 'in interface verify partials'
        #self.sumb.verifydrdxsfile()
        self.sumb.verifydcfdx(1)
        #self.sumb.solveradjoint()

        return
    
    

    def setupGradientRHSVolume(self,objective):
        '''
        setup the rhs partial for the mesh derivatives
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing RHS Partial for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        self.sumb.setupgradientrhsvolume(self.level,SUmbCostfunctions[possibleObjectives[objective]])
        #self.sumb.setupgradientrhsspatial(self.level,SUmbCostfunctions[possibleObjectives[objective]],self.sps)
        #endfor

        return

    
    def setupGradientRHSFlow(self,objective):
        '''
        setup the rhs partial for the mesh derivatives
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing RHS Partial for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif

        self.sumb.setupgradientrhsflow(self.level,SUmbCostfunctions[possibleObjectives[objective]])
        #self.sumb.setupgradientrhsextra(self.level,SUmbCostfunctions[possibleObjectives[objective]],self.sps)
        #endfor

        return
    
        
    def setupGradientMatrixVolume(self):

        """Set up the residual sensitivity w.r.t. spatial design variables."""

        self.sumb.setupgradientmatrixspatial(self.level)

	if (self.myid == 0):
            print "ADjoint: Spatial residual sensitivity set up successfully."

        return

    def setupGradientMatrixFlow(self):

        """Set up the residual sensitivity w.r.t. spatial design variables."""

        self.sumb.setupgradientmatrixextra(self.level)

	if (self.myid == 0):
            print "ADjoint: Extra Vars residual sensitivity set up successfully."

        return

    def setupVolumeSurfaceDerivatives(self):
        """Set up the derivative of the volume mesh  w.r.t. the
        CFD Surface...(meshwarpingderivatives)"""
        
        self.sumb.setupvolumesurfacederivativesdv()
        
	if (self.myid == 0):
            print "Meshwarping derivatives set up successfully."
        #endif
        return
    
    def computeTotalSurfaceDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing total Surface mesh derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        self.sumb.computeadjointgradientsurfacedv(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return
    
    def computeTotalVolumeDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing total mesh derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        self.sumb.computeadjointgradientspatial(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return
        
    def getTotalVolumeDerivatives(self,objective):

        """
        Get the  sensitivities from Self.Sumb.
                
	"""
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        grad = numpy.zeros((self.sumb.adjointvars.ndesignspatial),float)
        grad[:] = self.sumb.adjointvars.functiongradspatial[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for item in objective:
        
        #for i in xrange(self.sumb.adjointvars.ndesignspatial):
            #grad[i] = self.sumb.adjointvars.functiongradspatial[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
           
        #endfor
        #endfor
        
        return grad
    def getTotalSurfaceDerivatives(self,objective):

        """
        Get the  sensitivities from Self.Sumb.
                
	"""
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        grad = numpy.zeros((3*self.sumb.mddata.mdnsurfnodescompact),float)
        grad[:] = self.sumb.adjointvars.functiongradsurfacedv[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
   
        
        return grad

    def computeTotalFlowDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing total flow derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        self.sumb.computeadjointgradientextra(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return

    def getTotalFlowDerivatives(self,objective):

        """
        Get the  sensitivities from Self.Sumb.
                
	"""
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        grad = numpy.zeros((self.sumb.adjointvars.ndesignextra),float)
        grad[:] = self.sumb.adjointvars.functiongrad[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for item in objective:
        #for i in xrange(self.sumb.adjointvars.ndesignextra):
            #grad[i] = self.sumb.adjointvars.functiongrad[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
        #endfor
        #endfor
        
        return grad

    def computeAeroCouplingDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        #for item in objective:
        if self.myid==0:
            print 'Computing Aero Coupling derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
            print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
        #endif
        
        self.sumb.computeaerocoupling(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

        return
        
    def getAeroCouplingDerivatives(self,objective):

        """
        Get the  sensitivities from Self.Sumb.
                
	"""
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        grad = numpy.zeros((self.sumb.adjointvars.ndesignspatial),float)
        #for item in objective:
        grad[:] = self.sumb.adjointvars.functiongradcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for i in xrange(self.sumb.adjointvars.ndesignspatial):
            #grad[i] = self.sumb.adjointvars.functiongradcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
            
        #endfor
        #endfor
        
        return grad

    def computeAeroExplicitCouplingDerivative(self,objective):
        '''
        compute the total mesh derivatives
        '''
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        try:
            #for item in objective:
            if self.myid==0:
                print 'Computing Aero Coupling derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
                print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
            #endif

            self.sumb.computeaeroexpcoupling(SUmbCostfunctions[possibleObjectives[objective]])
        except:
            print 'not an aerodynamic cost function'
        #end

        return
        
    def getAeroExplicitCouplingDerivatives(self,objective):

        """
        Get the  sensitivities from Self.Sumb.
                
	"""
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }

        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        grad = numpy.zeros((self.sumb.adjointvars.ndesignspatial),float)
        #for item in objective:
        try:
            grad[:] = self.sumb.adjointvars.functiongradexpcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
            #for i in xrange(self.sumb.adjointvars.ndesignspatial):
                #grad[i] = self.sumb.adjointvars.functiongradexpcoupling[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
            #endfor
        except:
            print 'not an aerodynamic cost function'
        #endtry
        
        return grad
    
    def getGlobalNodesLocal(self,blocknum,il,jl,kl):
        #get global node ordering from sumb
        if(self.myid==0):print 'in sumbInterface'
        globalNodes = self.sumb.getglobalnodes(blocknum,il,jl,kl)
        
        return globalNodes

    def getFunctionValues(self):
        '''
        retrieve the solution values from SUmb
        '''
        print 'interface getting solution'
        # Map cost functions
        self.sumb.getsolution()

        #print 'solution mapped'
        print 'sumb value:',self.sumb.adjointvars.costfuncliftcoef-1

        SUmbsolutions = {'cl':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncliftcoef-1],\
                         'cd':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncdragcoef-1],\
                         'cFx':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncforcexcoef-1],\
                         'cFy':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncforceycoef-1],\
                         'cFz':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncforcezcoef-1],\
                         'cMx':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncmomxcoef-1],\
                         'cMy':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncmomycoef-1],\
                         'cMz':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncmomzcoef-1],\
                         'cMzAlpha':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccmzalpha-1],\
                         'cM0':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccm0-1],\
                         'clAlpha':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfuncclalpha-1],\
                         'cl0':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccl0-1],\
                         'cdAlpha':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccdalpha-1],\
                         'cd0':self.sumb.adjointvars.functionvalue[self.sumb.adjointvars.costfunccd0-1]
                         }

        return SUmbsolutions

    def augmentADjointRHS(self,objective,structAdjoint):
        '''
        run the routines to augment the RHS of the ADjoint
        '''

        self.sumb.setupcouplingmatrixstruct(1)
        self.sumb.setupadjointrhsstruct(structAdjoint)
        

        return

    def getAdjoint(self,objective):

        """
        Get the  sensitivities from Self.Sumb.
                
	"""
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }

        self.sumb.getadjoint(SUmbCostfunctions[possibleObjectives[objective]])
        nstate = self.sumb.flowvarrefstate.nw*self.sumb.adjointvars.ncellsglobal
        adjoint = numpy.zeros((nstate),float)
        #for item in objective:
        #This should work
        adjoint[:] = self.sumb.adjointvars.adjoint[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for i in xrange(nstate):
            # adjoint[i] = self.sumb.adjointvars.adjoint[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
        #endfor
       
        return adjoint

    def getTotalStructDerivatives(self,objective):

        """
        Get the  sensitivities from Self.Sumb.
                
	"""
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        
        grad = numpy.zeros((self.sumb.adjointvars.ndesignspatial),float)
        #for item in objective:
        grad[:] = self.sumb.adjointvars.functiongradstruct[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for i in xrange(self.sumb.adjointvars.ndesignspatial):
            #grad[i] = self.sumb.adjointvars.functiongradstruct[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
        #endfor
        
        return grad
    
    def aeroComputeTotalDerivatveStruct(self,objective,structAdjoint={}):
        '''
        compute the force portion of the total structural derivative
        based on the coupled structural adjoint
        '''
        
        self.sumb.setupcouplingtotalstruct(1)
        
        SUmbCostfunctions = {'cl':self.sumb.adjointvars.costfuncliftcoef,\
                             'cd':self.sumb.adjointvars.costfuncdragcoef,\
                             'cFx':self.sumb.adjointvars.costfuncforcexcoef,\
                             'cFy':self.sumb.adjointvars.costfuncforceycoef,\
                             'cFz':self.sumb.adjointvars.costfuncforcezcoef,\
                             'cMx':self.sumb.adjointvars.costfuncmomxcoef,\
                             'cMy':self.sumb.adjointvars.costfuncmomycoef,\
                             'cMz':self.sumb.adjointvars.costfuncmomzcoef,\
                             'cMzAlpha':self.sumb.adjointvars.costfunccmzalpha,\
                             'cM0':self.sumb.adjointvars.costfunccm0,\
                             'clAlpha':self.sumb.adjointvars.costfuncclalpha,\
                             'cl0':self.sumb.adjointvars.costfunccl0,\
                             'cdAlpha':self.sumb.adjointvars.costfunccdalpha,\
                             'cd0':self.sumb.adjointvars.costfunccd0
                             }
        
        possibleObjectives = { 'lift':'cl','Lift':'cl','CL':'cl','cl':'cl',\
                               'drag':'cd','Drag':'cd','CD':'cd','cd':'cd',\
                               'forcx':'cFx','xForce':'cFx','CFX':'cFx','cFx':'cFx',\
                               'forcey':'cFy','yForce':'cFy','CFY':'cFy','cFy':'cFy',\
                               'forcez':'cFz','zForce':'cFz','CFZ':'cFz','cFz':'cFz',\
                               'momentx':'cMx','xMoment':'cMx','CMX':'cMx','cMx':'cMx',\
                               'momenty':'cMy','yMoment':'cMy','CMY':'cMy','cMy':'cMy',\
                               'momentz':'cMz','zMoment':'cMz','CMZ':'cMz','cMz':'cMz',\
                               'cMzAlpha':'cMzAlpha',\
                               'cM0':'cM0',\
                               'clAlpha':'clAlpha',\
                               'cl0':'cl0',\
                               'cdAlpha':'cdAlpha',\
                               'cd0':'cd0'
                               }
        try:
            #for item in objective:
            if self.myid==0:
                print 'Computing Aero Coupling derivative for costfuntion: ',possibleObjectives[objective]#item#objective[item]
                print 'SUmb index:',SUmbCostfunctions[possibleObjectives[objective]]
            #endif
                
            self.sumb.setupadjointtotalstruct(structAdjoint,SUmbCostfunctions[possibleObjectives[objective]])
        except:
            print 'not an aerodynamic cost function'
        #end

        return

    def computeStabilityParameters(self):
        '''
        run the stability derivative driver to compute the stability parameters
        from the time spectral solution
        '''

        self.sumb.stabilityderivativedriver()

        return
    


