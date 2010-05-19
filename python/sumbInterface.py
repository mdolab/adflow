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
#    import mpi
    _parallel = True


    if mpi.COMM_WORLD.rank==0:
        print 'importingMPI',_parallel
    #endif

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

            if mpi.COMM_WORLD.rank==0:
                print "Warning: Running in an MPI environment, but failed"
                print "         to import parallel version of SUmb.  Proceeding"
                print "         with a sequential version."
            #endif

        except ImportError:
            if mpi.COMM_WORLD.rank==0:
                print "Error: Failed to import parallel or sequential version"
                print "       of SUmb."
            #endif
else:
    try:
        import sumb
    except ImportError:
        print "Error: Failed to import sequential version of SUmb."


# =============================================================================

class SUmbMesh(object):
    """Represents a SUmb mesh.
 
    """
 
    def __init__(self,comm):
        """Initialize the object."""
        self._update_geom_info = False
        self._suggar_interface_initialized = False
        self.solid_warp_initialized=False
        self.comm = comm
        self.myid = comm.rank

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
        sumb.iteration.groundlevel = 1
        if (filename):
            sumb.inputio.newgridfile[:] = ''
            sumb.inputio.newgridfile[0:len(filename[0])] = filename[0]
        sumb.monitor.writegrid=True
        sumb.monitor.writevolume=True#False
        sumb.monitor.writesurface=True
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
        #print 'setting block coords',il,jl,kl,xyz.real.shape
        sumb.setblockcoords(blocknum,il,jl,kl,xyz.real)
        #print 'bloock coords set'
        #sumb.mdsetcoor(sps,blocknums,ranges,xyz)
        self._update_geom_info = True

    def _UpdateGeometryInfo(self):
        """Update the SUmb internal geometry info, if necessary."""
        if (self._update_geom_info):
            if (self.myid==0): print 'Updating Geometry...'
            #print 'updatecoords'
            sumb.updatecoordinatesalllevels()
            #print 'update_wall'
            sumb.updatewalldistancealllevels()
            #print 'update sliding'
            sumb.updateslidingalllevels()
            #print 'update metrics'
            sumb.updatemetricsalllevels()
            sumb.updategridvelocitiesalllevels()
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

    def getSingleState(self,blocknum,i,j,k,l):
        '''
        get the requested state value
        '''
        return sumb.getsinglestate(blocknum,i,j,k,l)

    def setSingleState(self,blocknum,i,j,k,l,state):
        '''
        set the requested state value
        '''
        sumb.setsinglestate(blocknum,i,j,k,l,state)
        
        return 
##     def internalWarping(self,new_cfd_surf,indices):
##         '''
##         Call the integrated warping routines...
##         '''
##         sumb.integratedwarp(new_cfd_surf,indices)

##         return

    def initializeInternalWarping(self):
        '''
        Initialize the required variables for the internal
        Meshwarping and derivatives
        '''
        if(sumb.cgnsgrid.cgnsnfamilies>0):
            famID = 1
        else:
            famID = 0
        #endif
        
        sumb.initializewarping(famID)

        self.nGlobalSurfNodes = sumb.mddata.mdnsurfnodescompact
        return

    def SetGlobalSurfaceCoordinates(self,xyz,reinitialize=True):
        """Set the surface coordinates for the mesh                                                                                
        Keyword arguments:
                                                                                
        xyz -- the new xyz coordinates, dimension (3,ncoords)
                                                                                        
        """
        sumb.iteration.groundlevel = 1
        #ncoords = xyz.shape[1]
        #print 'ncoords',ncoords,xyz.shape
        #sys.exit(0)
        #sumb.updatefacesglobal(ncoords,xyz)
        xyz = self.metricConversion*xyz


        sumb.updatefacesglobal(xyz,reinitialize)

    
        self._update_geom_info = True

        return

    def GetGlobalSurfaceCoordinates(self):
        """
        Get the surface coordinates for the mesh    
        """
        if self.sol_type == 'Steady' or self.sol_type == 'steady':
            return sumb.mddata.mdglobalsurfxx/self.metricConversion
        elif self.sol_type == 'Time Spectral' or self.sol_type == 'time spectral':
            return sumb.mddata.mdglobalsurfxx[:,:,0]/self.metricConversion
        else:
            print 'invalid solutions type for surface coords...Exiting'
            sys.exit(0)
        #endif
        return

    def warpMesh(self):
        '''
        run the internal meshwarping scheme
        '''

        sumb.warpmesh()

        return

    def warpMeshSolid(self,*args,**kwargs):
        '''
        run the Solid Mesh Warping scheme
       '''

        # Define several functions required for initialization
        def inBinarySearch(list,element):
            ileft = bisect.bisect_left(list,element)
            if ileft >= len(list):
                return False
            if list[ileft] == element:
                return True
            else:
                return False
        # end def
        
        def convertFaces(face_list): # Convert faces from pyPSG to SUmb
            new_face_list = zeros(6,face_list.dtype)
            new_face_list[0] = face_list[2] # ilo
            new_face_list[1] = face_list[3] # ihi
            new_face_list[2] = face_list[4] # jlo
            new_face_list[3] = face_list[5] # jhi
            new_face_list[4] = face_list[0] # klo
            new_face_list[5] = face_list[1] # khi
            return new_face_list

        def convertEdges(edge_list): # Convert Edges from pyPSG to SUmb
            new_edge_list = zeros(12,edge_list.dtype)
            new_edge_list[0] = edge_list[0]
            new_edge_list[1] = edge_list[1]
            new_edge_list[2] = edge_list[4]
            new_edge_list[3] = edge_list[5]
            new_edge_list[4] = edge_list[2]
            new_edge_list[5] = edge_list[3]
            new_edge_list[6] = edge_list[6]
            new_edge_list[7] = edge_list[7]
            new_edge_list[8] = edge_list[8]
            new_edge_list[9] = edge_list[9]
            new_edge_list[10] = edge_list[10]
            new_edge_list[11] = edge_list[11]
            return new_edge_list

        # We have initialization stuff to do ONLY on the root prco
        if not self.solid_warp_initialized:
            if self.myid == 0:
                mpiPrint('\nInitializating Solid Mesh Warping...')
                import pyspline # Direct access to the compiled spline library

                if 'file' in kwargs: # we're using a separate CGNS file
                    file_name = kwargs['file']
                elif 'n' in kwargs:
                    file_name = self.GetFilename() #we're using mesh file with constant 'n'
                    n = kwargs['n']
                elif 'topo' in kwargs:
                    file_name = self.GetFilename()
                else:
                    mpiPrint('Error: Keyword arguments \'n=integer\' \'file=<cgns_file>\' or \'topo=topo_file\' must be passed to warpMeshSolid',comm=comm)

                if 'sym' in kwargs:
                    if kwargs['sym'] == 'xy' or kwargs['sym'] == 'yx':
                        sym = [0,0,-1]
                    elif kwargs['sym'] == 'yz' or kwargs['sym'] == 'zy':
                        sym = [-1,0,0]
                    elif kwargs['sym'] == 'xz' or kwargs['sym'] == 'zx':
                        sym = [0,-1,0]
                    else:
                        mpiPrint('  ## Error: sym must be one of xy, yz or xz')
                        sys.exit(0)
                    # end if
                else:
                    mpiPrint('  ** Warning: sym is not specified.',comm=self.comm)
                    sym = [0,0,0]
                # end if

                cg,nzones = pyspline.open_cgns(file_name)
                sizes = []
                BCs = []
                corners = numpy.zeros((nzones,8,3))
                for i in xrange(nzones):
                    zoneshape = pyspline.read_cgns_zone_shape(cg,i+1)
                    X,faceBCs = pyspline.read_cgns_zone(cg,i+1,zoneshape[0],zoneshape[1],zoneshape[2])
                    sizes.append(zoneshape)
                    BCs.append(faceBCs)
                    corners[i][0] = X[0,0,0]
                    corners[i][1] = X[-1,0,0]
                    corners[i][2] = X[0,-1,0]
                    corners[i][3] = X[-1,-1,0]
                    corners[i][4] = X[0,0,-1]
                    corners[i][5] = X[-1,0,-1]
                    corners[i][6] = X[0,-1,-1]
                    corners[i][7] = X[-1,-1,-1]
                # end for
                pyspline.close_cgns(cg) # Done with cgns file
                sizes = array(sizes) # make sizes into a numpy array
                if 'n' in kwargs:
                    FE_topo = BlockTopology(corners) # Compute topology
                    sizes[:,:] = n
                    FE_topo.calcGlobalNumbering(sizes)
                elif 'file' in kwargs:
                    FE_topo = BlockTopology(corners) # Compute topology
                    FE_topo.calcGlobalNumbering(sizes)
                elif 'topo' in kwargs:
                    FE_topo = BlockTopology(file=kwargs['topo'])
                    FE_topo.calcGlobalNumbering()
                 # end if
                nuu,nus,l_index_flat,l_ptr,l_sizes = self._reOrderIndices(FE_topo,BCs,sym) # Re-order numbering to account for constrained dof

                self.FE_topo = FE_topo
                solidWarpData = [nuu,nus,array(l_index_flat),array(l_ptr),array(l_sizes)]
            else:
                self.FE_topo = None
                solidWarpData = None
            # end if (myID == 0 test)
                
            # Bcast data from root node to all
            self.FE_topo  = self.comm.bcast(self.FE_topo,root=0)
            solidWarpData = self.comm.bcast(solidWarpData,root=0)

            # Also produce the global to local mapping for the
            # mdsurfacenodescompact back to each of the blocks 
            
            data = self.comm.gather(sumb.mddatalocal.mdsurfglobalindlocal,root=0)
            if self.myid == 0:
                md_g_index = [[] for i in xrange(sumb.mddata.mdnsurfnodescompact)]
                # Now loop over data from each proc and each node on the proc
                for iproc in xrange(len(data)):
                    for ii in xrange(data[iproc].shape[1]):
                        i  = data[iproc][0,ii] # i index on this block
                        j  = data[iproc][1,ii] # j index on this block
                        k  = data[iproc][2,ii] # k index on this block
                        nn = data[iproc][3,ii] # domain number of process
                        id = data[iproc][4,ii] # global node id
                                        
                        md_g_index[id].append([i,j,k,iproc,nn])
                    # end for
                # end for
                     
                # Now flatten the g_index along with a pointer
                md_g_ptr = zeros((len(md_g_index)+1),'intc')
                md_g_ptr[0] = 0# Zerobased Here
                for i in xrange(len(md_g_index)):
                    md_g_ptr[i+1] = md_g_ptr[i] + len(md_g_index[i])*5
                # end for

                # Dont' ask..it works since we only have 1 level deep...
                md_g_index = array([item for sublist in md_g_index for item in sublist]).flatten()
            else:
                md_g_index = None
                md_g_ptr = None

            # Now bcast g_index and g_ptr back to everyone
            md_g_index = self.comm.bcast(md_g_index,root=0)
            md_g_ptr   = self.comm.bcast(md_g_ptr,root=0)
            
            # Set the required data in the module
            sumb.solidwarpmodule.nuu      = solidWarpData[0]
            sumb.solidwarpmodule.nus      = solidWarpData[1]
            sumb.solidwarpmodule.l_index  = solidWarpData[2]
            sumb.solidwarpmodule.lptr     = solidWarpData[3]
            sumb.solidwarpmodule.l_sizes  = solidWarpData[4]
            sumb.solidwarpmodule.nblock   = len(solidWarpData[2])
            sumb.solidwarpmodule.md_g_index = md_g_index
            sumb.solidwarpmodule.md_g_ptr   = md_g_ptr

            # Now run the fortran initialization
            sumb.initializewarpmeshsolid_parallel()

            # Initialization Complete
            mpiPrint('  -> Solid Mesh Warping Initialized.',comm=self.comm)
            self.solid_warp_initialized = True
        # end if Initialization

        # Now run the actual warping command
        sumb.warpmeshsolid_parallel()

        return

    def warpMeshSolidDeriv(self,*args,**kwargs):
        sumb.calculatesolidwarpderiv()

    def _reOrderIndices(self,FE_topo,faceBCs,sym):
        '''This funcion takes the order from self.FE_topo and reorders
        them according to the Boundary Condition types in each volume
        class. The sole purpose of this is to facilitate the
        structural mesh warping algorithim for matrix assembly.'''
        # We want the global indicies ordered according to:
        #[ freedof ]
        #[ constrained dof ]
        sym = array(sym)
        pt_dof = numpy.zeros((FE_topo.nGlobal,3),'intc')

        for ii in xrange(FE_topo.nGlobal):
            for jj in xrange(len(FE_topo.g_index[ii])):

                ivol = FE_topo.g_index[ii][jj][0]
                i    = FE_topo.g_index[ii][jj][1]
                j    = FE_topo.g_index[ii][jj][2]
                k    = FE_topo.g_index[ii][jj][3]

                N = FE_topo.l_index[ivol].shape[0]
                M = FE_topo.l_index[ivol].shape[1]
                L = FE_topo.l_index[ivol].shape[2]
              
                type,number,index1,index2 = \
                    indexPosition3D(i,j,k,N,M,L)
                
                checkFaces = []
                if type == 0:
                    pt_dof[ii] = [0,0,0]
                else: # Face
                    if type == 1: # Face
                        checkFaces.append(number)
                    elif type == 2: # Edge
                        if number in [0,1,2,3]:
                            checkFaces.append(0)
                        if number in [4,5,6,7]:
                            checkFaces.append(1)
                        if number in [2,6,8,10]:
                            checkFaces.append(2)
                        if number in [3,7,9,11]:
                            checkFaces.append(3)
                        if number in [0,4,8,9]:
                            checkFaces.append(4)
                        if number in [1,5,10,11]:
                            checkFaces.append(5)
                    elif type == 3: # Corner
                        if number == 0:
                            checkFaces.extend([0,2,4])
                        elif number == 1:
                            checkFaces.extend([0,3,4])
                        elif number == 2:
                            checkFaces.extend([0,2,5])
                        elif number == 3:
                            checkFaces.extend([0,3,5])
                        elif number == 4:
                            checkFaces.extend([1,2,4])
                        elif number == 5:
                            checkFaces.extend([1,3,4])
                        elif number == 6:
                            checkFaces.extend([1,2,5])
                        elif number == 7:
                            checkFaces.extend([1,3,5])
                    # end if
                    
                    # We now now all faces a point that belong to a
                    # pt, check each for boun dary conditions
                    for iii in xrange(len(checkFaces)):
                        iface = checkFaces[iii]
                        if faceBCs[ivol][iface] in [1,2]: # BC_wall OR Farfield/inflow
                            pt_dof[ii] = [1,1,1]
                        # end if
                        if faceBCs[ivol][iface] == 3:
                            # Only set it as a symmetry plane if nothing is already set
                            index_set = where(sym==-1)[0][0]
                            pt_dof[ii][index_set] = 1
                        # end if
                    # end for
                # end if
            # end for
        # end for
      
        nus = int(sum(sum(pt_dof)))
        nuu = FE_topo.nGlobal*3-nus
        mpiPrint('  -> Total DOF  : %d'%(FE_topo.nGlobal*3),comm=self.comm)
        mpiPrint('  -> Unknown DOF: %d'%(nuu),comm=self.comm)
        mpiPrint('  -> Known DOF  : %d'%(nus),comm=self.comm)

        # We will forgo the g_index reorganization...it is not
        # strictly necessary We want l_index[ivol] to be of size
        # (nu,nv,nw,3) with each entry pointing to the dof in the
        # global matrix

        free_dof_count = 0
        constr_dof_count = 0
        l_index = []
        for ivol in xrange(len(FE_topo.l_index)):
            l_index.append(zeros((FE_topo.l_index[ivol].shape[0],
                                  FE_topo.l_index[ivol].shape[1],
                                  FE_topo.l_index[ivol].shape[2],3),'intc'))
        # end for

        for ii in xrange(FE_topo.nGlobal):
            for iii in xrange(3):
                if pt_dof[ii][iii] == 0:
                    for jj in xrange(len(FE_topo.g_index[ii])):
                        ivol = FE_topo.g_index[ii][jj][0]
                        i    = FE_topo.g_index[ii][jj][1]
                        j    = FE_topo.g_index[ii][jj][2]
                        k    = FE_topo.g_index[ii][jj][3]
                        l_index[ivol][i,j,k,iii] = free_dof_count
                    # end for
                    free_dof_count += 1
                # end if
                if pt_dof[ii][iii] == 1:
                    for jj in xrange(len(FE_topo.g_index[ii])):
                        ivol = FE_topo.g_index[ii][jj][0]
                        i    = FE_topo.g_index[ii][jj][1]
                        j    = FE_topo.g_index[ii][jj][2]
                        k    = FE_topo.g_index[ii][jj][3]
                        l_index[ivol][i,j,k,iii] = nuu + constr_dof_count
                    # end for
                    constr_dof_count += 1
                # end if
            # end for (iii loop)
        # end for (ii loop)

        # Lastly, we need to flatten the l_index for fortran use

        l_index_flat = []
        l_ptr = [0] # -> Zero Based Here
        l_sizes = zeros((len(l_index),3),'intc')
        for i in xrange(len(l_index)):
            l_index_flat.extend(l_index[i].flatten())
            l_ptr.append(l_ptr[-1] + l_index[i].size)
            l_sizes[i] = [l_index[i].shape[0],l_index[i].shape[1],l_index[i].shape[2]]
        # end for

        return nuu,nus,l_index_flat,l_ptr,l_sizes

    def GetMeshQuality(self,file_name):
        # Return a list of the quality of all elements
        
        # Step 1: Determine the number of elements on this processor
        nElem = 0
        for i in xrange(sumb.block.ndom):
            size = self.getBlockDimensions(i+1)
            nElem += (size[0]-1)*(size[1]-1)*(size[2]-1)
        # end for
        quality = sumb.getquality(nElem)

        # Dump it out to a file
        f = open(file_name+'_%d'%(self.myid),'w')
        quality.tofile(f,sep="\n",format="%20.16g")
        
        return     


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

        # The very first thing --> Set the MPI Communicators
        if 'comm' in kwargs:
            sumb.communication.sumb_comm_world = kwargs['comm'].py2f()
            sumb.communication.sumb_comm_self  = mpi.COMM_SELF.py2f()
            self.sumb_comm_world = kwargs['comm']
        else:
            sumb.communication.sumb_comm_world = mpi.COMM_WORLD.py2f()
            sumb.communication.sumb_comm_self  = mpi.COMM_SELF.py2f()
            self.sumb_comm_world = mpi.COMM_WORLD
        # end if
        
        if 'init_petsc' in kwargs:
            if kwargs['init_petsc']:
                sumb.initializepetsc()
        # end if

        # Setup the mesh object with sumb_comm_world
        self.Mesh = SUmbMesh(self.sumb_comm_world)

        # Determine the rank sumb_comm_world size
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
            sys.exit(1)
        # end try

        # Set the SUmb module value of standalonemode to false and
        # the value of deforming_grid to the input value.
        sumb.iteration.standalonemode = False
        sumb.iteration.deforming_grid = deforming_mesh = False

        # Write the intro message
        sumb.writeintromessage()

        # Set the frompython flag to true
        sumb.killsignals.frompython=True

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
        if kwargs['solver_options']:
            self.probName = kwargs['solver_options']['probName']
            self.OutputDir = kwargs['solver_options']['OutputDir']
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
        
        sumb.inputio.paramfile[0:len(startfile)] = startfile

        # Read the parameter file
        sumb.readparamfile()

        #This is just to flip the -1 to 1 possibly a memory issue?
        sumb.inputio.storeconvinneriter=abs(sumb.inputio.storeconvinneriter)

        if(self.myid ==0):print ' -> Partitioning and Reading Grid'
        sumb.partitionandreadgrid()

        if(self.myid==0):print ' -> Preprocessing'
        sumb.preprocessing()

        if(self.myid==0):print ' -> Initializing flow'
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

        # Create Surface Node list
        if(self.myid==0): print ' -> Creating Surface Node List'
        sumb.mdcreatensurfnodes()

        # Reduce the total number of blocks
        self.Mesh.nmeshblocks = self.sumb_comm_world.allreduce(
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

        #if (self.myid==0):print 'Alpha',aero_problem._flows.alpha*(pi/180.0),aero_problem._flows.alpha,velDir,liftDir,dragDir
        #update the flow vars
        sumb.updateflow()
        return
    def resetFlow(self):
        '''
        Reset the flow for the complex derivative calculation
        '''

        sumb.setuniformflow()

        return
    
    def generateInputFile(self,aero_problem,sol_type,grid_file,startfile,file_type='cgns',eqn_type='Euler',*args,**kwargs):
        ''' Code to generate an SUmb Input File on the fly'''
        if (self.myid==0): print ' -> Generating Input File'
        
        #Convert alpha and beta to a freestream vector
        [velDir,liftDir,dragDir]= sumb.adjustinflowangleadj((aero_problem._flows.alpha*(pi/180.0)),(aero_problem._flows.beta*(pi/180.0)),aero_problem._flows.liftIndex)
        
        #autofile = open("autogen.input",'w')
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
            autofile.write("                           Restart: %s\n"%(kwargs['solver_options']['sol_restart']))
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
        try: kwargs['solver_options']['FamilyRot']
        except KeyError:
            Rotating = False
        else:
            Rotating = True
        #endif
        #print 'Rotating',Rotating,kwargs['solver_options']['FamilyRot']
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
        autofile.write(  "#                         Reynolds: 100000\n")
        autofile.write(  "       Reynolds length (in meter): 1.0\n")
        #autofile.write(  "   Free stream velocity direction: 1.0 0.05 0.0\n")
        #autofile.write(  "                   Lift direction: -0.05 1.0 0.0\n")
        autofile.write(  "   Free stream velocity direction: %12.12e %12.12e %12.12e\n"%(velDir[0],velDir[1],velDir[2]))
        autofile.write(  "                   Lift direction: %12.12e %12.12e %12.12e\n"%(liftDir[0],liftDir[1],liftDir[2]))
        autofile.write(  "     # Default is normal to free stream without y-component\n")
        autofile.write(  "   Free stream temperature (in K): %12.12e\n"%(kwargs['solver_options']['Reference Temp.']))
        autofile.write(  " Free stream eddy viscosity ratio: 0.01\n")
        autofile.write(  "  Free stream turbulent intensity: 0.001\n")
        autofile.write(  "\n")
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Reference State\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "            Reference pressure (in Pa): %12.12e\n"%(kwargs['solver_options']['Reference Pressure']))
        autofile.write(  "         Reference density (in kg/m^3): 1.25\n")
        autofile.write(  "          Reference temperature (in K): %12.12e\n"%(kwargs['solver_options']['Reference Temp.']))
        autofile.write(  " Conversion factor grid units to meter: %6.4f\n"%(kwargs['solver_options']['MetricConversion']))
        self.Mesh.metricConversion = kwargs['solver_options']['MetricConversion']
        autofile.write( "\n")
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Geometrical Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "           Reference surface: %12.12e\n"%(aero_problem._refs.sref*kwargs['solver_options']['MetricConversion']**2))

        #autofile.write(  "           Reference surface: %2.1f\n"%(1.0))

        autofile.write(  "            Reference length: %12.12e\n"%(aero_problem._refs.cref*kwargs['solver_options']['MetricConversion']))
        autofile.write(  "    Moment reference point x:  %12.12e\n"%(aero_problem._refs.xref*kwargs['solver_options']['MetricConversion']))
        autofile.write(  "    Moment reference point y:  %12.12e\n"%(aero_problem._refs.yref*kwargs['solver_options']['MetricConversion']))
        autofile.write(  "    Moment reference point z:  %12.12e\n"%(aero_problem._refs.zref*kwargs['solver_options']['MetricConversion']))
        autofile.write( "\n")
        
        #! Write the keywords and default values for the discretization
        #! parameters.

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Fine Grid Discretization Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "       Discretization scheme: %s\n"%(kwargs['solver_options']['Discretization']))
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
        autofile.write(  "   Exponent dissipation scaling: %4.2f\n"%(kwargs['solver_options']['Dissipation Scaling Exponent']))
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
            autofile.write("             Number time intervals spectral: %d\n"%(kwargs['solver_options']['Time Intervals']))
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

        autofile.write(  "                 Number of multigrid cycles: 200\n")
        autofile.write(  "   Number of single grid startup iterations: 0\n")
        autofile.write(  "                                 Save every: 10\n")
        autofile.write(  "                         Save surface every: 10\n")
        autofile.write(  "                                 CFL number: %2.1f\n"%(kwargs['solver_options']['CFL']))
        autofile.write( "\n")
        if eqn_type=='RANS':
            autofile.write(  "                       Turbulent relaxation: Explixit\n")
            autofile.write(  "                            # Possibilities: Explicit\n")
            autofile.write(  "                            #              : Implicit\n")
            autofile.write(  "                     Alpha turbulent DD-ADI: 0.8\n")
            autofile.write(  "                      Beta turbulent DD-ADI: -1  # Same as alpha\n")
        #endif
        autofile.write(  "           Relative L2 norm for convergence: %3.2e\n"%(kwargs['solver_options']['L2Convergence']))
        autofile.write( "\n")

        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "     Multigrid Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "      Number of multigrid cycles coarse grid:  -1  # -1 Means same as on fine grid\n")
        autofile.write(  "                      CFL number coarse grid: -1  # -1 Means same as on fine grid\n")

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
        autofile.write(  "                 Multigrid cycle strategy: %s\n"%(kwargs['solver_options']['MGCycle']))
        autofile.write( "\n")

#
# ADjoint Parameters
#
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "      ADjoint Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "                       solve ADjoint :  %s\n"%(kwargs['solver_options']['solveADjoint']))
        autofile.write(  "                  # Other possibility: no\n")
        autofile.write(  "	                 set monitor   :  %s\n"%(kwargs['solver_options']['set Monitor']))
        autofile.write(  "                  # Other possibility: no\n")
        autofile.write(  "       Use Approximate Preconditioner: %s\n"%(kwargs['solver_options']['Approx PC']))
        autofile.write(  "                  # Other possibility: yes/no\n")
        autofile.write(  " 	            Adjoint solver type: %s\n"%(kwargs['solver_options']['Adjoint solver type']))
        autofile.write(  "                # Other possibilities: BiCGStab\n")
        autofile.write(  "                #                      CG\n")
        autofile.write(  "                #                      GMRES\n")
        autofile.write(  "                #                      FGMRES\n")
        autofile.write(  "        adjoint relative tolerance   : %3.2e\n"%(kwargs['solver_options']['adjoint relative tolerance']))
        autofile.write(  "        adjoint absolute tolerance   : %3.2e\n"%(kwargs['solver_options']['adjoint absolute tolerance']))
        autofile.write(  "        adjoint divergence tolerance : 1e5\n")
        autofile.write(  "        adjoint max iterations       : %d\n"%kwargs['solver_options']['adjoint max iterations'])
        autofile.write(  "        adjoint restart iteration    : %d\n"%kwargs['solver_options']['adjoint restart iteration'])
        autofile.write(  "        adjoint monitor step         : %d\n"%kwargs['solver_options']['adjoint monitor step'])
        autofile.write(  "        dissipation lumping parameter: %d\n"%kwargs['solver_options']['dissipation lumping parameter'])
        
        autofile.write(  "                  Preconditioner Side: %s\n"%(kwargs['solver_options']['Preconditioner Side']))
        autofile.write(  "                # Other possibilities: Left\n")
        autofile.write(  "                #                      Right\n")
        autofile.write(  "	             Matrix Ordering   : %s\n"%(kwargs['solver_options']['Matrix Ordering']))
        autofile.write(  "	               #                 ReverseCuthillMckee\n")
        autofile.write(  "	               #                 Natural\n")
        autofile.write(  "	               #                 NestedDissection\n")
        autofile.write(  "                     #                 OnewayDissection\n")
        autofile.write(  "                     #                 QuotientMinimumDegree\n")
        autofile.write(  "          Global Preconditioner Type : %s\n"%(kwargs['solver_options']['Global Preconditioner Type']))
        autofile.write(  "          #                            Jacobi\n")
        autofile.write(  "          #                            Block Jacobi\n")
        autofile.write(  "          #                            Additive Schwartz\n")
        autofile.write(  "                         ASM Overlap : 5\n")
        autofile.write(  "            Local Preconditioner Type: %s\n"%(kwargs['solver_options']['Local Preconditioner Type']))
        autofile.write(  "            #                          ILU\n")
        autofile.write(  "            #                          ICC\n")
        autofile.write(  "            #                          LU\n")
        autofile.write(  "            #                          Cholesky\n")
        autofile.write(  "                    ILU Fill Levels  : %d\n"%kwargs['solver_options']['ILU Fill Levels'])
        autofile.write(  "            Jacobi Scale Factor Type : RowMax\n")
        autofile.write(  "            #                          RowMax\n")
        autofile.write(  "            #                          RowSum\n")
        autofile.write(  "            #                          RowAbs\n")

        #Write the options for the Time Spectral stability derivatives to the input file
        
        autofile.write(  "-------------------------------------------------------------------------------\n")
        autofile.write(  "TS Stability Derivative Parameters\n")
        autofile.write(  "-------------------------------------------------------------------------------\n")
	autofile.write(  "compute TS stability derivatives : %s\n"%(kwargs['solver_options']['TS Stability']))
	autofile.write(  "#Other Possibilities:               no\n")
        if kwargs['solver_options']['TS Stability']== 'yes':
            autofile.write(  "TS Alpha mode: %s\n"%(kwargs['solver_options']['Alpha Mode']))
            autofile.write(  "TS Beta mode: %s\n"%(kwargs['solver_options']['Beta Mode']))
            autofile.write(  "TS p mode: %s\n"%(kwargs['solver_options']['p Mode']))
            autofile.write(  "TS q mode: %s\n"%(kwargs['solver_options']['q Mode']))
            autofile.write(  "TS r mode: %s\n"%(kwargs['solver_options']['r Mode']))
            autofile.write(  "TS Mach number mode: %s\n"%(kwargs['solver_options']['Mach Mode']))
            autofile.write(  "TS Altitude mode: %s\n"%(kwargs['solver_options']['Altitude Mode']))
        #endif

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

        if sol_type=='Time Spectral':
##        ! Write the keywords and example values for the grid motion.
            
            autofile.write("-------------------------------------------------------------------------------\n")
            autofile.write("     Grid motion Parameters\n")
            autofile.write( "-------------------------------------------------------------------------------\n")

            autofile.write( "     Rotation point body (x,y,z): %12.12e %12.12e %12.12e\n"%(kwargs['solver_options']['rotCenter'][0],kwargs['solver_options']['rotCenter'][1],kwargs['solver_options']['rotCenter'][2]))
            autofile.write("\n" )

            # inputs for TS stability derivatives
            autofile.write( "Degree polynomial Alpha: 0\n")
            autofile.write( "Degree polynomial Beta: 0\n")
            autofile.write( "Degree polynomial Mach: 0\n")
            
            autofile.write( "Polynomial coefficients Alpha: 0.0\n")
            autofile.write( "Polynomial coefficients Beta: 0.0\n")
            autofile.write( "Polynomial coefficients Mach: 0.0\n")

            if kwargs['solver_options']['Alpha Mode'] =='yes':
                autofile.write( "Degree fourier Alpha: 1\n")
                autofile.write( "Degree fourier Beta: 0\n")
                autofile.write( "Degree fourier Mach: 0#1\n")
            elif kwargs['solver_options']['Beta Mode'] =='yes':
                autofile.write( "Degree fourier Alpha: 0\n")
                autofile.write( "Degree fourier Beta: 1\n")
                autofile.write( "Degree fourier Mach: 0\n")
            elif kwargs['solver_options']['Mach Mode'] =='yes':
                autofile.write( "Degree fourier Alpha: 0\n")
                autofile.write( "Degree fourier Beta: 0\n")
                autofile.write( "Degree fourier Mach: 1\n")
            #endif

            if kwargs['solver_options']['Alpha Mode'] =='yes':
                autofile.write( "Omega fourier Alpha: %f\n"%(kwargs['solver_options']['Omega fourier']))
                autofile.write( "Omega fourier Beta: 0.0\n")
                autofile.write( "Omega fourier Mach: 0.0\n")
            elif kwargs['solver_options']['Beta Mode'] =='yes':
                autofile.write( "Omega fourier Alpha: 0.0\n")
                autofile.write( "Omega fourier Beta: %f\n"%(kwargs['solver_options']['Omega fourier']))
                autofile.write( "Omega fourier Mach: 0.0\n")
            elif kwargs['solver_options']['Mach Mode'] =='yes':
                autofile.write( "Omega fourier Alpha: 0.0\n")
                autofile.write( "Omega fourier Beta: 0.0\n")
                autofile.write( "Omega fourier Mach: %f\n"%(kwargs['solver_options']['Omega fourier']))
            #endif

            if kwargs['solver_options']['Alpha Mode'] =='yes':
                autofile.write( "Fourier cosine coefficients Alpha: 0.0 0.0\n")
                autofile.write( "Fourier sine coefficients Alpha: %12.12e\n"%(kwargs['solver_options']['Fourier sine coefficient']))
            elif kwargs['solver_options']['Beta Mode'] =='yes':
                autofile.write( "Fourier cosine coefficients Beta: 0.0 0.0\n")
                autofile.write( "Fourier sine coefficients Beta: %12.12e\n"%(kwargs['solver_options']['Fourier sine coefficient']))
            elif kwargs['solver_options']['Mach Mode'] =='yes':
                autofile.write( "Fourier cosine coefficients Mach: 0.0 0.0\n")
                autofile.write( "Fourier sine coefficients Mach: %12.12e\n"%(kwargs['solver_options']['Fourier sine coefficient']))
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
        autofile.write( "                Monitoring variables: resrho_cl_cd_cmx_cmy_cmz\n")
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
            autofile.write( "Rotating family %s : %6.6f %6.6f %6.6f    0.e+0 0.e+0 0.e+0    \n"%(kwargs['solver_options']['FamilyRot'],kwargs['solver_options']['rotCenter'][0],kwargs['solver_options']['rotCenter'][1],kwargs['solver_options']['rotCenter'][2]))
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
        
#         if(self.myid ==0):
#             print 'Iterations...',sumb.monitor.niterold,sumb.monitor.nitercur,\
#                 sumb.iteration.itertot
#         # end if
        sumb.monitor.niterold = self.sumb_comm_world.bcast(sumb.monitor.niterold,root=0)

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
        sumb.killsignals.routinefailed=False

        if sol_type=='Steady' or sol_type=='steady' or sol_type=='Time Spectral' or sol_type=='time spectral':
            #if(self.myid ==0):print 'steady',sumb.inputio.storeconvinneriter,True,False,sumb.inputio.storeconvinneriter==True,sumb.inputio.storeconvinneriter==False
            #if(self.myid ==0):sumb.inputio.storeconvinneriter=True
            #if(self.myid ==0):print 'steady',sumb.inputio.storeconvinneriter,True,False,sumb.inputio.storeconvinneriter==True,sumb.inputio.storeconvinneriter==False
            #set the number of cycles for this call
            sumb.inputiteration.ncycles = ncycles
            
            if (sumb.monitor.niterold == 0 and \
                sumb.monitor.nitercur == 0 and \
                sumb.iteration.itertot == 0):
                # No iterations have been done
                if (sumb.inputio.storeconvinneriter):
                    nn = sumb.inputiteration.nsgstartup+sumb.inputiteration.ncycles
                    if(self.myid==0):
                        # sumb.monitor.convarray = None
                       sumb.deallocconvarrays()
                       sumb.allocconvarrays(nn)
                       #print 'convarray',sumb.monitor.convarray
                    #endif
                #endif

            elif(sumb.monitor.nitercur == 0 and\
                 sumb.iteration.itertot == 0):
                
                # Reallocate convergence history array and
                # time array with new size, storing old values from restart
                if (self.myid == 0):
                    # number of time steps from restart
                    ntimestepsrestart = sumb.monitor.ntimestepsrestart
                                        
                    if (sumb.inputio.storeconvinneriter):
                        # number of iterations from restart
                        niterold = sumb.monitor.niterold#[0]
                        if storeHistory:
                            # store restart convergence history and deallocate array
                            temp = copy.deepcopy(sumb.monitor.convarray[:niterold+1,:])
                            sumb.deallocconvarrays()
                            # allocate convergence history array with new extended size
                            sumb.allocconvarrays(temp.shape[0]
                                                 +sumb.inputiteration.ncycles-1)
                            # recover values from restart and deallocate temporary array
                            sumb.monitor.convarray[:temp.shape[0],:temp.shape[1]] = temp
                            temp = None
                        else:
                            temp=copy.deepcopy(sumb.monitor.convarray[0,:,:])
                            sumb.deallocconvarrays()
                            # allocate convergence history array with new extended size
                            sumb.allocconvarrays(sumb.inputiteration.ncycles+1+niterold)
                            sumb.monitor.convarray[0,:,:] = temp
                            temp = None
                        #endif
                    #endif
                #endif
            else:
                # More Time Steps / Iterations in the same session
            
                # Reallocate convergence history array and
                # time array with new size, storing old values from previous runs
                if (self.myid == 0):
                    if (sumb.inputio.storeconvinneriter):
                        if storeHistory:
                            # store previous convergence history and deallocate array
                            temp = copy.deepcopy(sumb.monitor.convarray)

                            sumb.deallocconvarrays()
                            # allocate convergence history array with new extended size
                            nn = sumb.inputiteration.nsgstartup+sumb.inputiteration.ncycles
                            print 'alloc'
                            sumb.allocconvarrays(temp.shape[0]+nn-1)
                        
                            # recover values from previous runs and deallocate temporary array
                            sumb.monitor.convarray[:temp.shape[0],:] = copy.deepcopy(temp)
                            
                            temp = None
                        else:
                            temp=copy.deepcopy(sumb.monitor.convarray[0,:,:])
                            sumb.deallocconvarrays()
                            sumb.allocconvarrays(sumb.inputiteration.ncycles+1)
                            sumb.monitor.convarray[0,:,:] = temp
                            temp = None
                        #endif
                    #endif
                #endif
                # re-initialize iteration variables
                sumb.inputiteration.mgstartlevel = 1
                if storeHistory:
                    sumb.monitor.niterold  = sumb.monitor.nitercur
                else:
                    sumb.monitor.niterold  = 1#0#sumb.monitor.nitercur
                #endif
                sumb.monitor.nitercur  = 0#1
                sumb.iteration.itertot = 0#1
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
        elif sol_type=='unsteady':
            print 'unsteady not implemented yet...'
            sys.exit(0)
        else:
            print 'unknown sol_type:',sol_type,'exiting...'
            sys.exit(0)
        #endif

        #if self.myid ==0:print 'setupfinished'
        self.GetMesh()._UpdateGeometryInfo()
        #print 'killsignal',sumb.killsignals.routinefailed,True,self.myid
        self.routineFailed = self.sumb_comm_world.allreduce(sumb.killsignals.routinefailed,mpi.MIN)
        #print 'routinefailed',self.routineFailed,self.myid
        if (abs(self.routineFailed)==True):
            if self.myid ==0:print 'Error raise in updateGeometry'
            raise ValueError
            #print 'error raised'
            #return
        #endif
    
        #if self.myid ==0: print 'calling solver'
        sumb.solver()
        #Check to see whether we have a valid solution
        print 'killsignal',sumb.killsignals.routinefailed,True,self.myid
        #self.routineFailed = self.sumb_comm_world.allreduce(sumb.killsignals.routinefailed,mpi.MIN)
        #in this case routineFailed will be triggered on all processors
        self.routineFailed = self.sumb_comm_world.allreduce(abs(sumb.killsignals.routinefailed),mpi.SUM)/mpi.COMM_WORLD.size
        #therefore sum up and devide by nProc
        print 'routinefailed',self.routineFailed,self.myid
        if (abs(self.routineFailed)==True):
            #print 'raising error'
            if self.myid ==0: print 'Error Raised in solver'
            raise ValueError
            #print 'error raised'
            #return
        #endif
        if self.myid ==0:print 'solver called'

        return
        sys.exit(0)
        
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
                      if(self.myid==0):print 'setting ncycles'
                      sumb.inputiteration.ncycles = ncycles[1]
                # Reallocate convergence history array and
                # time array with new size
                if (self.myid == 0):
                    print 'sumb.monitor.timearray'#,sumb.monitor.timearray 
                    sumb.monitor.timearray = None
                    print 'sumb.monitor.timearray',sumb.monitor.timearray 
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
                    ntimestepsrestart = sumb.monitor.ntimestepsrestart#[0]
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
        if self.myid ==0:print 'setupfinished'
        #self.Mesh.WriteMeshFile('newmesh.cgns')
        if self.myid ==0: print 'calling solver'
        self.GetMesh()._UpdateGeometryInfo()

        if sumb.killsignals.routinefailed==True:
            raise ValueError
            return    
        #endif
        
        #print 'mesh info updated'
        #self.Mesh.WriteMeshFile('newmesh0.cgns')
        #print 'new file written'
        sumb.solver()
        if self.myid ==0:print 'solver called'
        
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
	        
        #Set the mesh level and timespectral instance for this
        #computation
        self.level = 1
        self.sps = 1
        
        sumb.iteration.currentlevel=1
        sumb.iteration.groundlevel=1

        #Check to see if initialization has already been performed
        if(self.adjointInitialized):
            return
        
        #Run the preprocessing routine. Sets the node numbering and
        #allocates memory.
        if(self.myid==0): print 'preprocessing adjoint'
        sumb.preprocessingadjoint(self.level)
        
        #Initialize the design variable and function storage
        if(self.myid==0):print 'Before design init'
        sumb.designinit()

        #initalize PETSc
        if(self.myid==0):print 'before petsc'
        sumb.initializepetsc()

        #create the neccesary PETSc objects
        if(self.myid==0):print 'before createpetsecars'
        sumb.createpetscvars()

        #mark the ADjoint as initialized
        self.adjointInitialized = True
        if(self.myid==0):
            print 'ADjoint Initialized Succesfully...'
        #endif
        if(self.myid==0):print 'before nspatial'
        self.nSpatial = sumb.adjointvars.ndesignspatial

        return

    def setupADjointMatrix(self):
        '''
        Setup the ADjoint dRdw matrix and create the PETSc
        Solution KSP object
        '''
        #create the neccesary PETSc objects
        if(self.myid==0):print 'before createpetscvars'
        #sumb.createpetscvars()
        #sumb.createpetscmat()
        #sumb.setupadjointmatrix(self.level)
        sumb.setupadjointmatrixtranspose(self.level)

        sumb.setuppetscksp(self.level)

        return

    def releaseAdjointMemeory(self):
        '''
        release the KSP memory...
        '''
        #sumb.destroypetscksp()
        #sumb.destroypetscvars()
        #sumb.destroypetscmat()
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
                             }

        if self.myid==0:
            print 'Analyzing ADjoint for costfuntion: ',objective
            print 'SUmb index:',SUmbCostfunctions[objective]
    
        #print SUmbCostfunctions[objective]

        sumb.setupadjointrhs(self.level,SUmbCostfunctions[objective])

        return

    def solveADjointPETSc(self):
        '''
        Solve the ADjoint system using PETSc
        '''
        #sumb.solveadjointpetsc()
        sumb.solveadjointtransposepetsc()

        return

    def verifyPartials(self):
        '''
        Run solverADjoint to verify the partial derivatives in the ADjoint
        '''
        print 'in interface verify partials'
        #sumb.verifydrdxsfile()
        sumb.verifydcfdx(1)
        #sumb.solveradjoint()

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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
        
        sumb.setupgradientrhsvolume(self.level,SUmbCostfunctions[possibleObjectives[objective]])
        #sumb.setupgradientrhsspatial(self.level,SUmbCostfunctions[possibleObjectives[objective]],self.sps)
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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

        sumb.setupgradientrhsflow(self.level,SUmbCostfunctions[possibleObjectives[objective]])
        #sumb.setupgradientrhsextra(self.level,SUmbCostfunctions[possibleObjectives[objective]],self.sps)
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

    def setupVolumeSurfaceDerivatives(self):
        """Set up the derivative of the volume mesh  w.r.t. the
        CFD Surface...(meshwarpingderivatives)"""
        
        sumb.setupvolumesurfacederivativesdv()
        
	if (self.myid == 0):
            print "Meshwarping derivatives set up successfully."
        #endif
        return
    
    def computeTotalSurfaceDerivative(self,objective):
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
        
        sumb.computeadjointgradientsurfacedv(SUmbCostfunctions[possibleObjectives[objective]])
        #endfor

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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
        
        grad = numpy.zeros((sumb.adjointvars.ndesignspatial),float)
        grad[:] = sumb.adjointvars.functiongradspatial[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
        #for item in objective:
        
        #for i in xrange(sumb.adjointvars.ndesignspatial):
            #grad[i] = sumb.adjointvars.functiongradspatial[SUmbCostfunctions[possibleObjectives[objective]]-1,i]
           
        #endfor
        #endfor
        
        return grad
    def getTotalSurfaceDerivatives(self,objective):

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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
        
        grad = numpy.zeros((3*sumb.mddata.mdnsurfnodescompact),float)
        grad[:] = sumb.adjointvars.functiongradsurfacedv[SUmbCostfunctions[possibleObjectives[objective]]-1,:]
   
        
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
        if(self.myid==0):print 'in sumbInterface'
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
        print 'sumb value:',sumb.adjointvars.costfuncliftcoef-1

        SUmbsolutions = {'cl':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncliftcoef-1],\
                         'cd':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncdragcoef-1],\
                         'cFx':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncforcexcoef-1],\
                         'cFy':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncforceycoef-1],\
                         'cFz':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncforcezcoef-1],\
                         'cMx':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncmomxcoef-1],\
                         'cMy':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncmomycoef-1],\
                         'cMz':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncmomzcoef-1],\
                         'cMzAlpha':sumb.adjointvars.functionvalue[sumb.adjointvars.costfunccmzalpha-1],\
                         'cM0':sumb.adjointvars.functionvalue[sumb.adjointvars.costfunccm0-1],\
                         'clAlpha':sumb.adjointvars.functionvalue[sumb.adjointvars.costfuncclalpha-1],\
                         'cl0':sumb.adjointvars.functionvalue[sumb.adjointvars.costfunccl0-1],\
                         'cdAlpha':sumb.adjointvars.functionvalue[sumb.adjointvars.costfunccdalpha-1],\
                         'cd0':sumb.adjointvars.functionvalue[sumb.adjointvars.costfunccd0-1]
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                             'cMzAlpha':sumb.adjointvars.costfunccmzalpha,\
                             'cM0':sumb.adjointvars.costfunccm0,\
                             'clAlpha':sumb.adjointvars.costfuncclalpha,\
                             'cl0':sumb.adjointvars.costfunccl0,\
                             'cdAlpha':sumb.adjointvars.costfunccdalpha,\
                             'cd0':sumb.adjointvars.costfunccd0
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
                
            sumb.setupadjointtotalstruct(structAdjoint,SUmbCostfunctions[possibleObjectives[objective]])
        except:
            print 'not an aerodynamic cost function'
        #end

        return

    def computeStabilityParameters(self):
        '''
        run the stability derivative driver to compute the stability parameters
        from the time spectral solution
        '''

        sumb.stabilityderivativedriver()

        return
    









# Only python version to compute ifaceptb and iedgept...not complete...new version 
# in fortran renders this unnecessary
#   # Now produce ifaceptb,iedgeptb using the topology
#                 #information and the face boundary condition data

#                 # Step 1: Loop over all faces on all blocks, for each
#                 # face, add the ID of each of the four nodes to the
#                 # explicitly perturbed corner list

#                 explicitly_perturbed_corners = []
#                 for ivol in xrange(FE_topo.nVol):
#                     for iface in xrange(6):
#                         if BCs[ivol][iface] == 1: # BCwallViscous
#                             corners = nodesFromFace(iface)
#                             for icorner in corners:
#                                 explicitly_perturbed_corners.append(FE_topo.node_link[ivol][icorner])
#                             # end for
#                         # end if
#                     # end for
#                 # end for

#                 # Step 2: Unique-ify the explicitly_perturbed_corners list
#                 # and sort them in the process

#                 explicitly_perturbed_corners = unique(explicitly_perturbed_corners)
#                 explicitly_perturbed_corners.sort() # It should alreadly be sorted, but just in case

#                 # Step 3: Loop over all edges on all volumes. If they have
#                 # ONE explictly perturbed corner they are IMPLICTLY
#                 # perturbed, if they have TWO explictly perturbed corners,
#                 # they are EXPLICTLY perturbed. Note we call
#                 # inBinarySearch here. This is simpilar to python's 'in'
#                 # function. However since we have a sorted list, we can do
#                 # a binary search instead of a linear search which is MUCH
#                 # faster. The return is True or False as to whether
#                 # 'element' is in 'list'

#                 iedgeptb = zeros((FE_topo.nVol,12),'intc')

#                 for ivol in xrange(FE_topo.nVol):
#                     for iedge in xrange(12):
#                         edge_number = FE_topo.edge_link[ivol][iedge]
#                         corners = [FE_topo.edges[edge_number].n1,FE_topo.edges[edge_number].n2]
#                         explicit = array([False,False])
#                         for icorner in xrange(2):
#                             explicit[icorner] = inBinarySearch(\
#                                 explicitly_perturbed_corners,corners[icorner])
#                         # end for

#                         if explicit.all():
#                             iedgeptb[ivol][iedge] = 2
#                         elif explicit.any():
#                             iedgeptb[ivol][iedge] = 1
#                         else:
#                             iedgeptb[ivol][iedge] = 0
#                         # end if
#                     # end for
#                 # end for

#                 # Step 4: Loop over all faces on all volumes. An EXPLICTLY
#                 # perturbed face will have ALL corners explictly
#                 # perturbed. An IMPLICLTY perturbed face will have at
#                 # least EXPLICTLY perturbed corner. A face will have value
#                 # of 0 if there are NO EXPLICTLY perturbed corners.

#                 ifaceptb = zeros((FE_topo.nVol,6),'intc')

#                 for ivol in xrange(FE_topo.nVol):
#                     for iface in xrange(6):
#                         corners = nodesFromFace(iface)
#                         explicit = array([False,False,False,False])
#                         for icorner in xrange(4):
#                             node = FE_topo.node_link[ivol][corners[icorner]]
#                             explicit[icorner] = inBinarySearch(\
#                                 explicitly_perturbed_corners,node)
#                         # end for
#                         if explicit.all():
#                             ifaceptb[ivol][iface] = 2
#                         elif explicit.any():
#                             ifaceptb[ivol][iface] = 1
#                         else:
#                             ifaceptb[ivol][iface] = 0
#                         # end if
#                     # end for
#                 # end for
                            
#                 # Step 5: Finally permute the faces and edges from pyPSG
#                 # format to to SUmb format
#                 for ivol in xrange(FE_topo.nVol):
#                     ifaceptb[ivol] = convertFaces(ifaceptb[ivol])
#                     iedgeptb[ivol] = convertEdges(iedgeptb[ivol])
#                 # end for
