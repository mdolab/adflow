#! /usr/bin/env python
#
# File: sumbsolution.py
# Authors: Patrick LeGresley, Andre Marta, Steve Repsher
# Last Modified: 9/11/2005

__version__ = "1.03"

# =============================================================================
# Standard Python modules

import sys
import os
import copy
import string
import types

# =============================================================================
# Extension modules

import Numeric
from flowsolution import FlowSolution
from mesh import Mesh

# Try to import the mpi module
try:
    import mpi
    _parallel = True
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


def CheckIfParallel():
    """Return true if we are running with MPI."""
    return _parallel

# =============================================================================

class SUmbMesh(Mesh):
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

    def DeallocateSurfaceCoordinates(self):
        """Deallocate memory used for the surface coordinates."""
        sumb.mddeletesurfcoorlist()

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
                                             
    def DeallocateSurfaceIndices(self):
        """Deallocate memory used for the surface indices."""
        sumb.mddeletesurfindlist()

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
        blocknums = Numeric.zeros((0))
        ranges = Numeric.zeros((3,2,0))
        xyz = Numeric.zeros((3,0),'d')
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


# =============================================================================

class SUmbSolution(FlowSolution):
    """Represents a SUmb flow solution."""

    def __init__(self,startfile,communicator=None,deforming_mesh=False):
        """Initialize the object.

        Keyword arguments:

        startfile      -- the name of the SUmb input parameter file
        communicator   -- an MPI communicator for specifying which processors
                          to run on (optional)
        deforming_mesh -- set to True or False whether or not this is a 
                          deforming mesh case (optional, default is False)
                          

        """
        # Make sure the parameter file exists
        if not os.path.isfile(startfile):
            print 'Error: Could not find file %s' % startfile
            return None

        # Setup a mesh object
        self.Mesh = SUmbMesh()

        # Create a Python communicator to mirror the Fortran
        # SUmb_COMM_WORLD and set the Fortran SUmb communicator to the
        # Python one.
        if (communicator is not None):
            self.sumb_comm_world = communicator
        else:
            self.sumb_comm_world = mpi.COMM_WORLD.comm_create(
        				     mpi.COMM_WORLD[:])
        sumb.communication.sumb_comm_world = int(self.sumb_comm_world)
        self.Mesh.sumb_comm_world = self.sumb_comm_world

        # Store the name of the input file
        self.startfile = startfile
        sumb.inputio.paramfile[0:len(startfile)] = startfile

        # Determine the rank and number of processors inside the group
        # defined by sumb_comm_world.
        self.myid = sumb.communication.myid = self.sumb_comm_world.rank
        self.nproc = sumb.communication.nproc = self.sumb_comm_world.size

        # Allocate the memory for SENDREQUESTS and RECVREQUESTS.
        try:
            sumb.communication.sendrequests = Numeric.zeros(
        			  (self.sumb_comm_world.size))
            sumb.communication.recvrequests = Numeric.zeros(
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

        # Read the parameter file
        sumb.readparamfile()

        # Partition the blocks and read the grid
        sumb.partitionandreadgrid()
        
        # Perform the preprocessing task
        sumb.preprocessing()
        
        # Initialize the flow variables
        sumb.initflow()

        # Create dictionary of variables we are monitoring
        nmon = sumb.monitor.nmon[0]
        self.monnames = {}
        for i in range(nmon):
            self.monnames[string.strip(
        		   sumb.monitor.monnames[i].tostring())] = i
        
        # Create dictionary of the family names
        sumb.mdgetfamilynames()
        nfamilies = sumb.mddata.mdnfamilies[0]
        self.Mesh.families = {}
        for i in range(nfamilies):
            self.Mesh.families[string.strip(sumb.mddata.mdfamilynames[i]
        		       .tostring())] = i + 1
        sumb.mdcreatensurfnodes()

        # Determine the total number of blocks in the mesh and store it
        self.Mesh.nmeshblocks = self.sumb_comm_world.allreduce(
        			     sumb.block.ndom[0],mpi.SUM)

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
                        sumb.monitor.convarray = None
                        sumb.allocconvarrays(sumb.inputunsteady.ntimestepsfine*sumb.inputiteration.ncycles)

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
            if (ncycles):
                # Set new value of unsteady physical time steps to run
                if (ncycles[0] != sumb.inputunsteady.ntimestepsfine):
                    sumb.inputunsteady.ntimestepsfine = ncycles[0]
                # Set new value of MG cycles ro run
                if (len(ncycles) > 1):
                    if (ncycles[1] != sumb.inputiteration.ncycles):
                        sumb.inputiteration.ncycles = ncycles[1]

            # Reallocate convergence history array and
            # time array with new size, storing old values from previous runs
            if (self.myid == 0):
                # store previous time history and deallocate arrays
                temp_t = copy.deepcopy(sumb.monitor.timearray)
                sumb.monitor.timearray = None
                temp_td = copy.deepcopy(sumb.monitor.timedataarray)
                sumb.monitor.timedataarray = None
                # allocate time history arrays with new extended size
                sumb.alloctimearrays(temp_td.shape[0]+sumb.inputunsteady.ntimestepsfine)
                # recover values from previous runs and deallocate temporary arrays
                sumb.monitor.timearray[:temp_td.shape[0]] = temp_t
                sumb.monitor.timedataarray[:temp_td.shape[0],:temp_td.shape[1]] = temp_td
                temp_t = None
                temp_td = None

                if (sumb.inputio.storeconvinneriter):
                    # store previous convergence history and deallocate array
                    temp = copy.deepcopy(sumb.monitor.convarray)
                    sumb.monitor.convarray = None
                    # allocate convergence history array with new extended size
                    sumb.allocconvarrays(temp.shape[0]
                                         +sumb.inputunsteady.ntimestepsfine
                                         *sumb.inputiteration.ncycles-1)
                    # recover values from previous runs and deallocate temporary array
                    sumb.monitor.convarray[:temp.shape[0],:temp.shape[1]] = temp
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

        self.GetMesh()._UpdateGeometryInfo()
#        sumb.solver()
################################################################
# The code below reproduces solver.F90 in python
# See file solver.F90
################################################################

        if (sumb.inputdiscretization.wallbctreatment == sumb.inputdiscretization.normalmomentum):
           sumb.iteration.exchangepressureearly = sumb.eulerwallspresent()
        else:
           sumb.iteration.exchangepressureearly = False

        sumb.killsignals.localsignal == sumb.killsignals.nosignal
        try:
            sumb.connectsignals()
        except AttributeError:
            pass

        sumb.iteration.t0solver = sumb.su_wtime()

        sumb.monitor.timeunsteady = 0.0

        if (not sumb.iteration.pv3initialized):
            sumb.iteration.groundlevel = 1
            try:
                sumb.initalizepv3()
                sumb.iteration.pv3initialized = True
            except AttributeError:
                pass
        
        for sumb.iteration.groundlevel in range(sumb.inputiteration.mgstartlevel,0,-1):
            if (sumb.inputphysics.equationmode == sumb.inputphysics.steady or \
                sumb.inputphysics.equationmode == sumb.inputphysics.timespectral):
                sumb.solversteady()
            elif (sumb.inputphysics.equationmode == sumb.inputphysics.unsteady):
#                sumb.solverunsteady()
################################################################
# The code below reproduces solverUnsteady.F90 in python
# See file solverUnsteady.F90
################################################################
#               sumb.monitor.writevolume  = False
#               sumb.monitor.writesurface = False
#               sumb.monitor.writegrid    = False

                sumb.monitor.timestepunsteady = 0

                ntimesteps = sumb.inputunsteady.ntimestepscoarse
                if (sumb.iteration.groundlevel == 1):
                    ntimesteps =  sumb.inputunsteady.ntimestepsfine

                for iter in range(1,ntimesteps+1):
                    sumb.inittimesteppart1()
                    ##############################  
                    #sumb.additional_routine_here#
                    ##############################
                    sumb.inittimesteppart2()

                    sumb.solvestate()
#                   sumb.iteration.noldsolavail = sumb.iteration.noldsolavail + 1

#ifdef USE_PV3
#                   CALL PV_UPDATE(REAL(TIMESTEPUNSTEADY,REALPV3TYPE))
#endif

                    sumb.checkwriteunsteadyinloop()

                    if (sumb.killsignals.globalsignal == sumb.killsignals.signalwritequit):
                        break

                sumb.checkwriteunsteadyendloop()
################################################################
# end of solverUnsteady.F90
################################################################

            if (sumb.iteration.groundlevel > 1):

                sumb.iteration.currentlevel = sumb.iteration.groundlevel - 1

                if (sumb.communication.myid == 0):
                    print "#"
                    print "# Going down to grid level %i" % sumb.iteration.currentlevel
                    print "#"

                sumb.transfertofinegrid(False)

                if (sumb.inputphysics.equationmode == sumb.inputphysics.unsteady and \
                    sumb.iteration.changing_grid):
                    sumb.updatecoorfinemesh(sumb.monitor.timeunsteady, 1)

                sumb.iteration.noldsolavail = 1
################################################################
# end of solver.F90
################################################################

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
            sumb.inputio.surfacesolfile[0:len(filename)] = filename
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
