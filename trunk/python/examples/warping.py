#! /usr/bin/env python
#
# Demonstrates the linking of pywarp and pytflo.  This script will run serially
# with python or in parallel (at least the flow solver part) with pyMPI.
#
# By: Patrick LeGresley
# Last updated: 6/14/2004

# Import standard modules
import sys
import os
import copy
import math
import Numeric

# Try to import the mpi module
sys.path.append(os.path.abspath('..'))
try:
    import mpi
except ImportError:
    try:
        import dummy_mpi as mpi
    except ImportError:
        print "Error: Failed to import mpi or dummy_mpi."
my_id = mpi.rank

# Import class representing SUmb
from sumbsolution import SUmbSolution

# Get an instance of the SUmb class
myflow = SUmbSolution('startfile')

# Only run the mesh warping on one processor
if (my_id == 0):
    # Import class for mesh warping
    sys.path.append(os.path.abspath('../../../mdo-python/pywarp/python'))
    from multiblockmesh import MultiblockMesh
 
    # Get an instance of the mesh
    mymesh = MultiblockMesh(myflow.GetMesh().GetFilename())

# Run iterations.  With no arguments the number will be
# based on the input in the initialization file
myflow.RunIterations()

# Copy the original coordinates from block 97
if (my_id == 0):
    blnum = 97
    [blocknum,ijkrange,xyz_orig] = mymesh.GetCoordinates(blnum)
    xyz_orig = copy.deepcopy(xyz_orig)
    [imax,jmax,kmax] = mymesh.GetBlockDimensions(blnum)
    npts = imax*kmax
    xyz = Numeric.zeros((3,npts),'d')
    blocknums = blnum*Numeric.ones((npts))
    ijk = Numeric.zeros((3,npts))

# Keep looping and warping the mesh
for ntimes in range(1000):

    # Make some new coordinates to pass into warp
    if (my_id == 0):
        count = 0
        count_local = 0

        # Amplitude of the perturbation
        deltay = 0.75*(math.sin(float(ntimes)*0.1))

        # Take the long list of coordinates for the block
        # and make a shorter list with only coordinates for
        # the face on the solid surface.
        for k in range(kmax):
            for j in range(jmax):
                for i in range(imax):
                    if (j==0):
                        ijk[0,count_local] = i + 1
                        ijk[1,count_local] = j + 1
                        ijk[2,count_local] = k + 1
                        xyz[0,count_local] = xyz_orig[0,count]
                        xyz[1,count_local] = xyz_orig[1,count] + deltay
                        xyz[2,count_local] = xyz_orig[2,count]
                        count_local = count_local + 1
                    count = count + 1

        # The actual warping
        mymesh.Warp(xyz,blocknums,ijk)

    # Loop over the blocks and insert the coordinates into SUmb
    for n in range(1,myflow.GetMesh().GetNumberBlocks()+1):
        if (my_id == 0):
            # Get the new coordinates from warp
            [blocknum,ijkrange,xyz_new] = mymesh.GetCoordinates(n)
            # Give them to SUmb
            myflow.GetMesh().SetCoordinates(blocknum,ijkrange,xyz_new)
        else:
            # Make a dummy call on processors where warp is not running
            myflow.GetMesh().DummySetCoordinates()

    # Run more iterations
    myflow.RunIterations(1)
