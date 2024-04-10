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
import types
import numpy
import sys
from mpi4py import MPI
from baseclasses import AeroSolver, AeroProblem, getPy3SafeString
from baseclasses.utils import Error, CaseInsensitiveDict
from . import MExt
import hashlib
from collections import OrderedDict


class ADFLOWWarning(object):
    """
    Format a warning message
    """

    def __init__(self, message):
        msg = "\n+" + "-" * 78 + "+" + "\n" + "| pyADFLOW Warning: "
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78:  # Finish line and start new one
                msg += " " * (78 - i) + "|\n| " + word + " "
                i = 1 + len(word) + 1
            else:
                msg += word + " "
                i += len(word) + 1
        msg += " " * (78 - i) + "|\n" + "+" + "-" * 78 + "+" + "\n"
        print(msg)


# =============================================================================
# ADFLOW Class
# =============================================================================
class ADFLOW(AeroSolver):
    """
    Create the ADflow object.

    Parameters
    ----------
    comm : MPI intra comm
        The communicator on which to create ADflow. If not given, defaults
        to MPI.COMM_WORLD.
    options : dict
        The list of options to use with ADflow. This keyword argument
        is NOT OPTIONAL. It must always be provided. It must contain, at least
        the 'gridFile' entry for the filename of the grid to load
    debug : bool
        Set this flag to true when debugging with a symbolic
        debugger. The MExt module deletes the copied .so file when not
        required which causes issues debugging.
    dtype : str
        String type for float: 'd' or 'D'. Not needed to be used by user.
    """

    def __init__(self, comm=None, options=None, debug=False, dtype="d"):
        startInitTime = time.time()

        # Define the memory allocation flags ASAP (needed for __del__)
        # Matrix Setup Flag
        self.adjointSetup = False

        # Load the compiled module using MExt, allowing multiple
        # imports
        try:
            self.adflow
        except AttributeError:
            curDir = os.path.basename(os.path.dirname(os.path.realpath(__file__)))
            self.adflow = MExt.MExt("libadflow", curDir, debug=debug)._module

        libLoadTime = time.time()

        # Information for base class:
        name = "ADFLOW"
        category = "Three Dimensional CFD"
        informs = {}

        # Load all the option/objective/DV information:
        defaultOptions = self._getDefaultOptions()
        self.optionMap, self.moduleMap = self._getOptionMap()
        self.pythonOptions, deprecatedOptions, self.specialOptions = self._getSpecialOptionLists()
        immutableOptions = self._getImmutableOptions()
        self.rootChangedOptions = {}

        self.possibleAeroDVs, self.possibleBCDvs, self.basicCostFunctions = self._getObjectivesAndDVs()

        # Now add the group for each of the "basic" cost functions:
        self.adflowCostFunctions = OrderedDict()
        for key in self.basicCostFunctions:
            self.adflowCostFunctions[key] = [None, key]

        # Separate list of the supplied supplied functions
        self.adflowUserCostFunctions = OrderedDict()

        # This is the real solver so dtype is 'd'
        self.dtype = dtype

        # Next set the MPI Communicators and associated info
        if comm is None:
            comm = MPI.COMM_WORLD

        self.adflow.communication.adflow_comm_world = comm.py2f()
        self.adflow.communication.adflow_comm_self = MPI.COMM_SELF.py2f()
        self.adflow.communication.sendrequests = numpy.zeros(comm.size)
        self.adflow.communication.recvrequests = numpy.zeros(comm.size)
        self.myid = self.adflow.communication.myid = comm.rank
        self.adflow.communication.nproc = comm.size

        if options is None:
            raise Error(
                "The 'options' keyword argument must be passed "
                "adflow. The options dictionary must contain (at least) "
                "the gridFile entry for the grid."
            )

        # Set all internal adflow default options before we set anything from python
        self.adflow.inputparamroutines.setdefaultvalues()

        defSetupTime = time.time()

        # Initialize the inherited AeroSolver
        super().__init__(
            name,
            category,
            defaultOptions=defaultOptions,
            options=options,
            immutableOptions=immutableOptions,
            deprecatedOptions=deprecatedOptions,
            comm=comm,
            informs=informs,
        )

        baseClassTime = time.time()
        # Update turbresscale depending on the turbulence model specified
        self._updateTurbResScale()

        # Initialize PETSc in case the user has not already
        self.adflow.adjointutils.initializepetsc()

        # Set the stand-alone adflow flag to false...this changes how
        # terminate calls are handled.
        self.adflow.iteration.standalonemode = False

        # Set the frompython flag to true... this also changes how
        # terminate calls are handled
        self.adflow.killsignals.frompython = True

        # Dictionary of design variables and their index
        self.aeroDVs = []

        # Default counters
        self.updateTime = 0.0
        self.nSlice = 0
        self.nLiftDist = 0

        # Set default values
        self.adflow.inputio.autoparameterupdate = False
        self._updateGeomInfo = True

        # By Default we don't have an external mesh object or a
        # geometric manipulation object
        self.mesh = None
        self.DVGeo = None
        self.zipperCreated = False
        self.userSurfaceIntegrationsFinalized = False
        self.hasIntegrationSurfaces = False

        # Write the intro message
        if self.getOption("printIntro"):
            self.adflow.utils.writeintromessage()

        # Remind the user of all the adflow options:
        if self.getOption("printAllOptions"):
            self.printOptions()

        # Do the remainder of the operations that would have been done
        # had we read in a param file
        self.adflow.iteration.deforming_grid = True

        # In order to properly initialize we need to have mach number
        # and a few other things set. Just create a dummy aeroproblem,
        # use it, and then it will be deleted.

        dummyAP = AeroProblem(
            name="dummy",
            mach=0.5,
            altitude=10000.0,
            areaRef=1.0,
            chordRef=1.0,
            alpha=0.0,
            degreePol=0,
            coefPol=[0, 0],
            degreeFourier=1,
            omegaFourier=6.28,
            sinCoefFourier=[0, 0],
            cosCoefFourier=[0, 0],
        )

        self._setAeroProblemData(dummyAP, firstCall=True)
        # Now set it back to None so the user is none the wiser
        self.curAP = None

        # Finally complete loading
        self.adflow.inputparamroutines.dummyreadparamfile()
        if self.getOption("partitionOnly"):
            self.adflow.partitioning.partitionandreadgrid(True)
            return

        introTime = time.time()

        self.adflow.partitioning.partitionandreadgrid(False)

        partitioningTime = time.time()

        self.adflow.preprocessingapi.preprocessing()

        preprocessingTime = time.time()

        # Do explict hole cutting.
        ncells = self.adflow.adjointvars.ncellslocal[0]
        ntime = self.adflow.inputtimespectral.ntimeintervalsspectral
        n = ncells * ntime

        # Normally we would be able to just use the following...but
        # sometimes f2py is grumpy and we get
        # ValueError: data type must provide an itemsize
        # So we do it he hard way...sigh
        nFam = self.adflow.surfacefamilies.getnfam()
        famList = []
        for i in range(nFam):
            famList.append(getPy3SafeString(self.adflow.surfacefamilies.getfam(i + 1).strip()))

        # Add the initial families that already exist in the CGNS
        # file.
        for i in range(len(famList)):
            self.families[famList[i]] = [i + 1]

        # Add the special "all surfaces" family.
        self.allFamilies = "allSurfaces"
        self.addFamilyGroup(self.allFamilies, famList)

        # We also need to know about which surfaces are walls. Pull
        # out that information from fortran and set the special
        # familyGroup.
        wallList, nWalls = self.adflow.surfaceutils.getwalllist(len(famList))
        wallList = wallList[0:nWalls]
        wallListFam = []
        for i in range(len(wallList)):
            wallListFam.append(famList[wallList[i] - 1])
        self.allWallsGroup = "allWalls"
        self.addFamilyGroup(self.allWallsGroup, wallListFam)

        # Set the design families if given, otherwise default to all
        # walls
        self.designFamilyGroup = self.getOption("designSurfaceFamily")
        if self.designFamilyGroup is None:
            self.designFamilyGroup = self.allWallsGroup

        # Set the mesh families if given, otherwise default to all
        # walls
        self.meshFamilyGroup = self.getOption("meshSurfaceFamily")
        if self.meshFamilyGroup is None:
            self.meshFamilyGroup = self.allWallsGroup

        familySetupTime = time.time()

        # Get the CGNS zone name info for the user
        self.nZones = self.adflow.utils.getncgnszones()
        self.CGNSZoneNameIDs = {}
        for i in range(self.nZones):
            # index needs to go in as fortran numbering, so add 1
            name = getPy3SafeString(self.adflow.utils.getcgnszonename(i + 1).strip())
            self.CGNSZoneNameIDs[name] = i + 1

        # Call the user supplied callback if necessary
        cutCallBack = self.getOption("cutCallBack")
        flag = numpy.zeros(n)
        if cutCallBack is not None:
            xCen = self.adflow.utils.getcellcenters(1, n).T
            cellIDs = self.adflow.utils.getcellcgnsblockids(1, n)
            cutCallBack(xCen, self.CGNSZoneNameIDs, cellIDs, flag)

        cutCallBackTime = time.time()

        # exclude the cells inside closed surfaces if we are provided with them
        explicitSurfaceCallback = self.getOption("explicitSurfaceCallback")
        if explicitSurfaceCallback is not None:
            # the user wants to exclude cells that lie within a list of surfaces.

            # first, call the callback function with cgns zone name IDs.
            # this need to return us a dictionary with the surface mesh information,
            # as well as which blocks in the cgns mesh to include in the search
            surfDict = explicitSurfaceCallback(self.CGNSZoneNameIDs)

            # loop over the surfaces
            for surf in surfDict:
                if self.comm.rank == 0:
                    print(f"Explicitly blanking surface: {surf}")

                # this is the plot3d surface that defines the closed volume
                surfFile = surfDict[surf]["surfFile"]
                # the indices of cgns blocks that we want to consider when blanking inside the surface
                blockIDs = surfDict[surf]["blockIDs"]
                # the fortran lookup expects this list in increasing order
                blockIDs.sort()

                # check if there is a kMin provided
                if "kMin" in surfDict[surf]:
                    kMin = surfDict[surf]["kMin"]
                else:
                    kMin = -1

                # optional coordinate transformation to do general manipulation of the coordinates
                if "coordXfer" in surfDict[surf]:
                    coordXfer = surfDict[surf]["coordXfer"]
                else:
                    coordXfer = None

                # read the plot3d surface
                pts, conn = self._readPlot3DSurfFile(surfFile, convertToTris=False, coordXfer=coordXfer)

                # get a new flag array
                surfFlag = numpy.zeros(n, "intc")

                # call the fortran routine to determine if the cells are inside or outside.
                # this code is very similar to the actuator zone creation.
                self.adflow.oversetapi.flagcellsinsurface(pts.T, conn.T, surfFlag, blockIDs, kMin)

                # update the flag array with the new info
                flag = numpy.any([flag, surfFlag], axis=0)

        explicitSurfaceCutTime = time.time()

        # Need to reset the oversetPriority option since the CGNSGrid
        # structure wasn't available before;
        self.setOption("oversetPriority", self.getOption("oversetPriority"))

        # Set the closed surface families if given. Otherwise
        # default to all walls.
        famList = self.getOption("closedSurfaceFamilies")
        if famList is None:
            self.closedFamilyGroup = self.allWallsGroup
        else:
            self.addFamilyGroup("closedSurface", famList)
            self.closedFamilyGroup = "closedSurface"

        famList = self._getFamilyList(self.closedFamilyGroup)
        self.adflow.preprocessingapi.preprocessingoverset(flag, famList)

        oversetPreTime = time.time()

        self.adflow.tecplotio.initializeliftdistributiondata()
        self.adflow.initializeflow.updatebcdataalllevels()
        self.adflow.initializeflow.initflow()

        initFlowTime = time.time()

        self.coords0 = self.getSurfaceCoordinates(self.allFamilies, includeZipper=False)

        finalInitTime = time.time()

        if self.getOption("printTiming") and self.comm.rank == 0:
            print("+--------------------------------------------------+")
            print("|")
            print("| Initialization Times:")
            print("|")
            print("| %-30s: %10.3f sec" % ("Library Load Time", libLoadTime - startInitTime))
            print("| %-30s: %10.3f sec" % ("Set Defaults Time", defSetupTime - libLoadTime))
            print("| %-30s: %10.3f sec" % ("Base Class Init Time", baseClassTime - defSetupTime))
            print("| %-30s: %10.3f sec" % ("Introductory Time", introTime - baseClassTime))
            print("| %-30s: %10.3f sec" % ("Partitioning Time", partitioningTime - introTime))
            print("| %-30s: %10.3f sec" % ("Preprocessing Time", preprocessingTime - partitioningTime))
            print("| %-30s: %10.3f sec" % ("Family Setup Time", familySetupTime - preprocessingTime))
            print("| %-30s: %10.3f sec" % ("Cut Callback Time", cutCallBackTime - familySetupTime))
            print("| %-30s: %10.3f sec" % ("Explicit Surface Cut Time", explicitSurfaceCutTime - cutCallBackTime))
            print("| %-30s: %10.3f sec" % ("Overset Preprocessing Time", oversetPreTime - explicitSurfaceCutTime))
            print("| %-30s: %10.3f sec" % ("Initialize Flow Time", initFlowTime - oversetPreTime))
            print("|")
            print("| %-30s: %10.3f sec" % ("Total Init Time", finalInitTime - startInitTime))
            print("+--------------------------------------------------+")

    def __del__(self):
        """
        Clean up allocated memory if necessary
        """
        # Release PETSc memory
        self.releaseAdjointMemory()
        if hasattr(self, "adflow"):
            self.adflow.nksolver.destroynksolver()
            self.adflow.anksolver.destroyanksolver()

            # Release Fortran memory
            # Check these fortran routines, they are not complete.
            # However, they do work and the left over memory
            # is not too large.
            self.adflow.utils.releasememorypart1()
            self.adflow.utils.releasememorypart2()

    def setMesh(self, mesh):
        """
        Set the mesh object to the aero_solver to do geometric deformations

        Parameters
        ----------
        mesh : MBMesh or USMesh object
            The mesh object for doing the warping
        """

        # Store a reference to the mesh
        self.mesh = mesh

        # Setup External Warping with volume indices
        meshInd = self.getSolverMeshIndices()
        self.mesh.setExternalMeshIndices(meshInd)

        # Set the surface the user has supplied:

        conn, faceSizes, cgnsBlockIDs = self.getSurfaceConnectivity(
            self.meshFamilyGroup, includeZipper=False, includeCGNS=True
        )

        pts = self.getSurfaceCoordinates(self.meshFamilyGroup, includeZipper=False)

        self.mesh.setSurfaceDefinition(pts, conn, faceSizes, cgnsBlockIDs)

    def getSolverMeshIndices(self):
        """
        Get the list of indices to pass to the mesh object for the
        volume mesh mapping
        """
        ndof_1_instance = self.adflow.adjointvars.nnodeslocal[0] * 3
        meshInd = self.adflow.warping.getcgnsmeshindices(ndof_1_instance)

        return meshInd

    def setDisplacements(self, aeroProblem, dispFile):
        """
        This function allows the user to perform aerodynamic
        analysis/optimization while using a fixed set of displacements
        computed from a previous structural analysis. Essentially this
        allows the jig shape to designed, but performing analysis on
        the flying shape.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
           The AP object that the displacements should be applied to.
        dispFile : str
           The file contaning the displacements. This file should have
           been obtained from TACS

        Notes
        -----
        The fixed set of displacements do not affect the sensitivities,
        since they are fixed and not affected by any DVs.

        Also, in the case where the current surface mesh was not used
        to generate the displacements file, a nearest neighbor search
        is used to apply the displacements.
        """
        self.setAeroProblem(aeroProblem)

        # All processors read the file
        f = open(dispFile, "r")
        n = int(f.readline())
        X = []
        D = []
        for _i in range(n):
            aux = f.readline().split()
            X.append([float(aux[0]), float(aux[1]), float(aux[2])])
            D.append([float(aux[3]), float(aux[4]), float(aux[5])])
        f.close()

        localX = self.getSurfaceCoordinates()
        localD = numpy.zeros_like(localX)
        # Now we need to search each localX in X to find the corresponding D
        try:
            from scipy.spatial import KDTree
        except ImportError:
            raise Error("scipy must be available to use setDisplacements")
        tree = KDTree(numpy.array(X))
        d, index = tree.query(localX)
        for j in range(len(localX)):
            localD[j] = D[index[j]]

        # Fianlly we have the local displacements for this processor we save
        self.curAP.adflowData.disp = localD

    def addLiftDistribution(self, nSegments, direction, groupName=None):
        """Add a lift distribution to the surface output.

        Parameters
        ----------
        nSegments : int
            Number of slices to use for the distribution. Typically
            150-250 is sufficient

        direction : str {'x', 'y', 'z'}
            The normal of the slice direction. If you have z-axis 'out the wing'
            you would use 'z'

        groupName: str
            The family group to use for the lift distribution. Default
            of None corresponds to all wall surfaces.
        """
        # Create the zipper mesh if not done so
        self._createZipperMesh()

        # Determine the families we want to use
        if groupName is None:
            groupName = self.allWallsGroup

        groupTag = "%s: " % groupName
        famList = self._getFamilyList(groupName)
        direction = direction.lower()
        if direction not in ["x", "y", "z"]:
            if direction not in ["x", "y", "z"]:
                raise Error("direction must be one of 'x', 'y', or 'z'")

        if direction == "x":
            dirVec = [1.0, 0.0, 0.0]
            dirInd = 1
        elif direction == "y":
            dirVec = [0.0, 1.0, 0.0]
            dirInd = 2
        else:
            dirVec = [0.0, 0.0, 1.0]
            dirInd = 3

        distName = "LiftDist_%2.2d %s: %s normal" % (self.nLiftDist + 1, groupTag, direction)
        self.nLiftDist += 1

        self.adflow.tecplotio.addliftdistribution(nSegments, dirVec, dirInd, distName, famList)

    def addSlices(self, direction, positions, sliceType="relative", groupName=None):
        """
        Add parametric slice positions. Slices are taken of the wing
        at time addParaSlices() is called and the parametric positions
        of the intersections on the surface mesh are stored. On
        subsequent output, the position of the slice moves as the mesh
        moves/deforms. This effectively tracks the same location on
        the wing.

        Parameters
        ----------
        direction : str {'x', 'y', 'z'}
            The normal of the slice direction. If you have z-axis 'out the wing'
            you would use 'z'
        positions : float or array
            The list of slice positions *along the axis given by direction*.
        sliceType : str {'relative', 'absolute'}
            Relative slices are 'sliced' at the beginning and then parametrically
            move as the geometry deforms. As a result, the slice through the
            geometry may not remain planar. An absolute slice is re-sliced for
            every output so is always exactly planar and always at the initial
            position the user indicated.
        groupName : str
             The family to use for the slices. Default is None corresponding to all
             wall groups.
        """

        # call the preprocessing routine that avoids some code duplication across 3 types of slices
        sliceType, famList = self._preprocessSliceInput(sliceType, groupName)

        direction = direction.lower()
        if direction not in ["x", "y", "z"]:
            raise Error("'direction' must be one of 'x', 'y', or 'z'")

        positions = numpy.atleast_1d(positions)
        N = len(positions)
        tmp = numpy.zeros((N, 3), self.dtype)
        if direction == "x":
            tmp[:, 0] = positions
            normal = [1.0, 0.0, 0.0]
        elif direction == "y":
            tmp[:, 1] = positions
            normal = [0.0, 1.0, 0.0]
        elif direction == "z":
            tmp[:, 2] = positions
            normal = [0.0, 0.0, 1.0]

        # for regular slices, we dont use the direction vector to pick a projection direction
        sliceDir = [1.0, 0.0, 0.0]
        useDir = False

        for i in range(len(positions)):
            # It is important to ensure each slice get a unique
            # name...so we will number sequentially from python
            j = self.nSlice + i + 1
            sliceName = f"Slice_{j:04d} {groupName} {sliceType.capitalize()} Regular {direction} = {positions[i]:.8f}"

            if sliceType == "relative":
                self.adflow.tecplotio.addparaslice(sliceName, tmp[i], normal, sliceDir, useDir, famList)
            else:
                self.adflow.tecplotio.addabsslice(sliceName, tmp[i], normal, sliceDir, useDir, famList)

        self.nSlice += N

    def addArbitrarySlices(self, normals, points, sliceType="relative", groupName=None, sliceDir=None):
        """
        Add slices that vary arbitrarily in space. This is a generalization
        of the :meth:`addSlices <.addSlices>` routine, where we have the user specify a list of "normals"
        and "points" that define slice planes. Rest of the code is the same.
        This way users can add slices that follow the dihedral of a wing for example.

        Parameters
        ----------
        normals : 3-d array or ndarray (nSlice, 3)
            The normals of the slice directions. If an array of size 3 is passed,
            we use the same normal for all slices. if an array of multiple normals
            are passed, we use the individual normals for each point. in this case,
            the numbers of points and normals must match.
        points : ndarray (nSlice, 3)
            Point coordinates that define a slicing plane along with the normals
        sliceDir : ndarray (nSlice, 3)
            If this is provided, only the intersections that is in this direction
            starting from the normal is kept. This is useful if you are slicing a
            closed geometry and only want to keep one side of it. It is similar to
            the functionality provided in :meth:`addCylindricalSlices <.addCylindricalSlices>`.
        sliceType : str {'relative', 'absolute'}
            Relative slices are 'sliced' at the beginning and then parametrically
            move as the geometry deforms. As a result, the slice through the
            geometry may not remain planar. An absolute slice is re-sliced for
            every output so is always exactly planar and always at the initial
            position the user indicated.
        groupName : str
             The family to use for the slices. Default is None corresponding to all
             wall groups.
        """

        # call the preprocessing routine that avoids some code duplication across 3 types of slices
        sliceType, famList = self._preprocessSliceInput(sliceType, groupName)

        normals = numpy.atleast_2d(normals)
        points = numpy.atleast_2d(points)
        nSlice = len(points)

        if len(normals) == 1:
            tmp = numpy.zeros((nSlice, 3), self.dtype)
            tmp[:] = normals
            normals = tmp

        if sliceDir is None:
            # if we dont have a direction vector to pick a projection direction, we can just set this
            # array to an arbitrary direction. Because the useDir flag is false, the code won't actually
            # use this array
            dummySliceDir = [1.0, 0.0, 0.0]
            useDir = False
        else:
            useDir = True

        for i in range(nSlice):
            # It is important to ensure each slice get a unique
            # name...so we will number sequentially from python
            j = self.nSlice + i + 1

            if useDir:
                direction = sliceDir[i]
            else:
                direction = dummySliceDir

            sliceName = (
                f"Slice_{j:04d} {groupName} {sliceType.capitalize()} Arbitrary "
                f"Point = ({points[i, 0]:.8f}, {points[i, 1]:.8f}, {points[i, 2]:.8f}) "
                f"Normal = ({normals[i, 0]:.8f}, {normals[i, 1]:.8f}, {normals[i, 2]:.8f})"
            )
            if sliceType == "relative":
                self.adflow.tecplotio.addparaslice(sliceName, points[i], normals[i], direction, useDir, famList)
            else:
                self.adflow.tecplotio.addabsslice(sliceName, points[i], normals[i], direction, useDir, famList)

        self.nSlice += nSlice

    def addCylindricalSlices(
        self, pt1, pt2, nSlice=25, sliceBeg=0.0, sliceEnd=360.0, sliceType="relative", groupName=None
    ):
        """
        Add cylindrical projection slices. The cylindrical projection axis is defined
        from pt1 to pt2. Then we start from the lift direction and rotate around this
        axis until we are back at the beginning. The rotation direction follows the
        right hand rule; we rotate the plane in the direction of vectors from pt1 to
        pt2 crossed with the lift direction. This routine will not work as is if the
        cylinder axis is exactly aligned with the lift direction. If that use case is
        needed, consider using :meth:`addArbitrarySlices <.addArbitrarySlices>`.

        Parameters
        ----------
        pt1 : numpy array, length 3
            beginning point of the vector that defines the rotation axis
        pt2 : numpy array, length 3
            end point of the vector that defines the rotation axis
        nSlice : int
            number of slices around a 360 degree full rotation.
        sliceBeg : float
            beginning of the slices in the rotation direction in degrees
        sliceEnd : float
            end of the slices in the rotation direction in degrees
        sliceType : str {'relative', 'absolute'}
            Relative slices are 'sliced' at the beginning and then parametrically
            move as the geometry deforms. As a result, the slice through the
            geometry may not remain planar. An absolute slice is re-sliced for
            every output so is always exactly planar and always at the initial
            position the user indicated.
        groupName : str
             The family to use for the slices. Default is None corresponding to all
             wall groups.
        """

        # call the preprocessing routine that avoids some code duplication across 3 types of slices
        sliceType, famList = self._preprocessSliceInput(sliceType, groupName)

        # get the angles that we are slicing in radians
        angles = numpy.linspace(sliceBeg, sliceEnd, nSlice) * numpy.pi / 180.0

        # unit vector in the lift direction
        liftDir = numpy.zeros(3)
        liftDir[self.getOption("liftIndex") - 1] = 1.0

        # the unit vector on the axis
        vec1 = pt2 - pt1
        vec1 /= numpy.linalg.norm(vec1)

        # loop over angles
        for ii, angle in enumerate(angles):
            sin = numpy.sin(angle)
            cos = numpy.cos(angle)
            omc = 1.0 - cos
            vx = vec1[0]
            vy = vec1[1]
            vz = vec1[2]
            # rotation matrix by theta about vec1
            rot_mat = numpy.array(
                [
                    [cos + vx * vx * omc, vx * vy * omc - vz * sin, vx * vz * omc + vy * sin],
                    [vy * vx * omc + vz * sin, cos + vy * vy * omc, vy * vz * omc - vx * sin],
                    [vz * vx * omc - vy * sin, vz * vy * omc + vx * sin, cos + vz * vz * omc],
                ]
            )

            # rotate the lift direction vector to get the slice direction
            sliceDir = rot_mat.dot(liftDir)

            # In general, the "sliceDir" does not need to be orthogonal to vec1, which is the center axis of the cylinder. So we make it orthogonal
            sliceDir -= vec1.dot(sliceDir) * vec1

            # normalize sliceDir
            sliceDir /= numpy.linalg.norm(sliceDir)

            # cross product vec1 and vec2 to get the normal vector of this plane
            sliceNormal = numpy.cross(vec1, sliceDir)
            # normalize for good measure
            sliceNormal /= numpy.linalg.norm(sliceNormal)

            # we can now define a plane for the slice using p1 and sliceNormal.
            # in the cylindrical version, we also pick a "direction". This is which direction we are going
            # from the center axis. We will reject cells directly if their centroid are behind this
            # direction (i.e. <pt1, centroid> dot sliceDir is negative) and only use the cells that are in
            # the same direction (i.e. dot product positive.)
            # this flag determines that we will use the direction vector as well to pick the slice dir.
            useDir = True

            # It is important to ensure each slice get a unique
            # name...so we will number sequentially from python
            jj = self.nSlice + ii + 1

            sliceName = (
                f"Slice_{jj:04d} {groupName} {sliceType.capitalize()} Cylindrical "
                f"Point = ({pt1[0]:.8f}, {pt1[1]:.8f}, {pt1[2]:.8f}) "
                f"Normal = ({sliceNormal[0]:.8f}, {sliceNormal[1]:.8f}, {sliceNormal[2]:.8f}) "
                f"Axis = ({vec1[0]:.8f}, {vec1[1]:.8f}, {vec1[2]:.8f}) "
                f"Theta = {angle * 180.0 / numpy.pi:.8f}"
            )
            if sliceType == "relative":
                self.adflow.tecplotio.addparaslice(sliceName, pt1, sliceNormal, sliceDir, useDir, famList)
            else:
                self.adflow.tecplotio.addabsslice(sliceName, pt1, sliceNormal, sliceDir, useDir, famList)

        self.nSlice += nSlice

    def _preprocessSliceInput(self, sliceType, groupName):
        """
        Preprocessing routine that holds some of the duplicated code required for 3 types of slice methods.
        """
        # Create the zipper mesh if not done so
        self._createZipperMesh()

        # Determine the families we want to use
        if groupName is None:
            groupName = self.allWallsGroup
        famList = self._getFamilyList(groupName)

        # check the sliceType
        sliceType = sliceType.lower()
        if sliceType not in ["relative", "absolute"]:
            raise Error("'sliceType' must be 'relative' or 'absolute'.")

        return sliceType, famList

    def addIntegrationSurface(self, fileName, familyName, isInflow=True, coordXfer=None):
        """Add a specific integration surface for performing massflow-like
        computations.

        Parameters
        ----------

        fileName : str
           Surface Plot 3D file (may have multiple zones) defining
           integration surface.

        familyName : str
           User supplied name to use for this family. It should not
           already be a name that is defined by the surfaces of the
           CGNS file.

        isInflow : bool
           Flag to treat momentum forces as if it is an inflow or outflow
           face. Default is True

        coordXfer : function
            A callback function that performs a coordinate transformation
            between the original reference frame and any other reference
            frame relevant to the current CFD case. This allows user to apply
            arbitrary modifications to the loaded plot3d surface. The call
            signature is documented in DVGeometry's :meth:`addPointset <pygeo:pygeo.parameterization.DVGeo.DVGeometry.addPointSet>` method.
        """

        self.hasIntegrationSurfaces = True
        # Check that the family name is not already defined:
        if familyName.lower() in self.families:
            raise Error(
                "Cannot add integration surface with family name '%s'"
                "because the name it already exists." % familyName
            )

        # Need to add an additional family so first figure out what
        # the max family index is:
        maxInd = 0
        for fam in self.families:
            if len(self.families[fam]) > 0:
                maxInd = max(maxInd, numpy.max(self.families[fam]))
        famID = maxInd + 1
        self.families[familyName.lower()] = [famID]

        # We are going to operate under the assuption that the number
        # of nodes/elements of the defined surface are sufficiently
        # small that we do not have to worry about parallelization.

        pts, conn = self._readPlot3DSurfFile(fileName, coordXfer=coordXfer)

        self.adflow.usersurfaceintegrations.addintegrationsurface(pts.T, conn.T, familyName, famID, isInflow)

    def addActuatorRegion(
        self,
        fileName,
        axis1,
        axis2,
        familyName,
        thrust=0.0,
        torque=0.0,
        heat=0.0,
        relaxStart=None,
        relaxEnd=None,
        coordXfer=None,
    ):
        """
        Add an actuator disk zone defined by the supplied closed surface in the
        plot3d file "fileName". This surface defines the physical extent of the
        region over which to apply the source terms. The plot3d file may be
        multi-block but all the surface normals must point outside and no
        additional surfaces can be inside. Internally, we find all of the CFD
        volume cells that have their center inside this closed surface and
        apply the source terms over these cells. This surface is only used with
        the original mesh coordinates to mark the internal CFD cells, and we
        keep using these cells even if the geometric design and the mesh
        coordinates change. For example, the marked cells can lie outside of
        the original closed surface after a large design change, but we will
        still use the cells that were inside the surface with the baseline
        design.

        axis1 and axis2 define the vector that we use to determine the
        direction of the thrust addition. Internally, we compute a
        vector by $axisVec = axis2-axis1$ and then we normalize this
        vector. When adding the thrust terms, the direction of the thrust
        is obtained by multiplying the thrust magnitude by this vector.

        Optionally, the source terms in the actuator zone can be
        gradually ramped up as the solution converges. This continuation
        approach can be more robust but the users should be careful with
        the side effects of partially converged solutions. This behavior
        can be controlled by relaxStart and relaxEnd parameters. By
        default, the full magnitude of the source terms are added.
        relaxStart controls when the continuation starts ramping up.
        The value represents the relative convergence in a log10 basis.
        So relaxStart = 2 means the source terms will be inactive until
        the residual is decreased by a factor of 100 compared to free
        stream conditions. The source terms are ramped to the full
        magnitude at relaxEnd. E.g., a relaxEnd value of 4 would
        result in the full source terms after a relative reduction of
        1e4 in the total residuals. If relaxStart is not provided, but
        relaxEnd is provided, then the relaxStart is assumed to be 0.
        If both are not provided, we do not do ramping and just activate
        the full source terms from the beginning. When this continuation
        is used, we internally ramp up the magnitude of the source terms
        monotonically to prevent oscillations; i.e., decrease in the total
        residuals increase the source term magnitudes, but an increase
        in the residuals do not reduce the source terms back down.

        Parameters
        ----------

        fileName : str
           Surface Plot 3D file (multiblock ascii) defining the closed
           region over which the integration is to be applied.

        axis1 : numpy array, length 3
           The physical location of the start of the axis

        axis2 : numpy array, length 4
           The physical location of the end of the axis

        familyName : str
           The name to be associated with the functions defined on this
           region.

        thrust : scalar, float
           The total amount of axial force to apply to this region, in the direction
           of axis1 -> axis2

        torque : scalar, float
           The total amount of torque to apply to the region, about the
           specified axis.

        heat : scalar, float
            The total amount of head added in the actuator zone with source terms

        relaxStart : scalar, float
            The start of the relaxation in terms of
            orders of magnitude of relative convergence

        relaxEnd : scalar, float
            The end of the relaxation in terms of
            orders of magnitude of relative convergence

        coordXfer : function
            A callback function that performs a coordinate transformation
            between the original reference frame and any other reference
            frame relevant to the current CFD case. This allows user to apply
            arbitrary modifications to the loaded plot3d surface. The call
            signature is documented in DVGeometry's :meth:`addPointset <pygeo:pygeo.parameterization.DVGeo.DVGeometry.addPointSet>` method.

        """
        # ActuatorDiskRegions cannot be used in timeSpectralMode
        if self.getOption("equationMode").lower() == "time spectral":
            raise Error("ActuatorRegions cannot be used in Time Spectral Mode.")

        # Load in the user supplied surface file.
        pts, conn = self._readPlot3DSurfFile(fileName, convertToTris=False, coordXfer=coordXfer)

        # We should do a topology check to ensure that the surface the
        # user supplied is actually closed.
        # self._closedTopoCheck(pts, conn)

        # Check if the family name given is already used for something
        # else.
        if familyName.lower() in self.families:
            raise Error(
                "Cannot add ActuatorDiskRegion with family name '%s'" "because the name it already exists." % familyName
            )

        # Need to add an additional family so first figure out what
        # the max family index is:
        maxInd = 0

        for fam in self.families:
            if len(self.families[fam]) > 0:
                maxInd = max(maxInd, numpy.max(self.families[fam]))
        famID = maxInd + 1
        self.families[familyName.lower()] = [famID]

        if relaxStart is None and relaxEnd is None:
            # No relaxation at all
            relaxStart = -1.0
            relaxEnd = -1.0

        if relaxStart is None and relaxEnd is not None:
            # Start at 0 orders if start is not given
            relaxStart = 0.0

        if relaxEnd is None and relaxStart is not None:
            raise Error("relaxEnd must be given is relaxStart is specified")

        # Now continue to fortran were we setup the actual actuator region.
        self.adflow.actuatorregion.addactuatorregion(
            pts.T, conn.T, axis1, axis2, familyName, famID, thrust, torque, heat, relaxStart, relaxEnd
        )

    def writeActuatorRegions(self, fileName, outputDir=None):
        """
        Debug method that writes the cells included in actuator regions
        to a tecplot file. This routine should be called on all of the
        procs in self.comm. The output can then be used to verify that
        the actuator zones are defined correctly.

        Parameters
        ----------
        fileName : str
            Name of the output file

        outputDir : str
            output directory. If not provided, defaults to the output
            directory defined with the aero_options
        """

        if outputDir is None:
            outputDir = self.getOption("outputDirectory")

        # Join to get the actual filename root
        fileName = os.path.join(outputDir, fileName)

        # Ensure extension is .plt even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".plt"

        # just call the underlying fortran routine
        self.adflow.actuatorregion.writeactuatorregions(fileName)

    def addUserFunction(self, funcName, functions, callBack):
        """Add a new function to ADflow by combining existing functions in a
        user-supplied way. The allows the user to define a function
        such as L/D while only requiring a single adjoint solution.

        Parameters
        ----------
        funcName : str
            The name of user-supplied function
        functions : list of strings
            The exisitng function strings which are the input to the
            callBack function
        callBack : python function handle
            The user supplied function that will produce the user function.
            This routine must be complex-step safe for the purpose of
            computing sensitivities.

        Examples
        --------
        >>> ap = AeroProblem(....,evalFuncs=['L/D'])
        >>> CFDSolver = ADFLOW(options=....)
        >>> def LoverD(funcs):
                funcs['L/D'] = funcs['cl']/funcs['cd']
                return funcs
        >>> CFDSolver.addUserFunction('L/D', ['cl','cd'], LoverD)
        """
        funcName = funcName.lower()

        # First check if the function name supplied is already used:
        if funcName in self.adflowCostFunctions:
            raise Error("The supplied funcName %s, is already used" % funcName)

        # Now check that we've supplied at least 2 or more functions
        if len(functions) == 1:
            raise Error("addUserFunction requries at least 2 existing functions" " to be provided")
        # Check that each of the functions are valid:
        for func in functions:
            if func not in self.adflowCostFunctions:
                raise Error("Supplied  function %s to addUserFunction " "not know to ADflow" % func)
        self.adflowUserCostFunctions[funcName] = adflowUserFunc(funcName, functions, callBack)

        return funcName

    def addFunction(self, funcName, groupName, name=None):
        """Add a "new" function to ADflow by restricting the integration of an
        existing ADflow function by a section of the mesh defined by
        'groupName'. The function will be named 'funcName_groupName'
        provided that the 'name' keyword argument is not given. It is
        is, the function will be use that name. If necessary, this
        routine may be used to "change the name" of a function. For
        example,

        >>> addFunction('cd', None, 'super_cd')

        will add a function that is the same as 'cd', but called
        'super_cd'.

        Parameters
        ----------
        funcName : str or list
            The name of the built-in adflow function
        groupName : str or list
            The family or group of families to use for the function
        Name : str or list
            An overwrite name.

        Returns
        -------
        names : list
           The names if the functions that were added

        """

        # Just call the addFunctions() routine with lists
        return self.addFunctions([funcName], [groupName], [name])[0]

    def addFunctions(self, funcNames, groupNames, names=None):
        """Add a series of new functions to ADflow. This is a vector version of
        the addFunction() routine. See that routine for more documentation."""

        if names is None:
            names = [None] * len(funcNames)

        if len(funcNames) != len(groupNames) or len(funcNames) != len(names):
            raise Error("funcNames, groupNames, and names all have to be " "lists of the same length")

        newFuncNames = []
        for i in range(len(funcNames)):
            funcName = funcNames[i]
            groupName = groupNames[i]
            name = names[i]

            # First make sure the supplied function is already known to adflow
            if funcName.lower() not in self.basicCostFunctions:
                raise Error("Supplied function name %s is not known to ADflow" % funcName.lower())

            # Check if the goupName has been added to the mesh.
            if groupName is not None:
                if groupName not in self.families:
                    raise Error(
                        "'%s' is not a family in the CGNS file or has not been added"
                        " as a combination of families" % (groupName)
                    )

            if name is None:
                adflowFuncName = "%s_%s" % (funcName, groupName)
            else:
                adflowFuncName = name

            # Now register the function into adflowCostFunctions
            self.adflowCostFunctions[adflowFuncName] = [groupName, funcName]

            newFuncNames.append(adflowFuncName)

        return newFuncNames

    def setRotationRate(self, rotCenter, rotRate, cgnsBlocks=None):
        """
        Set the rotational rate for the grid:

        rotCenter: 3-vectorThe center of rotation for steady motion
        rotRate: 3-vector or aeroProblem: If it is a 3-vector, take
        the rotations to about x-y-z, if it is an aeroProblem with
        p,q,r, defined, use that to compute rotations.  cgnsBlocks:
        The list of blocks to set. NOTE: This must be in 1-based ordering!
        """
        if isinstance(rotRate, AeroProblem):
            aeroProblem = rotRate
            a = numpy.sqrt(
                self.adflow.flowvarrefstate.gammainf
                * self.adflow.flowvarrefstate.pinfdim
                / self.adflow.flowvarrefstate.rhoinfdim
            )
            V = (self.adflow.inputphysics.machgrid + self.adflow.inputphysics.mach) * a

            p = aeroProblem.phat * V / aeroProblem.spanRef
            q = aeroProblem.qhat * V / aeroProblem.chordRef
            r = aeroProblem.rhat * V / aeroProblem.spanRef
            if aeroProblem.liftIndex == 2:
                rotations = [p, r, q]
            elif aeroProblem.liftIndex == 3:
                rotations = [p, q, r]
            else:
                raise Error(
                    "Invalid lift direction. Must be 2 or 3 for " "steady rotations and specifying an aeroProblem"
                )

        else:
            rotations = rotRate

        if cgnsBlocks is None:
            # By Default set all the blocks:
            cgnsBlocks = numpy.arange(1, self.adflow.cgnsgrid.cgnsndom + 1)

        self.adflow.preprocessingapi.updaterotationrate(rotCenter, rotations, cgnsBlocks)
        self._updateVelInfo = True

    def checkPartitioning(self, nprocs):
        """This function determines the potential load balancing for
        nprocs. The intent is this function can be run in serial
        to determine the best number of procs for load balancing. The
        grid is never actually loaded so this function can be run with
        VERY large grids without issue.

        Parameters
        ----------
        nProcs : int
            The number of procs check partitioning on

        Returns
        -------
        loadInbalance : float
            The fraction of inbalance for cells. 0 is good. 1.0 is really bad
        faceInbalance : float
            The fraction of inbalance for faces. This gives an idea of
            communication code. 0 is god. 1.0 is really bad.
        """
        loadInbalance, faceInbalance = self.adflow.partitioning.checkpartitioning(nprocs)

        return loadInbalance, faceInbalance

    def __call__(self, aeroProblem, **kwargs):
        """
        Common routine for solving the problem specified in
        aeroProblem.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The complete description of the problem to solve
        mdCallBack : python function
            Call back for doing unsteady aeroelastic. Currently
            not supported.
        writeSolution : bool
            Flag to override any solution writing parameters. This
            is used in a multidisciplinary environment when the outer
            solver can suppress all I/O during intermediate solves.
        """
        self.rootChangedOptions = self.comm.bcast(self.rootChangedOptions, root=0)
        for key, val in self.rootChangedOptions.items():
            self.setOption(key, val)
        startCallTime = time.time()

        # Get option about adjoint memory
        releaseAdjointMemory = kwargs.pop("relaseAdjointMemory", True)

        # Set the aeroProblem
        self.setAeroProblem(aeroProblem, releaseAdjointMemory)

        # Set filenames for possible forced write during solution
        self._setForcedFileNames()

        # Possibly release adjoint memory if not already done
        # so. Don't remove this. If switching problems, this is done
        # in setAeroProblem, but for when not switching it must be
        # done here regardless.
        if releaseAdjointMemory:
            self.releaseAdjointMemory()

        # Clear out any saved adjoint RHS since they are now out of
        # data. Also increment the counter for this case.
        self.curAP.adflowData.adjointRHS = OrderedDict()
        self.curAP.adflowData.callCounter += 1

        # --------------------------------------------------------------
        # Setup iteration arrays ---- don't touch this unless you
        # REALLY REALLY know what you're doing!

        # iterTot is non-zero which means we already did a solution,
        # so we can run on the finest level.
        if self.adflow.iteration.itertot != 0:
            self.adflow.inputiteration.mgstartlevel = 1

        self.adflow.iteration.itertot = 0
        desiredSize = self.adflow.inputiteration.nsgstartup + self.adflow.inputiteration.ncycles
        self.adflow.utils.allocconvarrays(desiredSize)
        # --------------------------------------------------------------

        if self.getOption("equationMode").lower() == "unsteady":
            self.adflow.utils.alloctimearrays(self.getOption("nTimeStepsFine"))
            self._setUnsteadyFileParameters()

        # Mesh warp may have already failed:
        if not self.adflow.killsignals.fatalfail:
            if (
                self.getOption("equationMode").lower() == "steady"
                or self.getOption("equationMode").lower() == "time spectral"
            ):
                self.updateGeometryInfo()

                # Check to see if the above update routines failed.
                self.adflow.killsignals.fatalfail = self.comm.allreduce(
                    bool(self.adflow.killsignals.fatalfail), op=MPI.LOR
                )

        if self.adflow.killsignals.fatalfail:
            numDigits = self.getOption("writeSolutionDigits")
            fileName = f"failed_mesh_{self.curAP.name}_{self.curAP.adflowData.callCounter:0{numDigits}d}.cgns"
            self.pp(f"Fatal failure during mesh warp! Bad mesh is written in output directory as {fileName}")
            self.writeMeshFile(os.path.join(self.getOption("outputDirectory"), fileName))
            self.curAP.fatalFail = True
            self.curAP.solveFailed = True
            return

        # We now know the mesh warping was ok so reset the flags:
        self.adflow.killsignals.routinefailed = False
        self.adflow.killsignals.fatalfail = False
        sys.stdout.flush()
        t1 = time.time()  # Signal check time

        # Solve the equations in appropriate mode
        mode = self.getOption("equationMode").lower()
        if mode in ["steady", "time spectral"]:
            # In steady mode, the wall temperature has to be set after the solver
            # switches to current AP, otherwise the data will be overwritten.
            wallTemperature = kwargs.pop("wallTemperature", None)
            groupName = kwargs.pop("groupName", None)
            if wallTemperature is not None:
                self.setWallTemperature(wallTemperature, groupName=groupName)

            self.adflow.solvers.solver()
        elif mode == "unsteady":
            # Initialization specific for unsteady simulation
            self.adflow.solvers.solverunsteadyinit()

            # If ADflow is to be coupled with other solvers, we are done
            if self.getOption("coupledsolution"):
                return
            # Otherwise, the problem is solved within this function
            else:
                # Get the callback functions
                surfaceMeshCallback = kwargs.pop("surfaceMeshCallback", None)
                volumeMeshCallback = kwargs.pop("volumeMeshCallback", None)

                # Call solver wrapper for unsteady simulations
                self._solverunsteady(surfaceMeshCallback=surfaceMeshCallback, volumeMeshCallback=volumeMeshCallback)
        else:
            raise Error("Mode {0} not supported".format(mode))

        # Save the states into the aeroProblem
        self.curAP.adflowData.stateInfo = self._getInfo()

        # Save the winf vector. (just the flow part)
        self.curAP.adflowData.oldWinf = self.adflow.flowvarrefstate.winf[0:5].copy()

        # Save the current ANK CFL in the ap data
        self.curAP.adflowData.ank_cfl = self.adflow.anksolver.ank_cfl

        # Assign Fail Flags
        self.curAP.solveFailed = bool(self.adflow.killsignals.routinefailed)
        self.curAP.fatalFail = bool(self.adflow.killsignals.fatalfail)

        # Reset Flow if there's a fatal fail reset and return;
        # --> Do not write solution
        if self.curAP.fatalFail:
            self.resetFlow(aeroProblem)
            return

        t2 = time.time()
        solTime = t2 - t1

        # Post-Processing
        # --------------------------------------------------------------
        # Solution is written if writeSolution argument true or does not exist
        if kwargs.pop("writeSolution", True):
            self.writeSolution()

        writeSolutionTime = time.time()

        if self.getOption("TSStability"):
            self.computeStabilityParameters()  # should this be in evalFuncs?

        stabilityParameterTime = time.time()

        if self.getOption("printTiming") and self.comm.rank == 0:
            print("+-------------------------------------------------+")
            print("|")
            print("| Solution Timings:")
            print("|")
            print("| %-30s: %10.3f sec" % ("Set AeroProblem Time", t1 - startCallTime))
            print("| %-30s: %10.3f sec" % ("Solution Time", solTime))
            print("| %-30s: %10.3f sec" % ("Write Solution Time", writeSolutionTime - t2))
            print("| %-30s: %10.3f sec" % ("Stability Parameter Time", stabilityParameterTime - writeSolutionTime))
            print("|")
            print("| %-30s: %10.3f sec" % ("Total Call Time", stabilityParameterTime - startCallTime))
            print("+--------------------------------------------------+")

    def _solverunsteady(self, surfaceMeshCallback=None, volumeMeshCallback=None):
        """
        This routine is intended for time accurate simulation, including
        - Unsteady problem with (possibly) prescribed mesh motion
        - Unsteady problem with warping mesh
        - Unsteady problem with rigidly translating/rotating mesh (N/A yet)
        The latter two require ALE.

        Parameters
        ----------
        surfaceMeshCallback : Python function handle
            Computes new surface coordinates for mesh warping

        volumeMeshCallback : Python function handle
            Computes new volume coordinates for mesh moving
        """

        # Store initial data if necessary
        if surfaceMeshCallback is not None:
            refCoor = self.getSurfaceCoordinates(self.designFamilyGroup)
        else:
            refCoor = None
        if volumeMeshCallback is not None:
            refGrid = self.mesh.getSolverGrid().reshape(-1, 3)
        else:
            refGrid = None

        # Time advancing
        nTimeStepsFine = self.getOption("ntimestepsfine")
        for _tdx in range(1, nTimeStepsFine + 1):
            # Increment counter
            curTime, curTimeStep = self.advanceTimeStepCounter()

            # Warp mesh if necessary
            if surfaceMeshCallback is not None:
                newCoor = surfaceMeshCallback(refCoor, curTime, curTimeStep)
                self.setSurfaceCoordinates(newCoor, self.designFamilyGroup)
                self.updateGeometryInfo()
            # Rigidly move mesh if necessary
            elif volumeMeshCallback is not None:
                newGrid = volumeMeshCallback(refGrid, curTime, curTimeStep).reshape(-1)
                self.adflow.warping.setgrid(newGrid)
                self._updateGeomInfo = True
                self.updateGeometryInfo(warpMesh=False)
            # Otherwise do naive unsteady simulation
            else:
                self.adflow.preprocessingapi.shiftcoorandvolumes()
                self.adflow.solvers.updateunsteadygeometry()

            # Solve current step
            self.solveTimeStep()
            self.writeSolution()

    def getConvergenceHistory(self, workUnitTime=None):
        """
        Retrieve the convergence history from the fortran level.

        This information is printed to the terminal during a run.
        It is iterTot, IterType, CFL, Step, Linear Res, and CPU time (if added as a monitor variable)
        and the data for each monitor variable.

        Parameters
        ----------
        workUnitTime : float
            The scaling factor specific to the processor.
            If provided and CPU time is a monitor variable (`showcpu` is true), the work units
            (a processor independent time unit) ~will be added to the returned dict too.

        Returns
        -------
        convergeDict : dict
            A dictionary of arrays and lists. The keys are the data types.
            The indices of the arrays are the major iteration numbers.

        Examples
        --------
        >>> CFDsolver(ap)
        >>> CFDsolver.evalFunctions(ap1, funcs, ['cl', 'cd'])
        >>> hist = CFDSolver.getConvergenceHistory()
        >>> if MPI.COMM_WORLD.rank == 0:
        >>>     with open(os.path.join(output_dir, "convergence_history.pkl"), "wb") as f:
        >>>         pickle.dump(hist, f)
        """

        # only the root proc has all the data and the data should be global.
        # so get the data from the root proc and then broadcast it.
        if self.comm.rank == 0:
            # --- get data ---
            convergeArray = self.adflow.monitor.convarray
            solverDataArray = self.adflow.monitor.solverdataarray

            # We can't access monnames directly, because f2py apparently can't
            # handle an allocatable array of strings.
            # related stackoverflow posts:
            # https://stackoverflow.com/questions/60141710/f2py-on-module-with-allocatable-string
            # https://stackoverflow.com/questions/64900066/dealing-with-character-arrays-in-f2py
            monNames = self.adflow.utils.getmonitorvariablenames(int(self.adflow.monitor.nmon))

            # similarly we need to use a special function to
            # retrive the additional solver type information
            typeArray = self.adflow.utils.getsolvertypearray(
                self.adflow.iteration.itertot, self.adflow.inputtimespectral.ntimeintervalsspectral
            )

            # --- format the data ---
            # trim the array
            solverDataArray = self._trimHistoryData(solverDataArray)
            convergeArray = self._trimHistoryData(convergeArray)
            monNames = self._convertFortranStringArrayToList(monNames)

            # convert the fortran sting array to a list of strings
            type_list = []
            for idx_sps in range(self.adflow.inputtimespectral.ntimeintervalsspectral):
                type_list.append(self._convertFortranStringArrayToList(typeArray[:, idx_sps, :]))

            # there is only one time spectral instance, so that dimension can be removed
            if self.adflow.inputtimespectral.ntimeintervalsspectral == 1:
                type_list = type_list[0]

            # --- create a dictionary with the labels and the values ---
            convergeDict = CaseInsensitiveDict(
                {
                    "total minor iters": numpy.array(solverDataArray[:, 0], dtype=int),
                    "CFL": solverDataArray[:, 1],
                    "step": solverDataArray[:, 2],
                    "linear res": solverDataArray[:, 3],
                    "iter type": type_list,
                }
            )

            if self.adflow.monitor.showcpu:
                convergeDict["wall time"] = solverDataArray[:, 4]
                if workUnitTime is not None:
                    convergeDict["work units"] = convergeDict["wall time"] * (self.comm.size) / workUnitTime

            for idx, var in enumerate(monNames):
                convergeDict[var] = convergeArray[:, idx]
        else:
            convergeDict = {}

        convergeDict = self.comm.bcast(convergeDict, root=0)

        return convergeDict

    def _trimHistoryData(self, data_array):
        """
        Trim the convergence history data to the number of iterations actually run.
        and remove the extra dimensions if only one time spectral instance was used.
        """

        # --- trim the array ---
        # how many iterations were actually run
        num_iters = self.adflow.iteration.itertot
        # add an iteration for the initial residual evaluation
        data_array = data_array[: (num_iters + 1)]

        # remove possibly unnecessary dimensions
        # if there is only one grid in the convergence history (no multigrid)
        if self.adflow.inputtimespectral.ntimeintervalsspectral == 1:
            data_array = data_array.reshape((data_array.shape[0], data_array.shape[2]))

        return data_array

    def advanceTimeStepCounter(self):
        """
        Advance one unit of timestep and physical time.
        """
        self.adflow.monitor.timestepunsteady = self.adflow.monitor.timestepunsteady + 1
        self.adflow.monitor.timeunsteady = self.adflow.monitor.timeunsteady + self.adflow.inputunsteady.deltat

        if self.myid == 0:
            self.adflow.utils.unsteadyheader()

        return self.adflow.monitor.timeunsteady, self.adflow.monitor.timestepunsteady

    def solveTimeStep(self):
        """
        Solve the current time step, and write solutions if necessary
        """
        t1 = time.time()
        self.adflow.solvers.solverunsteadystep()
        t2 = time.time()
        if self.getOption("printTiming") and self.comm.rank == 0:
            print("Timestep solution time: {0:10.3f} sec".format(t2 - t1))

    def evalFunctions(self, aeroProblem, funcs, evalFuncs=None, ignoreMissing=False):
        """
        Evaluate the desired functions given in iterable object,
        'evalFuncs' and add them to the dictionary 'funcs'. The keys
        in the funcs dictionary will be have an ``<ap.name>_`` prepended to
        them.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aerodynamic problem to to get the solution for

        funcs : dict
            Dictionary into which the functions are saved.

        evalFuncs : iterable object containing strings
          If not None, use these functions to evaluate.

        ignoreMissing : bool
            Flag to suppress checking for a valid function. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> CFDsolver(ap)
        >>> CFDsolver.evalFunctions(ap1, funcs, ['cl', 'cd'])
        >>> funcs
        >>> # Result will look like (if aeroProblem, ap1, has name of 'wing'):
        >>> # {'wing_cl':0.501, 'wing_cd':0.02750}
        """

        startEvalTime = time.time()

        # Set the AP
        self.setAeroProblem(aeroProblem)

        aeroProblemTime = time.time()

        if evalFuncs is None:
            evalFuncs = sorted(self.curAP.evalFuncs)

        # Make sure we have a list that has only lower-cased entries
        tmp = []
        for f in evalFuncs:
            tmp.append(f.lower())
        evalFuncs = sorted(tmp)

        # We need to determine how many different groupings we have,
        # since we will have to call getSolution for each *unique*
        # function grouping. We can also do the error checking
        groupMap = OrderedDict()

        def addToGroup(f):
            group = self.adflowCostFunctions[f][0]
            basicFunc = self.adflowCostFunctions[f][1]
            if group not in groupMap:
                groupMap[group] = [[basicFunc, f]]
            else:
                groupMap[group].append([basicFunc, f])

        callBackFuncs = OrderedDict()
        for f in evalFuncs:
            # Standard functions
            if f in self.adflowCostFunctions:
                addToGroup(f)

            # User supplied functions
            elif f in self.adflowUserCostFunctions:
                for sf in self.adflowUserCostFunctions[f].functions:
                    addToGroup(sf)
                    # This just flags the sub-function 'sf' as being
                    # required for callBacks.
                    callBackFuncs[sf] = 0.0
            else:
                if not ignoreMissing:
                    raise Error("Supplied function %s is not known to ADflow." % f)

        # Now we loop over the unique groups calling the required
        # getSolution, there may be just one. No need for error
        # checking here.
        for group in groupMap:
            # Call the lower level getSolution() function
            res = self.getSolution(groupName=group)

            # g contains the "basic function" (index 0) and the actual
            # function name (index 1)
            for g in groupMap[group]:
                key = self.curAP.name + "_%s" % g[1]
                self.curAP.funcNames[g[1]] = key
                if g[1].lower() in evalFuncs:
                    funcs[key] = res[g[0]]

                if g[1].lower() in callBackFuncs:
                    callBackFuncs[g[1]] = res[g[0]]

        getSolutionTime = time.time()

        # Execute the user supplied functions if there are any
        for f in self.adflowUserCostFunctions:
            if f in evalFuncs:
                self.adflowUserCostFunctions[f].evalFunctions(callBackFuncs)
                key = self.adflowUserCostFunctions[f].funcName
                value = callBackFuncs[key]
                funcs[self.curAP.name + "_%s" % key] = value

        userFuncTime = time.time()
        if self.getOption("printTiming") and self.comm.rank == 0:
            print("+---------------------------------------------------+")
            print("|")
            print("| Function Timings:")
            print("|")
            print("| %-30s: %10.3f sec" % ("Function AeroProblem Time", aeroProblemTime - startEvalTime))
            print("| %-30s: %10.3f sec" % ("Function Evaluation Time", getSolutionTime - aeroProblemTime))
            print("| %-30s: %10.3f sec" % ("User Function Evaluation Time", userFuncTime - getSolutionTime))
            print("|")
            print("| %-30s: %10.3f sec" % ("Total Function Evaluation Time", userFuncTime - startEvalTime))
            print("+--------------------------------------------------+")

    def _getFuncsBar(self, f):
        # Internal routine to return the funcsBar dictionary for the
        # given objective. For regular functions this just sets a 1.0
        # for f. For the user-supplied functions, it linearizes the
        # user-supplied function and sets all required seeds.

        if f.lower() in self.adflowCostFunctions:
            funcsBar = {f.lower(): 1.0}
        elif f in self.adflowUserCostFunctions:
            # Need to get the funcs-bar derivative from the
            # user-supplied function
            funcsBar = self.adflowUserCostFunctions[f].evalFunctionsSens()
        else:
            funcsBar = {f.lower(): 0.0}
        return funcsBar

    def evalFunctionsSens(self, aeroProblem, funcsSens, evalFuncs=None):
        """
        Evaluate the sensitivity of the desired functions given in
        iterable object, 'evalFuncs' and add them to the dictionary
        'funcSens'. The keys in the 'funcsSens' dictionary will be have an
        ``<ap.name>_`` prepended to them.

        Parameters
        ----------
        funcSens : dict
            Dictionary into which the function derivatives are saved.

        evalFuncs : iterable object containing strings
            The additional functions the user wants returned that are
            not already defined in the aeroProblem

        Examples
        --------
        >>> funcSens = {}
        >>> CFDsolver.evalFunctionsSens(ap1, funcSens, ['cl', 'cd'])
        """

        # This is the one and only gateway to the getting derivatives
        # out of adflow. If you want a derivative, it should come from
        # here.

        startEvalSensTime = time.time()

        self.setAeroProblem(aeroProblem)

        # we now need to reset the adjoint fail flag. This may be the
        # first time we are setting this for this AP, or we may have
        # a fail flag remaining from a previous call. Either way, we
        # want to at least try solving the adjoints here, regardless
        # of what happened in the previous iterations.
        self.curAP.adjointFailed = False

        if evalFuncs is None:
            evalFuncs = sorted(self.curAP.evalFuncs)

        # Make sure we have a list that has only lower-cased entries
        tmp = []
        for f in evalFuncs:
            tmp.append(f.lower())
        evalFuncs = tmp

        adjointStartTime = {}
        adjointEndTime = {}
        totalSensEndTime = {}

        # Do the functions one at a time:
        for f in evalFuncs:
            if f in self.adflowCostFunctions:
                pass
            elif f in self.adflowUserCostFunctions:
                pass
            else:
                raise Error("Supplied %s function is not known to ADflow." % f)

            adjointStartTime[f] = time.time()

            if self.comm.rank == 0:
                print("Solving adjoint: %s" % f)

            key = self.curAP.name + "_%s" % f

            # Set dict structure for this derivative
            funcsSens[key] = OrderedDict()

            # Solve adjoint equation
            self.solveAdjoint(aeroProblem, f)

            adjointEndTime[f] = time.time()
            # Now, due to the use of the super combined
            # computeJacobianVectorProductBwd() routine, we can complete
            # the total derivative computation in a single call. We
            # simply seed resBar with *NEGATIVE* of the adjoint we
            # just computed along with 1.0 for the objective and then
            # everything completely falls out. This saves doing 4
            # individual calls to: getdRdXvTPsi, getdIdx, getdIda,
            # getdIdaTPsi.

            # These are the reverse mode seeds
            psi = -self.getAdjoint(f)

            # Compute everything and update into the dictionary
            funcsSens[key].update(
                self.computeJacobianVectorProductBwd(resBar=psi, funcsBar=self._getFuncsBar(f), xDvDeriv=True)
            )

            totalSensEndTime[f] = time.time()

        finalEvalSensTime = time.time()

        if self.getOption("printTiming") and self.comm.rank == 0:
            print("+--------------------------------------------------+")
            print("|")
            print("| Adjoint Times:")
            print("|")
            for f in evalFuncs:
                print(
                    "| %-30s: %10.3f sec" % ("Adjoint Solve Time - %s" % (f), adjointEndTime[f] - adjointStartTime[f])
                )
                print(
                    "| %-30s: %10.3f sec"
                    % ("Total Sensitivity Time - %s" % (f), totalSensEndTime[f] - adjointEndTime[f])
                )
            print("|")
            print("| %-30s: %10.3f sec" % ("Complete Sensitivity Time", finalEvalSensTime - startEvalSensTime))
            print("+--------------------------------------------------+")

    def propagateUncertainty(self, aeroProblem, evalFuncs=None, UQDict=None):
        """
        Use the first order second moment method to predict output uncertainties for the current solution.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aerodynamic problem to solve

        evalFuncs : iterable object containing strings
            The functions the user wants the uncertainty of

        UQDict : dict
            Dictionary containing the mean and std dev. of the input parameters that are providing the uncertain input.
        """

        # Set the current aeroProblem
        self.setAeroProblem(aeroProblem)

        # Check for functions to propagate to
        if evalFuncs is None:
            evalFuncs = sorted(self.curAP.evalFuncs)

        # Make sure we have a list that has only lower-cased entries
        tmp = []
        for f in evalFuncs:
            tmp.append(f.lower())
        evalFuncs = tmp

        # get the function values
        funcs = {}
        self.evalFunctions(aeroProblem, funcs, evalFuncs)

        # Start by getting the total derivatives of the outputs of interest
        # with respect to the variables providing the uncertainty.

        derivs = {}
        self.evalFunctionsSens(aeroProblem, derivs, evalFuncs)

        UQOut = {}
        for outKey in evalFuncs:
            dictKey = aeroProblem.name + "_" + outKey
            UQOut[dictKey] = {}

            # compute sigma**2
            sigma2 = 0
            for key in UQDict.keys():
                varKey = key + "_" + aeroProblem.name
                for i in range(UQDict[key]["size"]):
                    if isinstance(derivs[dictKey][varKey], (list, tuple, numpy.ndarray)):
                        dodsigma = derivs[dictKey][varKey][i]
                        sigma = UQDict[key]["sigma"][i]
                    else:
                        dodsigma = derivs[dictKey][varKey]
                        sigma = UQDict[key]["sigma"]

                    sigma2 += (dodsigma * sigma) ** 2

            UQOut[dictKey]["mu"] = funcs[dictKey]
            UQOut[dictKey]["sigma"] = numpy.sqrt(sigma2)

        return UQOut

    def solveCL(
        self,
        aeroProblem,
        CLStar,
        alpha0=None,
        CLalphaGuess=None,
        delta=0.5,
        tol=1e-3,
        autoReset=False,
        maxIter=20,
        relaxCLa=1.0,
        relaxAlpha=1.0,
        L2ConvRel=None,
        stopOnStall=False,
        writeSolution=False,
        workUnitTime=None,
        updateCutoff=1e-16,
    ):
        """
        This is a simple secant method search for solving for a
        fixed CL. This really should only be used to determine the
        starting alpha for a lift constraint in an optimization.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aerodynamic problem to solve
        CLStar : float
            The desired target CL
        alpha0 : angle (deg)
            Initial guess for secant search (deg). If None, use the
            value in the aeroProblem
        CLalphaGuess : float or None
            The user can provide an estimate for the lift curve slope
            in order to accelerate convergence. If the user supply a
            value to this option, it will not use the delta value anymore
            to select the angle of attack of the second run. The value
            should be in 1/deg.
        delta : angle (deg)
            Initial step direction for secant search
        tol : float
            Desired tolerance for CL
        autoReset : bool
            Flag to reset flow between solves. The Euler NK method has
            issues when only the alpha is changed (Martois effect we
            think). This will reset the flow after each solve which
            solves this problem. Not necessary (or desired) when using
            the RK solver. The useCorrection option is the preferred
            way of addressing this issue with the NK/ANK methods.
            If solver still struggles after restarts with the correction,
            this option can be used.
        maxIter : int
            Maximum number of secant solver iterations.
        relaxCLa : float
            Relaxation factor for the CLalpha update. Value of 1.0
            will give the standard secant solver. A lower value will
            damp the update to the CL alpha. Useful when the CL output
            is expected to be noisy.
        relaxAlpha : float
            Relaxation factor for the alpha update. Value of 1.0
            will give the standard secant solver. A lower value will
            damp the alpha update. Useful to prevent large changes to
            alpha.
        L2ConvRel : float or list or None
            Temporary relative L2 convergence for each iteration of the
            CL solver. If the option is set to None, we don't modify the
            L2ConvergenceRel option. If a float is provided, we use this
            value as the target relative L2 convergence for each iteration.
            If a list of floats is provided, we iterate over the values
            and set the relative L2 convergence target separately for
            each iteration. If the list length is less than the maximum
            number of iterations we can perform, we keep using the last
            value of the list until we run out of iterations.
        stopOnStall : bool
            Flag to determine if we want to stop early if the solver
            fails to achieve the relative or total L2 convergence. This
            is useful when the solver is expected to converge to the target
            L2 tolerances at each call.
        writeSolution : bool
            Flag to enable a call to self.writeSolution as the CL solver
            is finished.
        workUnitTime : float or None
            Optional parameter that will be passed to getConvergenceHistory
            after each adflow call.
        updateCutoff : float
            Required change in alpha to trigger the clalpha update with the
            Secant method. If the change in alpha is smaller than this option
            we don't even bother changing clalpha and keep using the clalpha
            value from the last iteration. This prevents clalpha from getting
            bad updates due to noisy cl outputs that result from small alpha
            changes.

        Returns
        -------
        resultsDict : dictionary
            Dictionary that contains various results from the cl solve. The
            dictionary contains the following data:

            converged : bool
                Flag indicating cl solver convergence. True means both the
                CL solver target tolerance AND the overall target L2 convergence
                is achieved. False means either one or both of the tolerances
                are not reached.
            iterations : int
                Number of iterations ran.
            l2convergence : float
                Final relative L2 convergence of the solver w.r.t. free
                stream residual. Useful to check overall CFD convergence.
            alpha : float
                Final angle of attack value
            cl : float
                Final CL value.
            clstar : float
                The original target CL. returned for convenience.
            error : float
                Error in cl value.
            clalpha : float
                Estimate to clalpha used in the last solver iteration
            time : float
                Total time the solver needed, in seconds
            history : list
                List of solver convergence histories. Each entry
                in the list contains the convergence history of the
                solver from each CL Solve iteration. The history
                is returned using the "getConvergenceHistory" method.

        """

        # time the CL solve
        t1 = time.time()

        # pointer to the iteration module for faster access
        iterationModule = self.adflow.iteration

        # create the string template we want to print for each iteration
        iterString = (
            "\n"
            + "+--------------------------------------------------+\n"
            + "|\n"
            + "| CLSolve Iteration   {iIter}\n"
            + "| Elapsed Time        {curTime:.3f} sec\n"
            + "|\n"
            + "| L2 Convergence      {L2Conv}\n"
            + "| L2 Rel Convergence  {L2ConvRel}\n"
            + "| Alpha               {curAlpha}\n"
            + "| CL                  {CL}\n"
            + "| CLStar              {CLStar}\n"
            + "| Error               {err}\n"
            + "| \n"
            + "| CLAlpha             {clalpha}\n"
            + "| Delta Alpha         {deltaAlpha}\n"
            + "| New Alpha           {newAlpha}\n"
            + "+--------------------------------------------------+\n"
        )

        def checkConvergence(curError):
            # check if L2 Convergence is achieved, if not, we cant quit yet
            # this is not the same as checking if AP solve worked. Here, we
            # explicitly want to check if the final L2 target is reached,
            # i.e. the solution is a valid converged state.
            L2Conv = iterationModule.totalrfinal / iterationModule.totalr0
            L2ConvTarget = self.getOption("L2Convergence")

            # we check if the current error is below tolerance, and if we have also reached the L2 target
            if abs(curError) < tol and L2Conv <= L2ConvTarget:
                converged = True
            else:
                converged = False

            return converged

        def finalizeSolver(resultsDict, modifiedOptions):
            # exit the solver. we have a few housekeeping items:

            # put back the modified options
            for option, value in modifiedOptions.items():
                self.setOption(option, value)

            # write solution before returning if requested
            if writeSolution:
                self.writeSolution()

            # final printout
            if self.comm.rank == 0:
                # print all results but the history
                printDict = copy.deepcopy(resultsDict)
                printDict.pop("history")
                print("\n+--------------------------------------------------+")
                print("| CLSolve Results:")
                print("+--------------------------------------------------+")
                # print all but convergence history from every call
                self.pp(printDict)
                print("+--------------------------------------------------+\n", flush=True)

            return

        # list to keep track of convergence history
        convergenceHistory = []

        # Dictionary to keep track of all the options we have modified.
        # We keep the original values in this dictionary to be put back once we are done.
        modifiedOptions = {}

        # check if we have a custom relative l2 convergence.
        if L2ConvRel is not None:
            # save the original value
            modifiedOptions["L2ConvergenceRel"] = self.getOption("L2ConvergenceRel")

            # if a single value is provided, convert it into a list so that we can apply it at each iteration
            if isinstance(L2ConvRel, float):
                L2ConvRel = [L2ConvRel] * maxIter
            elif isinstance(L2ConvRel, list):
                # we may need to extend this list if its shorter than maxIter
                if len(L2ConvRel) < maxIter:
                    # pick the last value and add
                    L2ConvRel += [L2ConvRel[-1]] * (maxIter - len(L2ConvRel))
            else:
                raise Error("L2ConvRel needs to be a float or a list of floats")

            # set the first rel conv value.
            self.setOption("L2ConvergenceRel", L2ConvRel[0])

        # initialize the first alpha guess
        if alpha0 is not None:
            aeroProblem.alpha = alpha0
        else:
            alpha0 = aeroProblem.alpha

        # Set the starting value
        anm2 = aeroProblem.alpha

        # first analysis. We let the call counter get incremented here, but the following
        # calls do not increase adflow's call counter. As a result, the solution files will be
        # consistently named when multiple CL solutions are performed back to back.
        self.__call__(aeroProblem, writeSolution=False)
        convergenceHistory.append(self.getConvergenceHistory(workUnitTime=workUnitTime))
        sol = self.getSolution()
        fnm2 = sol["cl"] - CLStar

        if CLalphaGuess is None:
            # Use the delta option to define the next Aoa
            if fnm2 < 0:
                anm1 = alpha0 + abs(delta)
            else:
                anm1 = alpha0 - abs(delta)
            # set this to zero for the initial printout
            clalpha = 0.0
        else:
            # Use CLalphaGuess option to define the next Aoa
            anm1 = alpha0 - fnm2 / CLalphaGuess
            # initialize the clalpha working variable
            clalpha = CLalphaGuess

        # current L2 convergence is also needed for the print
        L2Conv = iterationModule.totalrfinal / iterationModule.totalr0

        # print iteration info
        if self.comm.rank == 0:
            print(
                iterString.format(
                    iIter=1,
                    curTime=time.time() - t1,
                    L2Conv=L2Conv,
                    L2ConvRel=iterationModule.totalrfinal / iterationModule.totalrstart,
                    curAlpha=aeroProblem.alpha,
                    CL=sol["cl"],
                    CLStar=CLStar,
                    err=fnm2,
                    clalpha=clalpha,
                    deltaAlpha=anm1 - aeroProblem.alpha,
                    newAlpha=anm1,
                )
            )

        # Check for convergence if the starting point is already ok.
        converged = checkConvergence(fnm2)

        # rest of the results
        CL = sol["cl"]
        err = sol["cl"] - CLStar
        t2 = time.time()
        resultsDict = {
            "converged": converged,
            "iterations": 1,
            "l2convergence": L2Conv,
            "alpha": aeroProblem.alpha,
            "cl": CL,
            "clstar": CLStar,
            "error": err,
            "clalpha": clalpha,
            "time": t2 - t1,
            "history": convergenceHistory,
        }

        if converged or (aeroProblem.solveFailed and stopOnStall) or aeroProblem.fatalFail:
            finalizeSolver(resultsDict, modifiedOptions)
            return resultsDict

        # Secant method iterations
        for _iIter in range(1, maxIter):
            if L2ConvRel is not None:
                self.setOption("L2ConvergenceRel", L2ConvRel[_iIter])

            # We may need to reset the flow since changing strictly
            # alpha leads to problems with the NK solver
            if autoReset:
                self.resetFlow(aeroProblem)

            # Set current alpha
            aeroProblem.alpha = anm1

            # Solve for n-1 value (anm1)
            self.__call__(aeroProblem, writeSolution=False)
            convergenceHistory.append(self.getConvergenceHistory(workUnitTime=workUnitTime))

            self.curAP.adflowData.callCounter -= 1
            sol = self.getSolution()
            fnm1 = sol["cl"] - CLStar

            # Secant Update
            if _iIter == 1 and CLalphaGuess is None:
                # we first need to compute the estimate to clalpha
                clalpha = (fnm1 - fnm2) / (anm1 - anm2)
            else:
                # only update clalpha if the change in alpha was greater than updateCutoff
                if numpy.abs(anm1 - anm2) > updateCutoff:
                    # we update the clalpha using a relaxed update
                    claNew = (fnm1 - fnm2) / (anm1 - anm2)
                    clalpha += relaxCLa * (claNew - clalpha)

            # similarly, the user might want to relax the alpha update
            anew = anm1 - (fnm1 / clalpha) * relaxAlpha

            # Shift n-1 values to n-2 values
            fnm2 = fnm1
            anm2 = anm1

            # Se the n-1 alpha value from update
            anm1 = anew

            # current L2 convergence is also needed for the print
            L2Conv = iterationModule.totalrfinal / iterationModule.totalr0

            # print iteration info
            if self.comm.rank == 0:
                print(
                    iterString.format(
                        iIter=_iIter + 1,
                        curTime=time.time() - t1,
                        L2Conv=L2Conv,
                        L2ConvRel=iterationModule.totalrfinal / iterationModule.totalrstart,
                        curAlpha=aeroProblem.alpha,
                        CL=sol["cl"],
                        CLStar=CLStar,
                        err=fnm1,
                        clalpha=clalpha,
                        deltaAlpha=anew - aeroProblem.alpha,
                        newAlpha=anew,
                    )
                )

            # Check for convergence if the starting point is already ok.
            converged = checkConvergence(fnm1)

            # rest of the results
            CL = sol["cl"]
            err = sol["cl"] - CLStar
            t2 = time.time()
            resultsDict = {
                "converged": converged,
                "iterations": _iIter + 1,
                "l2convergence": L2Conv,
                "alpha": aeroProblem.alpha,
                "cl": CL,
                "clstar": CLStar,
                "error": err,
                "clalpha": clalpha,
                "time": t2 - t1,
                "history": convergenceHistory,
            }

            if converged or (aeroProblem.solveFailed and stopOnStall) or aeroProblem.fatalFail:
                finalizeSolver(resultsDict, modifiedOptions)
                return resultsDict

        # we ran out of iterations and did not exit yet. so this is a failed case
        finalizeSolver(resultsDict, modifiedOptions)
        return resultsDict

    def solveTrimCL(
        self,
        aeroProblem,
        trimFunc,
        trimDV,
        dvIndex,
        CLStar,
        trimStar=0.0,
        alpha0=None,
        trim0=None,
        da=1e-3,
        deta=1e-2,
        tol=1e-4,
        nIter=10,
        Jac0=None,
        liftFunc="cl",
    ):
        """Solve the trim-Cl problem using a Broyden method.

        Parameters
        ----------
        ASProblem : ASProblem instance
            The aerostructural problem to be solved
        trimFunc : str
            Solution variable to use for trim. Usually 'cmy' or 'cmz'
        trimDV : str
            Dame of DVGeo design variable to control trim
        dvIndex : int
            Index of the trimDV function to use for trim
        CLStar : float
             Desired CL value
        trimStar : float
            Desired trimFunc value
        alpha0 : float or None
            Starting alpha. If None, use what is in the aeroProblem
        trim0 : float or None
            Starting trim value. If None, use what is in the DVGeo object
        da : float
            Initial alpha step for Jacobian
        deta : float
            Initial  stet in the 'eta' or trim dv function
        tol : float
            Tolerance for trimCL solve solution
        nIter : int
            Maximum number of iterations.
        Jac0 : 2x2 numpy array
            Initial guess for the trim-cl Jacobian. Usually obtained
            from a previous analysis and saves two function
            evaluations to produce the initial Jacobian.
        liftFunc : str
            Solution variable to use for lift. Usually 'cl' or a
            custom function created from cl.
        """

        self.setAeroProblem(aeroProblem)

        # Set defaults if they are not given.
        if alpha0 is None:
            alpha0 = self.curAP.alpha
        if trim0 is None:
            trim0 = 0.0

        def Func(Xn, CLstar):
            self.curAP.alpha = Xn[0]
            xGeo = self.DVGeo.getValues()
            xGeo[trimDV][dvIndex] = Xn[1]
            self.DVGeo.setDesignVars(xGeo)
            self.__call__(self.curAP, writeSolution=False)
            self.curAP.adflowData.callCounter -= 1
            funcs = {}
            self.evalFunctions(self.curAP, funcs, evalFuncs=[liftFunc, trimFunc])
            F = numpy.array(
                [funcs["%s_%s" % (self.curAP.name, liftFunc)] - CLstar, funcs["%s_%s" % (self.curAP.name, trimFunc)]]
            )
            return F

        # Generate initial point
        Xn = numpy.array([alpha0, trim0])
        Fn = Func(Xn, CLStar)

        # Next we generate the initial jacobian if we haven't been
        # provided with one.
        if Jac0 is not None:
            J = Jac0.copy()
        else:
            J = numpy.zeros((2, 2))

            # Perturb alpha
            Xn[0] += da
            Fpda = Func(Xn, CLStar)
            J[:, 0] = (Fpda - Fn) / da
            Xn[0] -= da

            # Perturb eta:
            Xn[1] += deta
            Fpdeta = Func(Xn, CLStar)
            J[:, 1] = (Fpdeta - Fn) / deta
            Xn[1] -= deta

        # Main iteration loop
        for jj in range(nIter):
            if self.comm.rank == 0:
                print("Fn:", Fn)

            # We now have a point and a jacobian...Newton's Method!
            Xnp1 = Xn - numpy.linalg.solve(J, Fn)
            if self.comm.rank == 0:
                print("Xnp1:", Xnp1)

            # Solve the new Xnp1
            Fnp1 = Func(Xnp1, CLStar)
            if self.comm.rank == 0:
                print("Fnp1:", Fnp1)

            # Update the jacobian using Broyden's method
            dx = Xnp1 - Xn
            dF = Fnp1 - Fn
            J = J + numpy.outer((dF - numpy.dot(J, dx)) / (numpy.linalg.norm(dx) ** 2), dx)

            if self.comm.rank == 0:
                print("New J:", J)

            # Shuffle the Fn and Xn backwards
            Fn = Fnp1.copy()
            Xn = Xnp1.copy()

            # Check for convergence
            if numpy.linalg.norm(Fn) < tol:
                if self.comm.rank == 0:
                    print("Converged!", jj, Fn, Xn)
                break
        return Xn

    def solveTargetFuncs(self, aeroProblem, funcDict, tol=1e-4, nIter=10, Jac0=None):
        """
        Solve the an arbitrary set of function-dv sets using a Broyden method.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aerodynamic problem to be solved
        funcDict : dict
            Dictionary of function DV pairs to solve:
            {'func':{'dv':str, 'dvIdx':idx,'target':val,'initVal':val,'initStep':val}}
            func : Name of function that is being solved for
            dv : design variable that has dominant control of this function value
            dvIdx : index into dv array if dv is not scalar
            target : target function value
            initVal : initial design variable value
            initStep : initial step for this dv in when generating finite difference starting jacobian
        tol : float
            Tolerance for the L2 norm of function error from target values
        nIter : int
            Maximum number of iterations.
        Jac0 : nxn numpy array
            Initial guess for the func-dv Jacobian. Usually obtained
            from a previous analysis and saves n function
            evaluations to produce the initial Jacobian.
        """

        self.setAeroProblem(aeroProblem)
        apDVList = aeroProblem.allVarFuncs
        # get a sorted list of keys so that we always access the dict in the same order
        funcKeys = sorted(funcDict.keys())

        # now parse everything out of the dictionary into lists
        dvKeys = []
        dvIndices = []
        targetVals = []
        initVals = []
        initStep = []
        for key in funcKeys:
            if "dv" in funcDict[key]:
                dvKey = funcDict[key]["dv"]
                if dvKey in apDVList:
                    dvKey += "_%s" % aeroProblem.name
                dvKeys.append(dvKey)
            else:
                raise Error("dv key not defined for function %s." % key)
                # Todo put in a check to see if the DVs selected are present?

            if "dvIdx" in funcDict[key]:
                dvIndices.append(funcDict[key]["dvIdx"])
            else:
                dvIndices.append(None)

            if "target" in funcDict[key]:
                targetVals.append(funcDict[key]["target"])
            else:
                targetVals.append(0.0)

            if "initVal" in funcDict[key]:
                initVals.append(funcDict[key]["initVal"])
            else:
                initVals.append(0.0)

            if "initStep" in funcDict[key]:
                initStep.append(funcDict[key]["initStep"])
            else:
                initStep.append(1e-1)

        nVars = len(funcKeys)

        def Func(Xn, targetVals):
            x = self.DVGeo.getValues()
            for i in range(nVars):
                if dvKeys[i] in x:
                    dvIdx = dvIndices[i]
                    if dvIdx is not None:
                        x[dvKeys[i]][dvIdx] = Xn[i]
                    else:
                        x[dvKeys[i]] = Xn[i]
                else:
                    x[dvKeys[i]] = Xn[i]

            self.curAP.setDesignVars(x)
            self.DVGeo.setDesignVars(x)

            self.__call__(self.curAP, writeSolution=False)
            self.curAP.adflowData.callCounter -= 1
            funcs = {}
            self.evalFunctions(self.curAP, funcs, evalFuncs=funcKeys)
            F = numpy.zeros([nVars])
            for i in range(nVars):
                F[i] = funcs["%s_%s" % (self.curAP.name, funcKeys[i])] - targetVals[i]

            return F

        # Generate initial point
        Xn = numpy.array(initVals)
        Fn = Func(Xn, targetVals)

        # Next we generate the initial jacobian if we haven't been
        # provided with one.
        if Jac0 is not None:
            J = Jac0.copy()
        else:
            J = numpy.zeros((nVars, nVars))

            # Perturb each var in succession
            for i in range(nVars):
                Xn[i] += initStep[i]
                Fpdelta = Func(Xn, targetVals)
                J[:, i] = (Fpdelta - Fn) / initStep[i]
                Xn[i] -= initStep[i]

        # Main iteration loop
        for jj in range(nIter):
            if self.comm.rank == 0:
                print("Fn:", Fn)

            # We now have a point and a jacobian...Newton's Method!
            Xnp1 = Xn - numpy.linalg.solve(J, Fn)
            if self.comm.rank == 0:
                print("Xnp1:", Xnp1)

            # Solve the new Xnp1
            Fnp1 = Func(Xnp1, targetVals)
            if self.comm.rank == 0:
                print("Fnp1:", Fnp1)

            # Update the jacobian using Broyden's method
            dx = Xnp1 - Xn
            dF = Fnp1 - Fn
            J = J + numpy.outer((dF - numpy.dot(J, dx)) / (numpy.linalg.norm(dx) ** 2), dx)

            if self.comm.rank == 0:
                print("New J:", J)

            # Shuffle the Fn and Xn backwards
            Fn = Fnp1.copy()
            Xn = Xnp1.copy()

            # Check for convergence
            if numpy.linalg.norm(Fn) < tol:
                if self.comm.rank == 0:
                    print("Converged!", "Iterations: ", jj, "Function Error:", Fn, "Variables: ", Xn)
                break
        return Xn

    def solveSep(
        self, aeroProblem, sepStar, nIter=10, alpha0=None, delta=0.1, tol=1e-3, expansionRatio=1.2, sepName=None
    ):
        """This is a safe-guarded secant search method to determine
        the alpha that yields a specified value of the separation
        sensor. Since this function is highly nonlinear we use a
        linear search to get the bounding range first.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aerodynamic problem to solve
        sepStar : float
            The desired target separation sensor value
        nIter : int
            Maximum number of iterations
        alpha0 : angle (deg) or None
            Initial guess. If none, use what is in aeroProblem.
        delta : angle (deg)
            Initial step. Only the magnitude is significant
        tol : float
            Desired tolerance for sepSensor

        sepName : str or None
            User supplied function to use for sep sensor. May be a
            user-added group function.

        Returns
        -------
        None, but the correct alpha is stored in the aeroProblem
        """
        self.setAeroProblem(aeroProblem)
        ap = aeroProblem
        # There are two main parts of the algorithm: First we need to
        # find the alpha range bracketing the desired function
        # value. Then we use a safe-guarded secant search to zero into
        # that vlaue.

        if sepName is None:
            sepName = "sepsensor"

        if alpha0 is None:
            alpha0 = ap.alpha

        # Name of function to use
        funcName = "%s_%s" % (ap.name, sepName)

        if not self.getOption("rkreset") and self.getOption("usenksolver"):
            ADFLOWWarning(
                "RKReset option is not set. It is usually necessary " "for solveSep() when NK solver is used."
            )

        # Solve first problem
        ap.alpha = alpha0
        self.__call__(ap, writeSolution=False)
        funcs = {}
        self.evalFunctions(ap, funcs, [sepName])
        f = funcs[funcName] - sepStar
        da = delta

        if self.comm.rank == 0:
            print("+----------------------------------------------+")
            print("| Presolve ")
            print("| Alpha      = %s" % ap.alpha)
            print("| Sep value  = %s" % funcs[funcName])
            print("| F          = %s" % f)
            print("| Sep*       = %s" % sepStar)
            print("+----------------------------------------------+")

        if self.comm.rank == 0:
            print("Searching for Correct Interval")

        # Now try to find the interval
        for i in range(nIter):
            if f > 0.0:
                da = -numpy.abs(da)
            else:
                da = numpy.abs(da)
            # Increment the alpha and solve again
            ap.alpha += da
            self.__call__(ap, writeSolution=False)
            self.evalFunctions(ap, funcs, [sepName])
            fnew = funcs[funcName] - sepStar

            if self.comm.rank == 0:
                print("+----------------------------------------------+")
                print("| Finished It= %d" % i)
                print("| Alpha      = %s" % ap.alpha)
                print("| Sep value  = %s" % funcs[funcName])
                print("| F          = %s" % fnew)
                print("| Sep*       = %s" % sepStar)
                print("+----------------------------------------------+")

            if numpy.sign(f) != numpy.sign(fnew):
                # We crossed the zero:
                break
            else:
                f = fnew
                # Expand the delta:
                da *= expansionRatio

        # Now we know anm2 and anm1 bracket the zero. We now start the
        # secant search but make sure that we stay inside of these known bounds
        anm1 = ap.alpha
        anm2 = ap.alpha - da
        fnm1 = fnew
        fnm2 = f
        lowAlpha = min(anm1, anm2)
        highAlpha = max(anm1, anm2)

        if self.comm.rank == 0:
            print("Switching to Secant Search")

        for iIter in range(nIter):
            if iIter != 0:
                ap.alpha = anm1
                self.__call__(ap, writeSolution=False)
                self.evalFunctions(ap, funcs, [sepName])
                fnm1 = funcs[funcName] - sepStar

            if self.comm.rank == 0:
                print("+----------------------------------------------+")
                print("| Finished It= %d" % i)
                print("| Alpha      = %s" % ap.alpha)
                print("| Sep value  = %s" % funcs[funcName])
                print("| F          = %s" % fnm1)
                print("| Sep*       = %s" % sepStar)
                print("+----------------------------------------------+")

            # Secant update
            anew = anm1 - fnm1 * (anm1 - anm2) / (fnm1 - fnm2)

            # Make sure the anew update doesn't go outside (lowAlpha,
            # highAlpha). Just use bisection in this case.
            if anew < lowAlpha:
                anew = 0.5 * (anm1 + lowAlpha)
            if anew > highAlpha:
                anew = 0.5 * (anm1 + highAlpha)

            # Shift the n-1 values to n-2
            fnm2 = fnm1
            anm2 = anm1

            # Set the new n-1 alpha value from update
            anm1 = anew

            # Finally, convergence check
            if abs(fnm1) < tol:
                break

    def writeSolution(self, outputDir=None, baseName=None, number=None, writeSlices=True, writeLift=True):
        """This is a generic shell function that potentially writes
        the various output files. The intent is that the user or
        calling program can call this file and ADflow write all the
        files that the user has defined. It is recommended that this
        function is used along with the associated logical flags in
        the options to determine the desired writing procedure

        Parameters
        ----------

        outputDir : str
            Use the supplied output directory
        baseName : str
            Use this supplied string for the base filename. Typically only used from an external solver.
        number : int
            Use the user supplied number to index solution. Again, only typically used from an external solver.
        writeSlices : bool
            Flag to determine if the slice files are written, if we have any slices added.
        writeLift : bool
            Flag to determine if the lift files are written, if we have any lift distributions added.
        """
        if outputDir is None:
            outputDir = self.getOption("outputDirectory")

        if baseName is None:
            baseName = self.curAP.name

        # If we are numbering solution, it saving the sequence of
        # calls, add the call number
        numDigits = self.getOption("writeSolutionDigits")
        if number is not None:
            # We need number based on the provided number:
            baseName = f"{baseName}_{number:0{numDigits}d}"
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
            if self.getOption("numberSolutions"):
                baseName = f"{baseName}_{self.curAP.adflowData.callCounter:0{numDigits}d}"

        # Join to get the actual filename root
        base = os.path.join(outputDir, baseName)

        # Get the timestep if this is time dependent solution, else we default
        # to the same values as the nsave values (which is 1) for steady
        # and time-spectral solution to be written,
        if self.getOption("equationMode").lower() == "unsteady":
            ts = self.adflow.monitor.timestepunsteady
        else:
            ts = 1

        # Now call each of the 4 routines with the appropriate file name:
        if self.getOption("writevolumesolution") and numpy.mod(ts, self.getOption("nsavevolume")) == 0:
            self.writeVolumeSolutionFile(base + "_vol.cgns")

        if self.getOption("writesurfacesolution") and numpy.mod(ts, self.getOption("nsavesurface")) == 0:
            self.writeSurfaceSolutionFile(base + "_surf.cgns")

        # ADD NSAVE TEST AROUND THESE AS WELL BUT
        # THIS IS SMALL COMPARED TO OTHER. REFACTOR
        if self.getOption("equationMode").lower() == "unsteady":
            liftName = base + "_lift_Timestep%4.4d.dat" % ts
            sliceName = base + "_slices_Timestep%4.4d.dat" % ts
            surfName = base + "_surf_Timestep%4.4d.plt" % ts
        else:
            liftName = base + "_lift.dat"
            sliceName = base + "_slices.dat"
            surfName = base + "_surf.plt"

        # Get the family list to write for the surface.
        famList = self._getFamilyList(self.getOption("outputSurfaceFamily"))

        # Flag to write the tecplot surface solution or not
        writeSurf = self.getOption("writeTecplotSurfaceSolution")

        # # Call fully compbined fortran routine.
        self.adflow.tecplotio.writetecplot(sliceName, writeSlices, liftName, writeLift, surfName, writeSurf, famList)

    def writeMeshFile(self, fileName):
        """Write the current mesh to a CGNS file. This call isn't used
        normally since the volume solution usually contains the grid

        Parameters
        ----------
        fileName : str
            Name of the mesh file

        """

        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".cgns"

        # Set Flags for writing
        self.adflow.monitor.writegrid = True
        self.adflow.monitor.writevolume = False
        self.adflow.monitor.writesurface = False

        # Set filename in adflow
        self.adflow.inputio.solfile = self._expandString(fileName)
        self.adflow.inputio.newgridfile = self._expandString(fileName)

        # Actual fortran write call. Family list doesn't matter.
        famList = self._getFamilyList(self.allFamilies)
        self.adflow.writesol(famList)

    def writeVolumeSolutionFile(self, fileName, writeGrid=True):
        """Write the current state of the volume flow solution to a CGNS
        file. This is a lower level routine; Normally one should call
        ``writeSolution()``.

        Parameters
        ----------
        fileName : str
            Name of the file. Should have .cgns extension.
        writeGrid : bool
            Flag specifying whether the grid should be included or if
            links should be used. Always writing the grid is
            recommended even in cases when it is not strictly necessary.
            Note that if ``writeGrid = False`` the volume files do not contain
            any grid coordinates rendering the file useless if a separate
            grid file was written out and is linked to it.
        """
        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".cgns"

        # Set Flags for writing
        self.adflow.monitor.writegrid = writeGrid
        self.adflow.monitor.writevolume = True
        self.adflow.monitor.writesurface = False

        # Set fileName in adflow
        self.adflow.inputio.solfile = self._expandString(fileName)
        # self.adflow.inputio.solfile[0:len(fileName)] = fileName
        self.adflow.inputio.newgridfile = self._expandString(fileName)
        # self.adflow.inputio.newgridfile[0:len(fileName)] = fileName

        # Actual fortran write call. family list doesn't matter
        famList = self._getFamilyList(self.allFamilies)
        self.adflow.writesol(famList)

    def writeSurfaceSolutionFile(self, fileName):
        """Write the current state of the surface flow solution to a CGNS file.
        Keyword arguments:

        Parameters
        ----------
        fileName : str
            Name of the file. Should have .cgns extension.
        """
        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".cgns"

        # Set Flags for writing
        self.adflow.monitor.writegrid = False
        self.adflow.monitor.writevolume = False
        self.adflow.monitor.writesurface = True

        # Set fileName in adflow
        self.adflow.inputio.surfacesolfile = self._expandString(fileName)

        # Actual fortran write call. Fam list matters.
        famList = self._getFamilyList(self.getOption("outputSurfaceFamily"))
        self.adflow.writesol(famList)

    def writeSurfaceSolutionFileTecplot(self, fileName):
        """Write the current state of the surface flow solution to a teclot file.

        Parameters
        ----------
        fileName : str
            Name of the file. Should have .plt extension.
        """
        # Ensure fileName is .dat even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".plt"

        # Actual write command
        sliceName = ""
        liftName = ""
        famList = self._getFamilyList(self.getOption("outputSurfaceFamily"))
        self.adflow.tecplotio.writetecplot(sliceName, False, liftName, False, fileName, True, famList)

    def writeLiftDistributionFile(self, fileName):
        """Evaluate and write the lift distibution to a tecplot file.

        Parameters
        ----------
        fileName : str
            File of lift distribution. Should have .dat extension.
        """
        # Ensure fileName is .dat even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".dat"

        # Actual write command. Family list doesn't matter.
        sliceName = ""
        surfName = ""
        self.adflow.tecplotio.writetecplot(sliceName, False, fileName, True, surfName, False, self.allFamilies)

    def writeSlicesFile(self, fileName):
        """Evaluate and write the defined slice information to a
        tecplot file.

        Parameters
        ----------
        fileName : str
            Slice file. Should have .dat extension.
        """

        # Ensure filename is .dat even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += ".dat"

        # Actual write command. Family list doesn't matter
        liftName = ""
        surfName = ""
        self.adflow.tecplotio.writetecplot(fileName, True, liftName, True, surfName, False, self.allFamilies)

    def writeForceFile(self, fileName, TS=0, groupName=None, cfdForcePts=None):
        """This function collects all the forces and locations and
        writes them to a file with each line having: X Y Z Fx Fy Fz.
        This can then be used to set a set of structural loads in TACS
        for structural only optimization

        Like the getForces() routine, an external set of forces may be
        passed in on which to evaluate the forces. This is only
        typically used in an aerostructural case.

        """
        if groupName is None:
            groupName = self.allWallsGroup

        # Now we need to gather the data:
        if cfdForcePts is None:
            pts = self.comm.gather(self.getSurfacePoints(groupName, True, TS), root=0)
        else:
            pts = self.comm.gather(cfdForcePts)

        # Get the forces and connectivity
        forces = self.comm.gather(self.getForces(groupName, TS=TS), root=0)
        conn, faceSize = self.getSurfaceConnectivity(groupName, True)
        conn = self.comm.gather(conn, root=0)

        # Write out Data only on root proc:
        if self.myid == 0:
            # First sum up the total number of nodes and elements:
            nPt = 0
            nCell = 0
            for iProc in range(len(pts)):
                nPt += len(pts[iProc])
                nCell += len(conn[iProc])

            # Open output file
            f = open(fileName, "w")

            # Write header with number of nodes and number of cells
            f.write("%d %d\n" % (nPt, nCell))

            # Now write out all the Nodes and Forces (or tractions)
            for iProc in range(len(pts)):
                for i in range(len(pts[iProc])):
                    f.write(
                        "%15.8g %15.8g %15.8g "
                        % (numpy.real(pts[iProc][i, 0]), numpy.real(pts[iProc][i, 1]), numpy.real(pts[iProc][i, 2]))
                    )
                    f.write(
                        "%15.8g %15.8g %15.8g\n"
                        % (
                            numpy.real(forces[iProc][i, 0]),
                            numpy.real(forces[iProc][i, 1]),
                            numpy.real(forces[iProc][i, 2]),
                        )
                    )

            # Now write out the connectivity information. We have to
            # be a little careful, since the connectivitiy is given
            # locally per proc. As we loop over the data from each
            # proc, we need to increment the connectivity by the
            # number of nodes we've "used up" so far
            nodeOffset = 0
            for iProc in range(len(conn)):
                for i in range(len(conn[iProc])):
                    f.write(
                        "%d %d %d %d\n"
                        % (
                            conn[iProc][i, 0] + nodeOffset,
                            conn[iProc][i, 1] + nodeOffset,
                            conn[iProc][i, 2] + nodeOffset,
                            conn[iProc][i, 3] + nodeOffset,
                        )
                    )

                nodeOffset += len(pts[iProc])
            f.close()
        # end if (root proc )

    def writeSurfaceSensitivity(self, fileName, func, groupName=None):
        """Write a tecplot file of the surface sensitivity. It is up to the use
        to make sure the adjoint already computed before calling this
        function.

        Parameters
        ----------
        fileName : str
           String for output filename. Should end in .dat for tecplot.
        func : str
           ADFlow objective string for the objective to write.
        groupName : str
           Family group to use. Default to all walls if not given (None)
        """

        if groupName is None:
            groupName = self.allWallsGroup

        func = func.lower()
        # Make sure that we have the adjoint we need:
        if func in self.curAP.adflowData.adjoints:
            # These are the reverse mode seeds
            psi = -self.getAdjoint(func)

            # Compute everything and update into the dictionary
            dXs = self.computeJacobianVectorProductBwd(resBar=psi, funcsBar=self._getFuncsBar(func), xSDeriv=True)
            dXs = self.mapVector(dXs, groupName, groupName, includeZipper=True)
        else:
            raise Error(
                "The adjoint for '%s' is not computed for the current "
                "aeroProblem. cannot write surface sensitivity." % (func)
            )

        # Be careful with conn...need to increment by the offset from
        # the points.
        pts = self.getSurfacePoints(groupName, includeZipper=True)

        ptSizes = self.comm.allgather(len(pts))
        offsets = numpy.zeros(len(ptSizes), "intc")
        offsets[1:] = numpy.cumsum(ptSizes)[:-1]

        # Now we need to gather the data to the root proc
        pts = self.comm.gather(pts)
        dXs = self.comm.gather(dXs)
        conn, faceSize = self.getSurfaceConnectivity(groupName, includeZipper=True)
        conn = self.comm.gather(conn + offsets[self.comm.rank])

        # Write out Data only on root proc:
        if self.myid == 0:
            # Make numpy arrays of all data
            pts = numpy.vstack(pts)
            dXs = numpy.vstack(dXs)
            conn = numpy.vstack(conn)

            # Next point reduce all nodes:
            uniquePts, link, nUnique = self.adflow.utils.pointreduce(pts.T, 1e-12)

            # Insert existing nodes into the new array. *Accumulate* dXs.
            newdXs = numpy.zeros((nUnique, 3))
            newPts = numpy.zeros((nUnique, 3))
            newConn = numpy.zeros_like(conn)
            for i in range(len(pts)):
                newPts[link[i] - 1] = pts[i]
                newdXs[link[i] - 1] += dXs[i]

            # Update conn. Since link is 1 based, we get 1 based order
            # which is what we need for tecplot
            for i in range(len(conn)):
                for j in range(4):
                    newConn[i, j] = link[conn[i, j]]

            pts = newPts
            dXs = newdXs
            conn = newConn

            # Open output file
            f = open(fileName, "w")

            # Write the variable headers
            f.write("Variables = CoordinateX CoordinateY CoordinateZ dX dY dZ\n")
            f.write(
                "ZONE Nodes=%d Elements=%d Zonetype=FEQuadrilateral \
            Datapacking=Point\n"
                % (len(pts), len(conn))
            )

            # Now write out all the Nodes and Forces (or tractions)
            for i in range(len(pts)):
                f.write("%15.8g %15.8g %15.8g " % (pts[i, 0], pts[i, 1], pts[i, 2]))
                f.write("%15.8g %15.8g %15.8g\n" % (dXs[i, 0], dXs[i, 1], dXs[i, 2]))

            for i in range(len(conn)):
                f.write("%d %d %d %d\n" % (conn[i, 0], conn[i, 1], conn[i, 2], conn[i, 3]))

            f.close()
        # end if (root proc )

    def resetAdjoint(self, obj):
        """
        Reset an adjoint 'obj' that is stored in the current
        aeroProblem. If the adjoint does not yet exist, nothing is done

        Parameters
        ----------
        obj : str
            String identifing the objective.
        """
        if obj in self.curAP.adflowData.adjoints:
            self.curAP.adflowData.adjoints[obj][:] = 0.0

    def resetFlow(self, aeroProblem, releaseAdjointMemory=True):
        """
        Reset the flow after a failure or for a complex step
        derivative operation.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aeroproblem with the flow information we would like
            to reset the flow to.
        """

        # The aeroProblem we're resetting to has an oldWinf in it, we
        # must invalidate it since it would try to use that the next
        # time this AP is used.
        try:
            aeroProblem.adflowData.oldWinf = None
        except AttributeError:
            pass

        self.setAeroProblem(aeroProblem, releaseAdjointMemory)
        self._resetFlow()
        # also reset the ank cfl for this AP
        self.resetANKCFL(aeroProblem)

    def _resetFlow(self):
        strLvl = self.getOption("MGStartLevel")
        nLevels = self.adflow.inputiteration.nmglevels
        if strLvl < 0 or strLvl > nLevels:
            strLvl = nLevels

        self.adflow.inputiteration.mgstartlevel = strLvl
        self.adflow.iteration.groundlevel = strLvl
        self.adflow.iteration.currentlevel = strLvl
        self.adflow.monitor.nitercur = 0
        self.adflow.iteration.itertot = 0
        self.adflow.initializeflow.setuniformflow()
        self.adflow.killsignals.routinefailed = False
        self.adflow.killsignals.fatalfail = False
        self.adflow.nksolver.freestreamresset = False

    def resetANKCFL(self, aeroProblem):
        """
        Resets the ANK CFL of the aeroProblem to the value set by "ANKCFL0" option.
        During the first ANK iteration that follows, this will be overwritten
        by the ANK solver itself with the CFL Min logic based on relative
        convergence. This is useful if keeping a really high ANK CFL is not desired.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aeroproblem whose ANK CFL will be reset.
        """
        aeroProblem.adflowData.ank_cfl = self.getOption("ANKCFL0")

    def getSolution(self, groupName=None):
        """Retrieve the basic solution variables from the solver. This will
        return all variables defined in basicCostFunctions for the
        specified group. This is a lower level function than
        evalFunctions() which should be used for optimization.

        Parameters
        ----------
        groupName : str
            The family group on which to evaluate the functions.

        """

        # Check that the user surface integrations are finalized
        if self.hasIntegrationSurfaces and (not self.userSurfaceIntegrationsFinalized):
            raise Error(
                "Integration surfaces have been added to the solver, but finalizeUserIntegrationSurfaces()"
                + " has not been called. It must be called before evaluating the funcs"
            )

        # Extract the familiy list we want to use for evaluation. We
        # explictly have just the one group in the call.
        famLists = self._expandGroupNames([groupName])

        # Run the underlying fortran routine
        costSize = self.adflow.constants.ncostfunction
        funcVals = self.adflow.surfaceintegrations.getsolutionwrap(famLists, costSize)

        # Build up the dictionary
        sol = {}
        for key in self.basicCostFunctions:
            sol[key] = funcVals[self.basicCostFunctions[key] - 1, 0]

        return sol

    def getIblankCheckSum(self, fileName=None):
        ncells = self.adflow.adjointvars.ncellslocal[0]
        nCellTotal = self.comm.allreduce(ncells)
        if self.myid != 0:
            nCellTotal = 1  # Should be zero, but f2py doesn't like
            # that

        blkList = self.adflow.oversetutilities.getoversetiblank(nCellTotal * 5)
        hexStr = None
        if self.myid == 0:
            hexStr = hashlib.md5(str(list(blkList[0::5]))).hexdigest()
            if fileName is not None:
                f = open(fileName, "w")
                for i in range(nCellTotal):
                    f.write(
                        "%d %d %d %d %d\n"
                        % (
                            blkList[5 * i],
                            blkList[5 * i + 1],
                            blkList[5 * i + 2],
                            blkList[5 * i + 3],
                            blkList[5 * i + 4],
                        )
                    )
                f.close()

        return self.comm.bcast(hexStr)

    # =========================================================================
    #   The following routines are public functions however, they should
    #   not need to be used by a user using this class directly. They are
    #   primarly used internally OR by a solver that uses this class
    #   i.e. an Aerostructural solver
    # =========================================================================

    def getSurfaceCoordinates(self, groupName=None, includeZipper=True, TS=0):
        # This is an alias for getSurfacePoints
        return self.getSurfacePoints(groupName, includeZipper, TS)

    def getPointSetName(self, apName):
        """
        Take the apName and return the mangled point set name.

        """
        return "adflow_%s_coords" % apName

    def setSurfaceCoordinates(self, coordinates, groupName=None):
        """
        Set the updated surface coordinates for a particular group.

        Parameters
        ----------
        coordinates : numpy array
            Numpy array of size Nx3, where N is the number of coordinates on this processor.
            This array must have the same shape as the array obtained with getSurfaceCoordinates()

        groupName : str
            Name of family or group of families for which to return coordinates for.

        """
        if self.mesh is None:
            return

        if groupName is None:
            groupName = self.allWallsGroup

        self._updateGeomInfo = True
        if self.mesh is None:
            raise Error("Cannot set new surface coordinate locations without a mesh" "warping object present.")

        # First get the surface coordinates of the meshFamily in case
        # the groupName is a subset, those values will remain unchanged.
        meshSurfCoords = self.getSurfaceCoordinates(self.meshFamilyGroup, includeZipper=False)
        meshSurfCoords = self.mapVector(
            coordinates, groupName, self.meshFamilyGroup, meshSurfCoords, includeZipper=False
        )
        self.mesh.setSurfaceCoordinates(meshSurfCoords)

    def setAeroProblem(self, aeroProblem, releaseAdjointMemory=True):
        """Set the supplied aeroProblem to be used in ADflow.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The supplied aeroproblem to be set.
        releaseAdjointMemory : bool, optional
            Flag to release the adjoint memory when setting a new aeroproblem, by default True
        """

        ptSetName = "adflow_%s_coords" % aeroProblem.name

        newAP = False
        # Tell the user if we are switching aeroProblems
        if self.curAP != aeroProblem:
            newAP = True
            self.pp("+" + "-" * 70 + "+")
            self.pp("|  Switching to Aero Problem: %-41s|" % aeroProblem.name)
            self.pp("+" + "-" * 70 + "+")

        # See if the aeroProblem has adflowData already, if not, create.
        try:
            aeroProblem.adflowData
        except AttributeError:
            aeroProblem.adflowData = adflowFlowCase()
            aeroProblem.ptSetName = ptSetName
            aeroProblem.surfMesh = self.getSurfaceCoordinates(self.designFamilyGroup)

        if self.curAP is not None:
            # If we have already solved something and are now
            # switching, save what we need:
            self.curAP.stateInfo = self._getInfo()
            self.curAP.surfMesh = self.getSurfaceCoordinates(self.designFamilyGroup)

            # Restore any options that the current aeroProblem
            # (self.curAP) modified. We have to be slightly careful
            # since a setOption() may have been called in between:

            for key in self.curAP.savedOptions["adflow"]:
                # Saved Val: This is the main option value when the
                # aeroProblem set it's own option
                # setVal: This is the value the aeroProblem set itself
                # curVal: This is the actual value currently set

                savedVal = self.curAP.savedOptions["adflow"][key]
                setVal = self.curAP.solverOptions["adflow"][key]
                curVal = self.getOption(key)

                if curVal == setVal:
                    # Restore the saved value, if it still what the
                    # aeroProblem had set
                    self.setOption(key, savedVal)

        # Now check if we have an DVGeo object to deal with:
        if self.DVGeo is not None:
            # DVGeo appeared and we have not embedded points!
            if ptSetName not in self.DVGeo.points:
                coords0 = self.mapVector(self.coords0, self.allFamilies, self.designFamilyGroup, includeZipper=False)
                self.DVGeo.addPointSet(coords0, ptSetName, **self.pointSetKwargs)

            # Check if our point-set is up to date:
            if not self.DVGeo.pointSetUpToDate(ptSetName) or aeroProblem.adflowData.disp is not None:
                coords = self.DVGeo.update(ptSetName, config=aeroProblem.name)

                # Potentially add a fixed set of displacements to it.
                if aeroProblem.adflowData.disp is not None:
                    coords += self.curAP.adflowData.disp

                self.setSurfaceCoordinates(coords, self.designFamilyGroup)

        self._setAeroProblemData(aeroProblem)

        # Remind the user of the modified options when switching AP
        if newAP and self.getOption("printIterations"):
            self.printModifiedOptions()

        # Reset the fail flags here since the updateGeometry info
        # updates the mesh which may result in a fatalFail.
        self.adflow.killsignals.routinefailed = False
        self.adflow.killsignals.fatalFail = False

        # Note that it is safe to call updateGeometryInfo since it
        # only updates the mesh if necessary

        self.updateGeometryInfo()
        # Create the zipper mesh if not done so
        self._createZipperMesh()

        failedMesh = False
        if self.adflow.killsignals.fatalfail:
            failedMesh = True
        # If the state info none, initialize to the supplied
        # restart file or the free-stream values by calling
        # resetFlow()
        if aeroProblem.adflowData.stateInfo is None:
            if self.getOption("restartFile") is not None:
                self.adflow.inputiteration.mgstartlevel = 1
                self.adflow.initializeflow.initflowrestart()
                # Save this state information
                aeroProblem.adflowData.stateInfo = self._getInfo()
            else:
                self._resetFlow()

            aeroProblem.adflowData.ank_cfl = self.getOption("ANKCFL0")

        if failedMesh:
            self.adflow.killsignals.fatalfail = True

        # Now set the data from the incomming aeroProblem:
        stateInfo = aeroProblem.adflowData.stateInfo
        if stateInfo is not None and newAP:
            self._setInfo(stateInfo)

        # Potentially correct the states based on the change in the alpha
        oldWinf = aeroProblem.adflowData.oldWinf
        if self.getOption("infChangeCorrection") and oldWinf is not None:
            correctionTol = self.getOption("infChangeCorrectionTol")
            correctionType = self.getOption("infChangeCorrectionType")
            self.adflow.initializeflow.infchangecorrection(oldWinf, correctionTol, correctionType)

        # We are now ready to associate self.curAP with the supplied AP
        self.curAP = aeroProblem
        self.curAP.adjointRHS = None

        # Destroy the ANK/NK solver and the adjoint memory
        if newAP:
            self.adflow.nksolver.destroynksolver()
            self.adflow.anksolver.destroyanksolver()
            if releaseAdjointMemory:
                self.releaseAdjointMemory()

        if not self.getOption("ANKCFLReset"):
            # also set the ANK CFL to the last value if we are restarting
            self.adflow.anksolver.ank_cfl = aeroProblem.adflowData.ank_cfl

    def _setAeroProblemData(self, aeroProblem, firstCall=False):
        """After an aeroProblem has been associated with self.curAP, set
        all the updated information in ADflow.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The current aeroProblem object.
        firstCall : bool, optional
            Flag that signifies this is being called for the first time, by default False
        """

        # Set any additional adflow options that may be defined in the
        # aeroproblem. While we do it we save the options that we've
        # modified if they are different than the current option.
        AP = aeroProblem
        try:
            AP.savedOptions
        except AttributeError:
            AP.savedOptions = {"adflow": {}}

        if "adflow" in AP.solverOptions:
            for key in AP.solverOptions["adflow"]:
                curVal = self.getOption(key)
                overwriteVal = AP.solverOptions["adflow"][key]
                if overwriteVal != curVal:
                    self.setOption(key, overwriteVal)
                    AP.savedOptions["adflow"][key] = curVal

        alpha = AP.alpha
        beta = AP.beta
        mach = AP.mach
        machRef = AP.machRef
        machGrid = AP.machGrid

        xRef = AP.xRef
        yRef = AP.yRef
        zRef = AP.zRef
        xRot = AP.xRot
        yRot = AP.yRot
        zRot = AP.zRot
        momentAxis = AP.momentAxis
        areaRef = AP.areaRef
        chordRef = AP.chordRef
        liftIndex = self.getOption("liftIndex")

        if AP.T is None or AP.P is None or AP.rho is None or AP.V is None or AP.mu is None:
            raise Error(
                "Insufficient information is given in the "
                "aeroProblem to determine physical state. "
                "See AeroProblem documentation for how to "
                "specify complete aerodynamic states."
            )

        if self.dtype == "d":
            mach = numpy.real(mach)

        if self.dtype == "d":
            T = numpy.real(AP.T)
            P = numpy.real(AP.P)
            rho = numpy.real(AP.rho)

            SSuthDim = numpy.real(AP.SSuthDim)
            muSuthDim = numpy.real(AP.muSuthDim)
            TSuthDim = numpy.real(AP.TSuthDim)
            RGasDim = numpy.real(AP.R)
            gammaConstant = numpy.real(AP.gamma)
            Pr = numpy.real(AP.Pr)
        else:
            T = AP.T
            P = AP.P
            rho = AP.rho

            SSuthDim = AP.SSuthDim
            muSuthDim = AP.muSuthDim
            TSuthDim = AP.TSuthDim
            RGasDim = AP.R
            gammaConstant = AP.gamma
            Pr = AP.Pr

        # Do some checking here for things that MUST be specified:
        if AP.mach is None:
            raise Error("'mach' number must be specified in the aeroProblem" " for ADflow.")
        if areaRef is None:
            raise Error("'areaRef' must be specified in aeroProblem" " for ADflow.")
        if chordRef is None:
            raise Error("'chordRef' must be specified in aeroProblem" " for ADflow.")

        # Now set defaults
        if alpha is None:
            alpha = 0.0
        if beta is None:
            beta = 0.0
        if xRef is None:
            xRef = 0.0
        if yRef is None:
            yRef = 0.0
        if zRef is None:
            zRef = 0.0
        if xRot is None:
            xRot = 0.0
        if yRot is None:
            yRot = 0.0
        if zRot is None:
            zRot = 0.0

        if momentAxis is None:  # Set the default to the x-axis through the origin
            axisX1 = 0.0
            axisX2 = 1.0
            axisY1 = 0.0
            axisY2 = 0.0
            axisZ1 = 0.0
            axisZ2 = 0.0
        else:
            axisX1 = momentAxis[0][0]
            axisX2 = momentAxis[1][0]
            axisY1 = momentAxis[0][1]
            axisY2 = momentAxis[1][1]
            axisZ1 = momentAxis[0][2]
            axisZ2 = momentAxis[1][2]

        # Set mach defaults if user did not specified any machRef or machGrid values

        # If the user is running time spectral but did not specify
        # machGrid then default it to be equal to mach.

        if machRef is None:
            machRef = mach

        if machGrid is None:
            if self.getOption("equationMode").lower() == "time spectral":
                machGrid = mach
            else:
                # Steady, unsteady
                machGrid = 0.0

        # 1. Angle of attack:
        dToR = numpy.pi / 180.0
        self.adflow.inputphysics.alpha = alpha * dToR
        self.adflow.inputphysics.beta = beta * dToR
        self.adflow.inputphysics.liftindex = liftIndex
        self.adflow.flowutils.adjustinflowangle()

        if self.getOption("printIterations") and self.comm.rank == 0:
            print("-> Alpha... %f " % numpy.real(alpha))

        # 2. Reference Points:
        self.adflow.inputphysics.pointref = [xRef, yRef, zRef]
        self.adflow.inputphysics.momentaxis = [[axisX1, axisX2], [axisY1, axisY2], [axisZ1, axisZ2]]
        self.adflow.inputmotion.rotpoint = [xRot, yRot, zRot]
        self.adflow.inputphysics.pointrefec = [0.0, 0.0, 0.0]

        # 3. Reference Areas
        self.adflow.inputphysics.surfaceref = areaRef
        self.adflow.inputphysics.lengthref = chordRef

        # 4. Set mach numbers
        # Mach number for time spectral needs to be set to set to zero
        # If time-spectral (TS) then mach = 0, machcoef = mach, machgrid = mach
        # If Steady-State (SS), time-accurate (TA) then mach = mach, machcoef = mach, machgrid = 0
        if self.getOption("equationMode").lower() == "time spectral":
            self.adflow.inputphysics.mach = 0.0
        else:
            self.adflow.inputphysics.mach = mach

        self.adflow.inputphysics.machcoef = machRef
        self.adflow.inputphysics.machgrid = machGrid

        # Set reference state information:

        # Set the dimensional free strema values
        self.adflow.flowvarrefstate.pinfdim = P
        self.adflow.flowvarrefstate.tinfdim = T
        self.adflow.flowvarrefstate.rhoinfdim = rho

        self.adflow.inputphysics.ssuthdim = SSuthDim
        self.adflow.inputphysics.musuthdim = muSuthDim
        self.adflow.inputphysics.tsuthdim = TSuthDim
        self.adflow.inputphysics.rgasdim = RGasDim
        self.adflow.inputphysics.prandtl = Pr

        # Update gamma only if it has changed from what currently is set
        if abs(self.adflow.inputphysics.gammaconstant - gammaConstant) > 1.0e-12:
            self.adflow.inputphysics.gammaconstant = gammaConstant
            self.adflow.flowutils.updategamma()  # NOTE! It is absolutely necessary to call this function, otherwise gamma is not properly updated.

        # 4. Periodic Parameters --- These are not checked/verified
        # and come directly from aeroProblem. Make sure you specify
        # them there properly!!
        if self.getOption("useExternalDynamicMesh"):
            # The inputs here are mainly used for TS stability problem.
            # For general aeroelastic problem, we rely on more general settings.
            # Currently, polynomial motion is not supported for aeroelastic problem.
            AP.degreePol = 0
            AP.coefPol = [0.0]
            #  no odd number of instance allowed!
            AP.degreeFourier = (self.adflow.inputtimespectral.ntimeintervalsspectral - 1) // 2
            AP.cosCoefFourier = [0.0]
            AP.sinCoefFourier = [0.0]
        # if use time spectral, we need one of the following
        # to be true to get a correct time period.
        # the number of time instance is set directly.
        if self.getOption("alphaMode"):
            self.adflow.inputmotion.degreepolalpha = int(AP.degreePol)
            self.adflow.inputmotion.coefpolalpha = AP.coefPol
            self.adflow.inputmotion.omegafouralpha = AP.omegaFourier
            self.adflow.inputmotion.degreefouralpha = AP.degreeFourier
            self.adflow.inputmotion.coscoeffouralpha = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffouralpha = AP.sinCoefFourier
        elif self.getOption("betaMode"):
            self.adflow.inputmotion.degreepolmach = int(AP.degreePol)
            self.adflow.inputmotion.coefpolmach = AP.coefPol
            self.adflow.inputmotion.omegafourbeta = AP.omegaFourier
            self.adflow.inputmotion.degreefourbeta = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourbeta = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourbeta = AP.sinCoefFourier
        elif self.getOption("machMode"):
            self.adflow.inputmotion.degreepolmach = int(AP.degreePol)
            self.adflow.inputmotion.coefpolmach = AP.coefPol
            self.adflow.inputmotion.omegafourmach = AP.omegaFourier
            self.adflow.inputmotion.degreefourmach = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourmach = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourmach = AP.sinCoefFourier
        elif self.getOption("pMode"):
            ### add in lift axis dependence
            self.adflow.inputmotion.degreepolxrot = int(AP.degreePol)
            self.adflow.inputmotion.coefpolxrot = AP.coefPol
            self.adflow.inputmotion.omegafourxrot = AP.omegaFourier
            self.adflow.inputmotion.degreefourxrot = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourxrot = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourxrot = AP.sinCoefFourier
        elif self.getOption("qMode"):
            self.adflow.inputmotion.degreepolzrot = int(AP.degreePol)
            self.adflow.inputmotion.coefpolzrot = AP.coefPol
            self.adflow.inputmotion.omegafourzrot = AP.omegaFourier
            self.adflow.inputmotion.degreefourzrot = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourzrot = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourzrot = AP.sinCoefFourier
        elif self.getOption("rMode"):
            self.adflow.inputmotion.degreepolyrot = int(AP.degreePol)
            self.adflow.inputmotion.coefpolyrot = AP.coefPol
            self.adflow.inputmotion.omegafouryrot = AP.omegaFourier
            self.adflow.inputmotion.degreefouryrot = AP.degreeFourier
            self.adflow.inputmotion.coscoeffouryrot = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffouryrot = AP.sinCoefFourier
        elif self.getOption("useExternalDynamicMesh"):
            # if it is an aeroelastic case
            self.adflow.inputtimespectral.omegafourier = AP.omegaFourier

        # Set any possible BC Data coming out of the aeroProblem
        nameArray, dataArray, groupArray, groupNames, empty = self._getBCDataFromAeroProblem(AP)
        if not empty:
            self.adflow.bcdata.setbcdata(nameArray, dataArray, groupArray, 1)

        if not firstCall:
            self.adflow.initializeflow.updatebcdataalllevels()

            # We need to do an all reduce here in case there's a BC failure
            # on any of the procs after we update the bc data.
            self.adflow.killsignals.fatalfail = self.comm.allreduce(bool(self.adflow.killsignals.fatalfail), op=MPI.LOR)

            if self.getOption("equationMode").lower() == "time spectral":
                self.adflow.preprocessingapi.updateperiodicinfoalllevels()
                self.adflow.preprocessingapi.updatemetricsalllevels()
            self.adflow.preprocessingapi.updategridvelocitiesalllevels()

    def _getBCDataFromAeroProblem(self, AP):
        variables = []
        dataArray = []
        groupNames = []

        for tmp in AP.bcVarData:
            varName, family = tmp
            variables.append(varName)
            dataArray.append(AP.bcVarData[tmp])
            groupNames.append(family)

        nameArray = self._createFortranStringArray(variables)
        groupArray = self._expandGroupNames(groupNames)
        if len(nameArray) > 0:
            return nameArray, dataArray, groupArray, groupNames, False
        else:
            # dummy data that doesn't matter
            return (self._createFortranStringArray(["Pressure"]), [0.0], [[1, 1]], groupNames, True)

    def getForces(self, groupName=None, TS=0):
        """Return the forces on this processor on the families defined by groupName.

        Parameters
        ----------
        groupName : str
            Group identifier to get only forces cooresponding to the
            desired group. The group must be a family or a user-supplied
            group of families. The default is None which corresponds to
            all wall-type surfaces.

        TS : int
            Spectral instance for which to get the forces

        Returns
        -------
        forces : array (N,3)
            Forces (or tractions depending on that forceAsTractions
            options) on this processor. Note that N may be 0, and an
            empty array of shape (0, 3) can be returned.
        """
        # Set the family to all walls group.
        npts, ncell = self._getSurfaceSize(self.allWallsGroup)

        forces = numpy.zeros((npts, 3), self.dtype)
        self.adflow.getforces(forces.T, TS + 1)
        if groupName is None:
            groupName = self.allWallsGroup

        # Finally map the vector as required.
        return self.mapVector(forces, self.allWallsGroup, groupName)

    def getHeatFluxes(self, groupName=None, TS=0):
        """Return the heat fluxes for isothermal walls on the families
        defined by group name on this processor.

        Parameters
        ----------
        groupName : str
            Group identifier to get only heat fluxes cooresponding to
            the desired group. The group must be a family or a
            user-supplied group of families. The default is None which
            corresponds to all wall-type surfaces.

        TS : int
            Spectral instance for which to get the fluxes.

        Returns
        -------
        heatFluxes : array (N)
            HeatFluxes on this processor. Note that N may be 0, and an
            empty array of shape (0) can be returned.

        """

        # Set the family to all walls group.
        npts, ncell = self._getSurfaceSize(self.allWallsGroup)

        fluxes = numpy.zeros(npts, self.dtype)
        self.adflow.getheatflux(fluxes, TS + 1)
        if groupName is None:
            groupName = self.allWallsGroup

        # Map vector expects and Nx3 array. So we will do just that.
        tmp = numpy.zeros((npts, 3))
        tmp[:, 0] = fluxes
        tmp = self.mapVector(tmp, self.allWallsGroup, groupName)
        fluxes = tmp[:, 0]
        return fluxes

    def setWallTemperature(self, temperature, groupName=None, TS=0):
        """Set the temperature of the isothermal walls.

        Parameters
        ----------
        temperature : numpy array

            Dimensional temperature to set for wall. This size must
            correpsond to the size of the heat flux obtained using the
            same groupName.

        groupName : str

            Group identifier to set only temperatures corresponding to
            the desired group. The group must be a family or a
            user-supplied group of families. The default is None which
            corresponds to all wall-type surfaces.

        TS : int
            Time spectral instance to set.
        """
        if groupName is None:
            groupName = self.allWallsGroup

        # For the mapVector to work correctly, we need to retrieve the
        # existing values and just overwrite the ones we've changed
        # using mapVector.
        npts, ncell = self._getSurfaceSize(self.allWallsGroup)
        fullTemp = self.adflow.gettnswall(npts, TS + 1)

        # Now map new values in and set.
        fullTemp = self.mapVector(temperature, groupName, self.allWallsGroup, fullTemp)
        self.adflow.settnswall(fullTemp, TS + 1)

    def setTargetCp(self, CpTargets, groupName=None, TS=0):
        """Set the CpTarget distribution for am inverse design problem.

        Parameters
        ----------
        CpSurf :  numpy array

            Array of CP targets to set for the surface. This size must
            correpsond to the size of the surface obtained using the
            same groupName.

        groupName : str

            Group identifier to set only CPtargets corresponding to
            the desired group. The group must be a family or a
            user-supplied group of families. The default is None which
            corresponds to all wall-type surfaces.

        TS : int
            Time spectral instance to set.
        """
        if groupName is None:
            groupName = self.allWallsGroup

        # For the mapVector to work correctly, we need to retrieve the
        # existing values and just overwrite the ones we've changed
        # using mapVector.
        npts, ncell = self._getSurfaceSize(self.allWallsGroup)
        fullCpTarget = numpy.atleast_2d(self.adflow.getcptargets(npts, TS + 1))

        # Now map new values in and set.
        fullCpTarget = self.mapVector(numpy.atleast_2d(CpTargets).T, groupName, self.allWallsGroup, fullCpTarget.T)
        fullCpTarget = fullCpTarget.T
        self.adflow.setcptargets(numpy.ravel(fullCpTarget), TS + 1)

    def getSurfacePoints(self, groupName=None, includeZipper=True, TS=0):
        """Return the coordinates for the surfaces defined by groupName.

        Parameters
        ----------
        groupName : str
            Group identifier to get only coordinates cooresponding to
            the desired group. The group must be a family or a
            user-supplied group of families. The default is None which
            corresponds to all wall-type surfaces.

        TS : int
            The time spectral instance to use for the forces.
        """

        if groupName is None:
            groupName = self.allWallsGroup

        if includeZipper:
            self._createZipperMesh()

        # Get the required size
        npts, ncell = self._getSurfaceSize(groupName, includeZipper)
        pts = numpy.zeros((npts, 3), self.dtype)
        famList = self._getFamilyList(groupName)
        if npts == 0:
            dummy = numpy.zeros((1, 3), self.dtype)
            self.adflow.surfaceutils.getsurfacepoints(dummy.T, TS + 1, famList, includeZipper)
        else:
            self.adflow.surfaceutils.getsurfacepoints(pts.T, TS + 1, famList, includeZipper)

        return pts

    def getSurfaceConnectivity(self, groupName=None, includeZipper=True, includeCGNS=False):
        """Return the connectivity dinates at which the forces (or tractions) are
        defined. This is the complement of getForces() which returns
        the forces at the locations returned in this routine.

        Parameters
        ----------
        groupName : str
            Group identifier to get only forces cooresponding to the
            desired group. The group must be a family or a user-supplied
            group of families. The default is None which corresponds to
            all wall-type surfaces.

        includeCGNS : bool
            Whether or not this function should return the indices of the CGNS blocks
            that each face patch belongs to. Zipper mesh patches will have cgnsBlockID = -1.

        """

        # Create the zipper mesh if we need it
        if includeZipper:
            self._createZipperMesh()

        if groupName is None:
            groupName = self.allWallsGroup

        # Set the list of surfaces this family requires
        famList = self._getFamilyList(groupName)
        npts, ncell = self._getSurfaceSize(groupName, includeZipper)
        conn = numpy.zeros((ncell, 4), dtype="intc")
        cgnsBlockID = numpy.zeros(max(1, ncell), dtype="intc")  # f2py will crash if we give a vector of length zero
        self.adflow.surfaceutils.getsurfaceconnectivity(
            numpy.ravel(conn), numpy.ravel(cgnsBlockID), famList, includeZipper
        )

        faceSizes = 4 * numpy.ones(len(conn), "intc")

        # Fix cgnsBlockID size if its length should be zero
        if ncell == 0:
            cgnsBlockID = numpy.zeros(ncell, dtype="intc")

        if includeCGNS:
            # Convert to 0-based ordering becuase we are in python
            return conn - 1, faceSizes, cgnsBlockID - 1
        else:
            # Convert to 0-based ordering becuase we are in python
            return conn - 1, faceSizes

    def _expandGroupNames(self, groupNames):
        """Take a list of family (group) names and return a 2D array of the
        fortran group numbers of the following form:

        [ 2 | y y x x x ]
        [ 1 | y x x x x ]
        [ 5 | y y y y y ]

        The first column is the number of significant family entries.
         The 'y' variables are the actual group numbers and the 'x'
         values are irrelevant.
        """

        nMax = 0
        for groupName in groupNames:
            nMax = max(nMax, len(self._getFamilyList(groupName)))

        groupArray = numpy.zeros((len(groupNames), nMax + 1))
        for i in range(len(groupNames)):
            famList = self._getFamilyList(groupNames[i])
            groupArray[i, 0] = len(famList)
            groupArray[i, 1 : 1 + len(famList)] = famList

        return groupArray

    def globalNKPreCon(self, inVec, outVec):
        """This function is ONLY used as a preconditioner to the
        global Aero-Structural system. This computes outVec =
        M^(-1)*inVec where M^(-1) is the approximate inverse
        application of the preconditing matrix.

        Parameters
        ----------
        inVec : array
            inVec must be size self.getStateSize()

        Returns
        -------
        outVec : array
            Preconditioned vector
        """
        return self.adflow.nksolver.applypc(inVec, outVec)

    def globalAdjointPreCon(self, inVec, outVec):
        """This function is ONLY used as a preconditioner to the
        global Aero-Structural ADJOINT system. This computes outVec =
        M^(-1)*inVec where M^(-1) is the approximate inverse
        application of the preconditing matrix.

        Parameters
        ----------
        inVec : array
            inVec must be size self.getAdjointStateSize()

        Returns
        -------
        outVec : array
            Preconditioned vector
        """
        return self.adflow.nksolver.applyadjointpc(inVec, outVec)

    def _addAeroDV(self, dv):
        """Add a single design variable that ADflow knows about.

        Parameters
        ----------
        dv : str
            dv name. Must be in the self.possibleAeroDVs list
        """
        dv = dv.lower()
        if dv not in self.possibleAeroDVs:
            raise Error(
                "%s was not one of the possible AeroDVs. "
                "The complete list of DVs for ADflow is %s. " % (dv, repr(set(self.possibleAeroDVs.keys())))
            )

        if dv not in self.aeroDVs:
            # A new DV add it:
            self.aeroDVs.append(dv)

    def _setupAdjoint(self, reform=False):
        """
        Setup the data structures required to solve the adjoint problem
        """

        # Destroy the NKsolver to free memory -- Call this even if the
        # solver is not used...a safeguard check is done in Fortran
        self.adflow.nksolver.destroynksolver()
        self.adflow.anksolver.destroyanksolver()
        self._setAeroDVs()

        if not self.adjointSetup or reform:
            # Create any PETSc variables if necessary
            self.adflow.adjointapi.createpetscvars()

            # Setup all required matrices in forward mode. (possibly none
            # of them)
            self.adflow.adjointapi.setupallresidualmatricesfwd()
            self.adflow.adjointapi.setuppetscksp()
            self.adjointSetup = True

    def releaseAdjointMemory(self):
        """
        release the PETSc Memory that have been allocated
        """
        if self.adjointSetup:
            self.adflow.adjointutils.destroypetscvars()
            self.adjointSetup = False

    def solveAdjoint(self, aeroProblem, objective, forcePoints=None, structAdjoint=None, groupName=None):
        # Remind the user they are using frozen turbulence.
        if self.getOption("frozenTurbulence") and self.myid == 0:
            self.getOption("equationType").lower() == "rans" and ADFLOWWarning(
                "Turbulence is frozen!!! DERIVATIVES WILL BE WRONG!!! " "USE AT OWN RISK!!!"
            )

        # May be switching aeroProblems here
        self.setAeroProblem(aeroProblem)

        # Possibly setup adjoint matrices/vector as required
        self._setupAdjoint()

        # Check to see if the RHS Partials have been computed
        if objective not in self.curAP.adflowData.adjointRHS:
            RHS = self.computeJacobianVectorProductBwd(funcsBar=self._getFuncsBar(objective), wDeriv=True)
            self.curAP.adflowData.adjointRHS[objective] = RHS.copy()
        else:
            RHS = self.curAP.adflowData.adjointRHS[objective].copy()

        # Check to see if we need to agument the RHS with a structural
        # adjoint:
        if structAdjoint is not None and groupName is not None:
            phi = self.mapVector(structAdjoint, groupName, self.allWallsGroup)

            agument = self.computeJacobianVectorProductBwd(fBar=phi, wDeriv=True)
            RHS -= agument

        # Check if objective is python 'allocated':
        if objective not in self.curAP.adflowData.adjoints:
            self.curAP.adflowData.adjoints[objective] = numpy.zeros(self.getAdjointStateSize(), float)

        # Initialize the fail flag in this AP if it doesn't exist
        if not hasattr(self.curAP, "adjointFailed"):
            self.curAP.adjointFailed = False

        # Check for any previous adjoint failure. If any adjoint has failed
        # on this AP, there is no point in solving the reset, so continue
        # with psi set as zero
        if not (self.curAP.adjointFailed and self.getOption("skipafterfailedadjoint")):
            if self.getOption("restartAdjoint"):
                # Use the previous solution as the initial guess
                psi = self.curAP.adflowData.adjoints[objective]
            else:
                # Use a zero initial guess
                psi = numpy.zeros_like(self.curAP.adflowData.adjoints[objective])

            # Actually solve the adjoint system
            # psi is updated with the new solution
            self.adflow.adjointapi.solveadjoint(RHS, psi, True)

            # Now set the flags and possibly reset adjoint
            if self.adflow.killsignals.adjointfailed:
                self.curAP.adjointFailed = True
                # Reset stored adjoint if we want to skip
                if self.getOption("skipafterfailedadjoint"):
                    if self.comm.rank == 0:
                        ADFLOWWarning(
                            "Current adjoint failed to converge, so the remaining adjoints will be skipped. Use the checkAdjointFailure method to get the correct fail flag for sensitivity evaluations. Finally, if you want to solve the remaining adjoints regardless of the current failure, set skipAfterFailedAdjoint to False."
                        )
                    self.curAP.adflowData.adjoints[objective][:] = 0.0
                else:
                    if self.comm.rank == 0:
                        ADFLOWWarning(
                            "Current adjoint failed to converge. The partially converged solution will be used to compute the total derivatives. It is up to the user to check for adjoint failures and pass the correct failure flag to the optimizer using the checkAdjointFailure method."
                        )
                    self.curAP.adflowData.adjoints[objective] = psi
            else:
                self.curAP.adflowData.adjoints[objective] = psi
                self.curAP.adjointFailed = False

    def _processAeroDerivatives(self, dIda, dIdBC):
        """This internal furncion is used to convert the raw array ouput from
        the matrix-free product bwd routine into the required
        dictionary format."""

        funcsSens = {}

        for dvName in self.curAP.DVs:
            key = self.curAP.DVs[dvName].key.lower()
            dvFam = self.curAP.DVs[dvName].family

            tmp = {}
            if key == "altitude":
                # This design variable is special. It combines changes
                # in temperature, pressure and density into a single
                # variable. Since we have derivatives for T, P and
                # rho, we simply chain rule it back to the the
                # altitude variable.
                self.curAP.evalFunctionsSens(tmp, ["P", "T", "rho"])

                # Extract the derivatives wrt the independent
                # parameters in ADflow
                dIdP = dIda[self.possibleAeroDVs["p"]]
                dIdT = dIda[self.possibleAeroDVs["t"]]
                dIdrho = dIda[self.possibleAeroDVs["rho"]]

                # Chain-rule to get the final derivative:
                funcsSens[dvName] = (
                    tmp[self.curAP["P"]][dvName] * dIdP
                    + tmp[self.curAP["T"]][dvName] * dIdT
                    + tmp[self.curAP["rho"]][dvName] * dIdrho
                )
            elif key == "mach":
                self.curAP.evalFunctionsSens(tmp, ["P", "rho"])
                # Simular story for Mach: It is technically possible
                # to use a mach number for a fixed RE simulation. For
                # the RE to stay fixed and change the mach number, the
                # 'P' and 'rho' must also change. We have to chain run
                # this dependence back through to the final mach
                # derivative. When Mach number is used with altitude
                # or P and T, this calc is unnecessary, but won't do
                # any harm.
                dIdP = dIda[self.possibleAeroDVs["p"]]
                dIdrho = dIda[self.possibleAeroDVs["rho"]]

                # Chain-rule to get the final derivative:
                funcsSens[dvName] = (
                    tmp[self.curAP["P"]][dvName] * dIdP
                    + tmp[self.curAP["rho"]][dvName] * dIdrho
                    + dIda[self.possibleAeroDVs["mach"]]
                )

            elif key in self.possibleAeroDVs:
                funcsSens[dvName] = dIda[self.possibleAeroDVs[key]]
                # Convert angle derivatives from 1/rad to 1/deg
                if key in ["alpha", "beta"]:
                    funcsSens[dvName] *= numpy.pi / 180.0

            elif key in self.possibleBCDvs:
                # We need to determine what the index is in dIdBC. For
                # now just do an efficient linear search:
                i = 0

                for tmp in self.curAP.bcVarData:
                    varName, family = tmp
                    if varName.lower() == key and family.lower() == dvFam.lower():
                        funcsSens[dvName] = dIdBC[i]
                    i += 1

        return funcsSens

    def _setAeroDVs(self):
        """Do everything that is required to deal with aerodynamic
        design variables in ADflow"""

        DVsRequired = list(self.curAP.DVs.keys())
        for dv in DVsRequired:
            key = self.curAP.DVs[dv].key.lower()
            if key in ["altitude"]:
                # All these variables need to be compined
                self._addAeroDV("P")
                self._addAeroDV("T")
                self._addAeroDV("rho")
            elif key in ["mach"]:
                self._addAeroDV("mach")
                self._addAeroDV("P")
                self._addAeroDV("rho")

            elif key in self.possibleAeroDVs:
                self._addAeroDV(key)
            elif key in self.possibleBCDvs:
                pass
            else:
                raise Error(
                    "The design variable '%s' as specified in the" " aeroProblem cannot be used with ADflow." % key
                )

    def solveAdjointForRHS(self, inVec, relTol=None):
        """
        Solve the adjoint system with an arbitary RHS vector.

        Parameters
        ----------
        inVec : numpy array
            Array of size w

        Returns
        -------
        outVec : numpy array
            Solution vector of size w
        """
        if relTol is None:
            relTol = self.getOption("adjointl2convergence")
        outVec = self.adflow.adjointapi.solveadjointforrhs(inVec, relTol)

        return outVec

    def solveDirectForRHS(self, inVec, relTol=None):
        """
        Solve the direct system with an arbitary RHS vector.

        Parameters
        ----------
        inVec : numpy array
            Array of size w

        Returns
        -------
        outVec : numpy array
            Solution vector of size w
        """
        if relTol is None:
            relTol = self.getOption("adjointl2convergence")
        outVec = self.adflow.adjointapi.solvedirectforrhs(inVec, relTol)

        return outVec

    def saveAdjointMatrix(self, baseFileName):
        """Save the adjoint matrix to a binary petsc file for
        possible future external testing

        Parameters
        ----------
        basefileName : str
            Filename to use. The Adjoint matrix, PC matrix(if it exists)
            and RHS  will be written
        """
        adjointMatrixName = baseFileName + "_drdw.bin"
        pcMatrixName = baseFileName + "_drdwPre.bin"
        self.adflow.adjointapi.saveadjointmatrix(adjointMatrixName)
        self.adflow.adjointapi.saveadjointpc(pcMatrixName)

    def saveAdjointRHS(self, baseFileName, objective):
        # Check to see if the RHS Partials have been computed
        if objective not in self.curAP.adflowData.adjointRHS:
            RHS = self.computeJacobianVectorProductBwd(funcsBar=self._getFuncsBar(objective), wDeriv=True)
            self.curAP.adflowData.adjointRHS[objective] = RHS.copy()
        else:
            RHS = self.curAP.adflowData.adjointRHS[objective].copy()

        rhsName = baseFileName + "_rhs_%s.bin" % objective
        self.adflow.adjointapi.saveadjointrhs(RHS, rhsName)

    def computeStabilityParameters(self):
        """
        run the stability derivative driver to compute the stability parameters
        from the time spectral solution
        """
        self.adflow.utils.stabilityderivativedriver()

    def updateGeometryInfo(self, warpMesh=True):
        """
        Update the ADflow internal geometry info.
        """

        # The mesh is modified
        if self._updateGeomInfo and self.mesh is not None:
            # If it is unsteady, and mesh is modified, then it has to be ALE
            if self.getOption("equationMode").lower() == "unsteady":
                self.adflow.preprocessingapi.shiftcoorandvolumes()
                self.adflow.aleutils.shiftlevelale()
            # Warp the mesh if surface coordinates are modified
            if warpMesh:
                timeA = time.time()
                self.mesh.warpMesh()
                newGrid = self.mesh.getSolverGrid()
                self.adflow.killsignals.routinefailed = False
                self.adflow.killsignals.fatalFail = False
                self.updateTime = time.time() - timeA
                if newGrid is not None and not self.getOption("useExternalDynamicMesh"):
                    # since updateGeometryInfo assumes one time slice
                    # when using ts with externally defined mesh, the grid is provided externally;
                    # updateGeometryInfo mainly serves to update the metrics and etc.
                    self.adflow.warping.setgrid(newGrid)
            # Update geometric data, depending on the type of simulation
            if self.getOption("equationMode").lower() == "unsteady":
                self.adflow.solvers.updateunsteadygeometry()
            else:
                self.adflow.preprocessingapi.updatecoordinatesalllevels()
                self.adflow.walldistance.updatewalldistancealllevels()
                self.adflow.preprocessingapi.updatemetricsalllevels()
                self.adflow.preprocessingapi.updategridvelocitiesalllevels()
                # Perform overset update
                ncells = self.adflow.adjointvars.ncellslocal[0]
                ntime = self.adflow.inputtimespectral.ntimeintervalsspectral
                n = ncells * ntime
                flag = numpy.zeros(n)

                # Only need to call the cutCallBack and regenerate the zipper mesh
                # if we're doing a full update.
                if self.getOption("oversetUpdateMode") == "full":
                    cutCallBack = self.getOption("cutCallBack")
                    if cutCallBack is not None:
                        xCen = self.adflow.utils.getcellcenters(1, n).T
                        cellIDs = self.adflow.utils.getcellcgnsblockids(1, n)
                        cutCallBack(xCen, self.CGNSZoneNameIDs, cellIDs, flag)

                    # Verify previous mesh failures
                    self.adflow.killsignals.routinefailed = self.comm.allreduce(
                        bool(self.adflow.killsignals.routinefailed), op=MPI.LOR
                    )
                    self.adflow.killsignals.fatalfail = self.adflow.killsignals.routinefailed

                    # We will need to update the zipper mesh later on.
                    # But we only do this if we have a valid mesh, otherwise the
                    # code might halt during zipper mesh regeneration
                    if not self.adflow.killsignals.fatalfail:
                        self.zipperCreated = False
                    else:  # We got mesh failure!
                        self.zipperCreated = True  # Set flag to true to skip zipper mesh call
                        if self.myid == 0:
                            print("ATTENTION: Zipper mesh will not be regenerated due to previous failures.")

                famList = self._getFamilyList(self.closedFamilyGroup)
                self.adflow.oversetapi.updateoverset(flag, famList)

            # Update flags
            self._updateGeomInfo = False
            self.adflow.killsignals.routinefailed = self.comm.allreduce(
                bool(self.adflow.killsignals.routinefailed), op=MPI.LOR
            )
            self.adflow.killsignals.fatalfail = self.adflow.killsignals.routinefailed

    def getAdjointResNorms(self):
        """
        Return the following adjoint residual norms:
        initRes Norm: Norm the adjoint RHS
        startRes Norm: Norm at the start of adjoint call (with possible non-zero restart)
        finalCFD Norm: Norm at the end of adjoint solve
        """
        initRes = self.adflow.adjointpetsc.adjresinit
        startRes = self.adflow.adjointpetsc.adjresstart
        finalRes = self.adflow.adjointpetsc.adjresfinal
        fail = self.adflow.killsignals.adjointfailed

        return initRes, startRes, finalRes, fail

    def getResNorms(self):
        """Return the initial, starting and final Res Norms. Typically
        used by an external solver."""
        return (
            numpy.real(self.adflow.iteration.totalr0),
            numpy.real(self.adflow.iteration.totalrstart),
            numpy.real(self.adflow.iteration.totalrfinal),
        )

    def setResNorms(self, initNorm=None, startNorm=None, finalNorm=None):
        """Set one of these norms if not None. Typlically used by an
        external solver"""
        if initNorm is not None:
            self.adflow.iteration.totalr0 = initNorm
        if startNorm is not None:
            self.adflow.iteration.totalrstart = startNorm
        if finalNorm is not None:
            self.adflow.iteration.finalNorm = finalNorm

    def _prescribedTSMotion(self):
        """Determine if we have prescribed motion timespectral analysis"""

        if (
            self.getOption("alphamode")
            or self.getOption("betamode")
            or self.getOption("machmode")
            or self.getOption("pmode")
            or self.getOption("qmode")
            or self.getOption("rmode")
            or self.getOption("altitudemode")
        ):
            return True
        else:
            return False

    def _isAeroObjective(self, objective):
        """
        This function takes in an external objective string and determines
        if the string is a valid function in ADflow. The function returns
        a boolean.

        Parameters
        ----------
        objective: string
            The objective to be tested.
        """
        if (
            objective.lower() in self.adflowCostFunctions.keys()
            or objective.lower() in self.adflowUserCostFunctions.keys()
        ):
            return True
        else:
            return False

    # =========================================================================
    #   The following routines two routines are the workhorse of the
    #   forward and reverse mode routines. These can compute *ANY*
    #   product that is possible from the solver.
    #   =========================================================================

    def computeJacobianVectorProductFwd(
        self,
        xDvDot=None,
        xSDot=None,
        xVDot=None,
        wDot=None,
        residualDeriv=False,
        funcDeriv=False,
        fDeriv=False,
        groupName=None,
        mode="AD",
        h=None,
        evalFuncs=None,
    ):
        """This the main python gateway for producing forward mode jacobian
        vector products. It is not generally called by the user by
        rather internally or from another solver. A DVGeo object and a
        mesh object must both be set for this routine.

        Parameters
        ----------
        xDvDot : dict
            Perturbation on the geometric design variables defined in DVGeo.
        xSDot : numpy array
            Perturbation on the surface
        xVDot : numpy array
            Perturbation on the volume
        wDot : numpy array
            Perturbation the state variables
        residualDeriv : bool
            Flag specifiying if the residualDerivative (dwDot) should be returned
        funcDeriv : bool
            Flag specifiying if the derviative of the cost functions
            (as defined in the current aeroproblem) should be returned.
        fderiv : bool
            Flag specifiying if the derviative of the surface forces (tractions)
            should be returned
        groupName : str
            Optional group name to use for evaluating functions. Defaults to all
            surfaces.
        mode : str ["AD", "FD", or "CS"]
            Specifies how the jacobian vector products will be computed.
        h : float
            Step sized used when the mode is "FD" or "CS

        Returns
        -------
        dwdot, funcsdot, fDot : array, dict, array
            One or more of the these are return depending on the ``*Deriv`` flags
        """

        if xDvDot is None and xSDot is None and xVDot is None and wDot is None:
            raise Error("computeJacobianVectorProductFwd: xDvDot, xSDot, xVDot and wDot cannot all be None")

        self._setAeroDVs()
        nTime = self.adflow.inputtimespectral.ntimeintervalsspectral

        # Default flags
        useState = False
        useSpatial = False

        # Process the Xs perturbation
        if xSDot is None:
            xsdot = numpy.zeros_like(self.coords0)
            xsdot = self.mapVector(xsdot, self.allFamilies, self.designFamilyGroup, includeZipper=False)
        else:
            xsdot = xSDot
            useSpatial = True

        # Process the Xv perturbation
        if xVDot is None:
            xvdot = numpy.zeros(self.getSpatialSize())
        else:
            xvdot = xVDot
            useSpatial = True

        # Process wDot perturbation
        if wDot is None:
            wdot = numpy.zeros(self.getStateSize())
        else:
            wdot = wDot
            useState = True

        # Process the extra variable perturbation....this comes from
        # xDvDot
        extradot = numpy.zeros(self.adflow.adjointvars.ndesignextra)
        bcDataNames, bcDataValues, bcDataFamLists, bcDataFams, bcVarsEmpty = self._getBCDataFromAeroProblem(self.curAP)
        bcDataValuesdot = numpy.zeros_like(bcDataValues)

        if xDvDot is not None:
            useSpatial = True
            for xKey in xDvDot:
                # We need to split out the family from the key.
                key = xKey.split("_")
                if len(key) == 1:
                    key = key[0].lower()
                    if key in self.possibleAeroDVs:
                        val = xDvDot[xKey]
                        if key.lower() == "alpha":
                            val *= numpy.pi / 180
                        extradot[self.possibleAeroDVs[key.lower()]] = val

                else:
                    fam = "_".join(key[1:])
                    key = key[0].lower()
                    if key in self.possibleBCDvs and not bcVarsEmpty:
                        # Figure out what index this should be:
                        for i in range(len(bcDataNames)):
                            if (
                                key.lower() == "".join(bcDataNames[i]).strip().lower()
                                and fam.lower() == bcDataFams[i].lower()
                            ):
                                bcDataValuesdot[i] = xDvDot[xKey]
        # For the geometric xDvDot perturbation we accumulate into the
        # already existing (and possibly nonzero) xsdot and xvdot
        if xDvDot is not None or xSDot is not None:
            if xDvDot is not None and self.DVGeo is not None:
                xsdot += self.DVGeo.totalSensitivityProd(xDvDot, self.curAP.ptSetName, config=self.curAP.name).reshape(
                    xsdot.shape
                )
            if self.mesh is not None:
                xsdot = self.mapVector(xsdot, self.meshFamilyGroup, self.designFamilyGroup, includeZipper=False)
                xvdot += self.mesh.warpDerivFwd(xsdot)
            useSpatial = True

        # Sizes for output arrays
        costSize = self.adflow.constants.ncostfunction
        fSize, nCell = self._getSurfaceSize(self.allWallsGroup, includeZipper=True)

        # process the functions and groups
        if evalFuncs is None:
            evalFuncs = sorted(self.curAP.evalFuncs)

        # Generate the list of families we need for the functions in curAP
        # We need to use a list to have the same ordering on every proc, but
        # we need 'set like' behavior so we don't include duplicate group names.
        groupNames = []
        for f in evalFuncs:
            fl = f.lower()
            if fl in self.adflowCostFunctions:
                groupName = self.adflowCostFunctions[fl][0]
                if groupName not in groupNames:
                    groupNames.append(groupName)
            if f in self.adflowUserCostFunctions:
                for sf in self.adflowUserCostFunctions[f].functions:
                    groupName = self.adflowCostFunctions[sf.lower()][0]
                    if groupName not in groupNames:
                        groupNames.append(groupName)

        if len(groupNames) == 0:
            famLists = self._expandGroupNames([self.allWallsGroup])
        else:
            famLists = self._expandGroupNames(groupNames)

        if mode == "AD":
            dwdot, tmp, fdot = self.adflow.adjointapi.computematrixfreeproductfwd(
                xvdot,
                extradot,
                wdot,
                bcDataValuesdot,
                useSpatial,
                useState,
                famLists,
                bcDataNames,
                bcDataValues,
                bcDataFamLists,
                bcVarsEmpty,
                costSize,
                max(1, fSize),
                nTime,
            )
        elif mode == "FD":
            if h is None:
                raise Error("if mode 'FD' is used, a stepsize must be specified using the kwarg 'h'")

            dwdot, tmp, fdot = self.adflow.adjointdebug.computematrixfreeproductfwdfd(
                xvdot,
                extradot,
                wdot,
                bcDataValuesdot,
                useSpatial,
                useState,
                famLists,
                bcDataNames,
                bcDataValues,
                bcDataFamLists,
                bcVarsEmpty,
                costSize,
                max(1, fSize),
                nTime,
                h,
            )
        elif mode == "CS":
            if h is None:
                raise Error("if mode 'CS' is used, a stepsize must be specified using the kwarg 'h'")

            if self.dtype == "D":
                dwdot, tmp, fdot = self.adflow.adjointdebug.computematrixfreeproductfwdcs(
                    xvdot,
                    extradot,
                    wdot,
                    bcDataValuesdot,
                    useSpatial,
                    useState,
                    famLists,
                    bcDataNames,
                    bcDataValues,
                    bcDataFamLists,
                    bcVarsEmpty,
                    costSize,
                    max(1, fSize),
                    nTime,
                    h,
                )
            else:
                raise Error("Complexified ADflow must be used to apply the complex step")
        else:
            raise Error(f"mode {mode} for computeJacobianVectorProductFwd not availiable")

        # Explictly put fdot to nothing if size is zero
        if fSize == 0:
            fdot = numpy.zeros((0, 3))

        # Process the derivative of the functions
        funcsDot = {}
        for f in evalFuncs:
            fl = f.lower()
            if fl in self.adflowCostFunctions:
                groupName = self.adflowCostFunctions[fl][0]
                basicFunc = self.adflowCostFunctions[fl][1]
                mapping = self.basicCostFunctions[basicFunc]
                # Get the group index:
                ind = groupNames.index(groupName)
                funcsDot[f] = tmp[mapping - 1, ind]
            elif f in self.adflowUserCostFunctions:
                # We must linearize the user-supplied function
                sens = self.adflowUserCostFunctions.evalFunctionsSens()
                funcsDot[f] = 0.0
                for sf in self.adflowUserCostFunctions[f].functions:
                    groupName = self.adflowCostFunctions[sf][0]
                    basicFunc = self.adflowCostFunctions[sf][1]
                    mapping = self.basicCostFunctions[basicFunc]
                    # Get the group index:
                    ind = groupNames.index(groupName)
                    funcsDot[f] += sens[sf] * tmp[mapping - 1, ind]

        # Assemble the returns
        returns = []
        if residualDeriv:
            returns.append(dwdot)
        if funcDeriv:
            returns.append(funcsDot)
        if fDeriv:
            returns.append(fdot.T)

        return tuple(returns) if len(returns) > 1 else returns[0]

    def computeJacobianVectorProductBwd(
        self,
        resBar=None,
        funcsBar=None,
        fBar=None,
        wDeriv=None,
        xVDeriv=None,
        xSDeriv=None,
        xDvDeriv=None,
        xDvDerivAero=None,
    ):
        """This the main python gateway for producing reverse mode jacobian
        vector products. It is not generally called by the user by
        rather internally or from another solver. A mesh object must
        be present for the ``xSDeriv=True`` flag and a mesh and DVGeo
        object must be present for ``xDvDeriv=True`` flag. Note that more
        than one of the specified return flags may be spcified. If
        more than one return is specified, the order of return is :
        ``(wDeriv, xVDeriv, XsDeriv, xDvDeriv, dXdvDerivAero)``.

        Parameters
        ----------
        resBar : numpy array
            Seed for the residuals (dwb in adflow)
        funcsBar : dict
            Dictionary of functions with reverse seeds. Only nonzero seeds
            need to be provided. All other seeds will be taken as zero.
        fBar : numpy array
            Seed for the forces (or tractions depending on the option value) to use
            in reverse mode.

        wDeriv : bool
            Flag specifiying if the state (w) derivative (wb) should be returned
        xVDeriv : bool
            Flag specifiying if the volume node (xV) derivative should be returned
        xSDeriv : bool
            Flag specifiying if the surface node (xS) derivative should be returned.
        xDvDeriv : bool
            Flag specifiying if the design variable (xDv) derviatives should
            be returned. This will include both geometric *and* aerodynamic derivatives
        xDvDerivAero : bool
            Flag to return *just* the aerodynamic derivatives. If this is True and
            xDvDeriv is False,*just* the aerodynamic derivatives are returned.

        Returns
        -------
        wbar, xvbar, xsbar, xdvbar, xdvaerobar : array, array, array, dict, dict
            One or more of these are returned depending on the ``*Deriv`` flags provided.

        """
        # Error Checking
        if resBar is None and funcsBar is None and fBar is None:
            raise Error(
                "computeJacobianVectorProductBwd: One of resBar, funcsBar and fBar"
                " must be given. resBar=%s, funcsBar=%s, fBar=%s" % (resBar, funcsBar, fBar)
            )
        if wDeriv is None and xVDeriv is None and xDvDeriv is None and xSDeriv is None and xDvDerivAero is None:
            raise Error(
                "computeJacobianVectorProductBwd: One of wDeriv, xVDeriv, "
                "xDvDeriv and xDvDerivAero must be given as True. "
                "wDeriv=%s, xVDeriv=%s, xDvDeriv=%s, xSDeriv=%s xDvDerivAero=%s."
                % (wDeriv, xVDeriv, xDvDeriv, xSDeriv, xDvDerivAero)
            )
        self._setAeroDVs()

        # ---------------------
        #  Check for resBar
        # ---------------------
        if resBar is None:
            if self.getOption("frozenTurbulence"):
                resBar = numpy.zeros(self.getAdjointStateSize())
            else:
                resBar = numpy.zeros(self.getStateSize())

        # -------------------------
        #  Check for fBar (forces)
        # ------------------------
        nTime = self.adflow.inputtimespectral.ntimeintervalsspectral
        nPts, nCell = self._getSurfaceSize(self.allWallsGroup, includeZipper=True)

        if fBar is None:
            fBar = numpy.zeros((nTime, nPts, 3))
        else:
            # Expand out to the sps direction in case there were only
            # 2 dimensions.
            fBar = fBar.reshape((nTime, nPts, 3))

        # ---------------------
        #  Check for funcsBar
        # ---------------------

        if funcsBar is None:
            funcsBar = {}

        # Generate the list of families:
        groupNames = []
        tmp = []

        for f in funcsBar:
            fl = f.lower()
            if fl in self.adflowCostFunctions:
                groupName = self.adflowCostFunctions[fl][0]
                basicFunc = self.adflowCostFunctions[fl][1]

                if groupName in groupNames:
                    ind = groupNames.index(groupName)
                else:
                    groupNames.append(groupName)
                    tmp.append(numpy.zeros(self.adflow.constants.ncostfunction))
                    ind = -1
                mapping = self.basicCostFunctions[basicFunc]
                tmp[ind][mapping - 1] += funcsBar[f]

        if len(tmp) == 0:
            # No adflow-functions given. Just set zeros. Same as
            # if funcsBar was None
            funcsBar = numpy.zeros((self.adflow.constants.ncostfunction, 1))
            famLists = self._expandGroupNames([self.allWallsGroup])
        else:
            famLists = self._expandGroupNames(groupNames)
            funcsBar = numpy.array(tmp).T

        # Determine if we can do a cheaper state-variable only
        # computation or if we need to include the spatial terms:
        useSpatial = False
        useState = False
        if wDeriv:
            useState = True
        if xDvDeriv or xVDeriv or xSDeriv or xDvDerivAero:
            useSpatial = True

        # Extract any possibly BC daa
        bcDataNames, bcDataValues, bcDataFamLists, bcDataFams, bcVarsEmpty = self._getBCDataFromAeroProblem(self.curAP)

        # Do actual Fortran call.
        xvbar, extrabar, wbar, bcdatavaluesbar = self.adflow.adjointapi.computematrixfreeproductbwd(
            resBar,
            funcsBar,
            fBar.T,
            useSpatial,
            useState,
            self.getSpatialSize(),
            self.adflow.adjointvars.ndesignextra,
            self.getAdjointStateSize(),
            famLists,
            bcDataNames,
            bcDataValues,
            bcDataFamLists,
            bcVarsEmpty,
        )

        # Assemble the possible returns the user has requested:
        returns = []
        if wDeriv:
            returns.append(wbar)
        if xVDeriv:
            returns.append(xvbar)

        # Process xVbar back to the xS or xDV (geometric variables) if necessary
        if xDvDeriv or xSDeriv:
            if self.mesh is not None:
                self.mesh.warpDeriv(xvbar)
                xsbar = self.mesh.getdXs()
                xsbar = self.mapVector(xsbar, self.meshFamilyGroup, self.designFamilyGroup, includeZipper=False)

                if xSDeriv:
                    returns.append(xsbar)
            else:
                # Only an error if xSDeriv is requested...since we
                # can't do it. xDVDeriv may be specified even when no
                # mesh is present.
                if xSDeriv:
                    raise Error("Could not complete requested xSDeriv " "derivatives since no mesh is present")

            # Process all the way back to the DVs:
            if xDvDeriv:
                xdvbar = {}
                if self.mesh is not None:  # Include geometric
                    # derivatives if mesh is
                    # present
                    if self.DVGeo is not None and self.DVGeo.getNDV() > 0:
                        xdvbar.update(
                            self.DVGeo.totalSensitivity(xsbar, self.curAP.ptSetName, self.comm, config=self.curAP.name)
                        )
                    else:
                        if self.comm.rank == 0:
                            ADFLOWWarning(
                                "No DVGeo object is present or it has no "
                                "design variables specified. No geometric "
                                "derivatives computed."
                            )
                else:
                    if self.comm.rank == 0:
                        ADFLOWWarning("No mesh object is present. No geometric " "derivatives computed.")

                # Include aero derivatives here:
                xdvbar.update(self._processAeroDerivatives(extrabar, bcdatavaluesbar))
                returns.append(xdvbar)

        # Include the aerodynamic variables if requested to do so and
        # we haven't already done so with xDvDeriv:
        if xDvDerivAero and not xDvDeriv:
            xdvaerobar = {}
            xdvaerobar.update(self._processAeroDerivatives(extrabar, bcdatavaluesbar))
            returns.append(xdvaerobar)

        # Single return (most frequent) is 'clean', otherwise a tuple.
        return tuple(returns) if len(returns) > 1 else returns[0]

    def computeJacobianVectorProductBwdFast(
        self,
        resBar=None,
    ):
        """
        This fast routine computes only the derivatives of the residuals with respect to the states.
        This is the operator used for the matrix-free solution of the adjoint system.

        Parameters
        ----------
        resBar : numpy array
            Seed for the residuals (dwb in adflow)

        Returns
        -------
        wbar: array
            state derivative seeds
        """
        wbar = self.adflow.adjointapi.computematrixfreeproductbwdfast(resBar, self.getAdjointStateSize())

        return wbar

    def mapVector(self, vec1, groupName1, groupName2, vec2=None, includeZipper=True):
        """This is the main workhorse routine of everything that deals with
        families in ADflow. The purpose of this routine is to convert a
        vector 'vec1' (of size Nx3) that was evaluated with
        'groupName1' and expand or contract it (and adjust the
        ordering) to produce 'vec2' evaluated on groupName2.

        A little ascii art might help. Consider the following "mesh"
        . Family 'fam1' has 9 points, 'fam2' has 10 pts and 'fam3' has
        5 points.  Consider that we have also also added two
        additional groups: 'f12' containing 'fam1' and 'fam2' and a
        group 'f23' that contains families 'fam2' and 'fam3'. The vector
        we want to map is 'vec1'. It is length 9+10. All the 'x's are
        significant values.

        The call: mapVector(vec1, 'f12', 'f23')

        will produce the "returned vec" array, containing the
        significant values from 'fam2', where the two groups overlap,
        and the new values from 'fam3' set to zero. The values from
        fam1 are lost. The returned vec has size 15.

        ::

                fam1     fam2      fam3
            |---------+----------+------|

            |xxxxxxxxx xxxxxxxxxx|        <- vec1
                      |xxxxxxxxxx 000000| <- returned vec (vec2)

        It is also possible to pass in vec2 into this routine. For
        that case, the existing values in the array will not be
        kept. In the previous examples, the values corresponding to
        fam3 will retain their original values.

        Parameters
        ----------
        vec1 : Numpy array
            Array of size Nx3 that will be mapped to a different family set.

        groupName1 : str
            The family group where the vector vec1 is currently defined

        groupName2 : str
            The family group where we want to the vector to mapped into

        vec2 : Numpy array or None
            Array containing existing values in the output vector we want to keep.
            If this vector is not given, the values will be filled with zeros.

        Returns
        -------
        vec2 : Numpy array
            The input vector maped to the families defined in groupName2.

        """

        if groupName1 not in self.families or groupName2 not in self.families:
            raise Error(
                "'%s' or '%s' is not a family in the CGNS file or has not been added"
                " as a combination of families" % (groupName1, groupName2)
            )

        # Why can't we do the short-cut listed below? Well, it may
        # happen that we have the same family group but one has the
        # zipper nodes and the other doesn't. In that case, vec2 must
        # be a differnet size and thus we can't immediately return
        # vec1 as being vec2.

        # if groupName1 == groupName2:
        #     return vec1

        if vec2 is None:
            npts, ncell = self._getSurfaceSize(groupName2, includeZipper)
            vec2 = numpy.zeros((npts, 3), self.dtype)

        famList1 = self.families[groupName1]
        famList2 = self.families[groupName2]

        self.adflow.surfaceutils.mapvector(vec1.T, famList1, vec2.T, famList2, includeZipper)

        return vec2

    def getStateSize(self):
        """Return the number of degrees of freedom (states) that are on this
        processor. This is (number of states)*(number of
        cells)*(number of spectral instances)

        """

        nstate = self.adflow.flowvarrefstate.nw
        ncells = self.adflow.adjointvars.ncellslocal[0]
        ntime = self.adflow.inputtimespectral.ntimeintervalsspectral

        return nstate * ncells * ntime

    def getAdjointStateSize(self):
        """Return the number of ADJOINT degrees of freedom (states)
        that are on this processor. The reason this is different from
        getStateSize() is that if frozenTurbulence is used for RANS,
        the nonlinear system has 5+neq turb states per cell, while the
        adjoint still has 5."""
        if self.getOption("frozenTurbulence"):
            nstate = self.adflow.flowvarrefstate.nwf
        else:
            nstate = self.adflow.flowvarrefstate.nw

        ncells = self.adflow.adjointvars.ncellslocal[0]
        ntime = self.adflow.inputtimespectral.ntimeintervalsspectral

        return nstate * ncells * ntime

    def getSpatialSize(self):
        """Return the number of degrees of spatial degrees of freedom
        on this processor. This is (number of nodes)*(number of
        spectral instances)*3"""

        nnodes = self.adflow.adjointvars.nnodeslocal[0]
        ntime = self.adflow.inputtimespectral.ntimeintervalsspectral

        return 3 * nnodes * ntime

    def getNCells(self):
        """Return the total number of cells in the mesh"""

        return self.adflow.adjointvars.ncellsglobal[0]

    def getStates(self):
        """Return the states on this processor. Used in aerostructural
        analysis"""

        return self.adflow.nksolver.getstates(self.getStateSize())

    def setStates(self, states):
        """Set the states on this processor. Used in aerostructural analysis
        and for switching aeroproblems
        """
        self.adflow.nksolver.setstates(states)

    def getSurfacePerturbation(self, seed=314):
        """This is is a debugging routine only. It is used only in regression
        tests when it is necessary to compute a consistent random
        surface perturbation seed that is independent of per-processor
        block distribution.

        Parameters
        ----------
        seed : integer
            Seed to use for random number. Only significant on root processor

        """
        nPts, nCell = self._getSurfaceSize(self.allWallsGroup, includeZipper=True)
        xRand = self.getSpatialPerturbation(seed)
        surfRand = numpy.zeros((nPts, 3))
        famList = self._getFamilyList(self.allWallsGroup)
        # Only coded for 1 spectal instance currently.
        self.adflow.warping.getsurfaceperturbation(xRand, numpy.ravel(surfRand), famList, 1)
        return surfRand

    def getStatePerturbation(self, seed=314):
        """This is is a debugging routine only. It is used only in regression
        tests when it is necessary to compute a consistent random
        state vector seed that is independent of per-processor block
        distribution. This routine is *not* memory scalable as a
        complete state vector is generated on each process.

        Parameters
        ----------
        seed : integer
            Seed to use for random number. Only significant on root processor
        """

        # Get the total number of DOF
        totalDOF = self.comm.reduce(self.getStateSize())
        numpy.random.seed(seed)
        randVec = None
        if self.comm.rank == 0:
            randVec = numpy.random.random(totalDOF)

        randVec = self.comm.bcast(randVec)

        return self.adflow.warping.getstateperturbation(randVec, self.getStateSize())

    def getSpatialPerturbation(self, seed=314):
        """This is is a debugging routine only. It is used only in regression
        tests when it is necessary to compute a consistent random
        spatial vector seed that is independent of per-processor block
        distribution.

        Parameters
        ----------
        seed : integer
            Seed to use for random number. Only significant on root processor
        """

        # This routine is *NOT SCALABLE*. It requires storing the full
        # mesh on each processor. It should only be used for
        # verification purposes.

        # Get all CGNS Mesh indices to all Proc.
        localIndices = self.adflow.warping.getcgnsmeshindices(self.getSpatialSize())
        cgnsIndices = numpy.hstack(self.comm.allgather(localIndices))

        # Gather all nodes to all procs.
        pts = self.adflow.warping.getgrid(self.getSpatialSize())
        allPts = numpy.hstack(self.comm.allgather(pts))

        # Also need the point offset.s.
        ptSizes = self.comm.allgather(len(pts))
        offsets = numpy.zeros(len(ptSizes), "intc")
        offsets[1:] = numpy.cumsum(ptSizes)[:-1]

        # Now Re-assemble the global CGNS vector.
        nDOFCGNS = numpy.max(cgnsIndices) + 1
        CGNSVec = numpy.zeros(nDOFCGNS, dtype=self.dtype)
        CGNSVec[cgnsIndices] = allPts
        CGNSVec = CGNSVec.reshape((len(CGNSVec) // 3, 3))

        # Run the pointReduce on the CGNS nodes
        uniquePts, linkTmp, nUnique = self.adflow.utils.pointreduce(CGNSVec.T, 1e-12)

        # Expand link out to the 3x the size and convert to 1 based ordering
        linkTmp -= 1
        link = numpy.zeros(len(linkTmp) * 3, "intc")
        link[0::3] = 3 * linkTmp
        link[1::3] = 3 * linkTmp + 1
        link[2::3] = 3 * linkTmp + 2

        # Set the seed and everyone produces the random vector for
        # nUnique pts.
        numpy.random.seed(seed)
        randVec = numpy.random.random(nUnique * 3)

        # Finally just extract out our own part:
        iStart = offsets[self.comm.rank]
        iEnd = iStart + self.getSpatialSize()
        indices = numpy.arange(iStart, iEnd)
        spatialPerturb = randVec[link[cgnsIndices[indices]]]

        return spatialPerturb

    def getUniqueSpatialPerturbationNorm(self, dXv):
        """This is is a debugging routine only. It is used only in regression
        tests when it is necessary to compute the norm of a spatial
        perturbuation on meshes that are split. This will unique-ify
        the nodes and accumulate onto the unique nodes thus giving the
        same norm independent of the block splits. Again, this routine
        is not memory scalable and should only be used for debugging
        purposes.

        Parameters
        ----------
        dXv : numpy vector
            Spatial perturbation of size getSpatialSize()
        """

        # Gather all nodes to the root proc:
        pts = self.adflow.warping.getgrid(self.getSpatialSize())
        allPts = numpy.hstack(self.comm.allgather(pts))
        dXv = numpy.hstack(self.comm.allgather(dXv))
        norm = None
        if self.myid == 0:
            allPts = allPts.reshape((len(allPts) // 3, 3))
            dXv = dXv.reshape((len(dXv) // 3, 3))
            # Run the pointReduce on all nodes
            uniquePts, link, nUnique = self.adflow.utils.pointreduce(allPts.T, 1e-12)
            uniquePtsBar = numpy.zeros((nUnique, 3))
            link = link - 1  # Convert to zero-based for python:

            for i in range(len(link)):
                uniquePtsBar[link[i]] += dXv[i]

            # You might be tempted to vectorize the loop above as:
            # uniquePtsBar[link] += dXv
            # But you would be wrong. It does not work when link
            # references multiple indices more than once.

            # Now just take the norm of uniquePtsBar. The flatten
            # isn't strictly necessary.
            norm = numpy.linalg.norm(uniquePtsBar.flatten())

        return self.comm.bcast(norm)

    def _getInfo(self):
        """Get the haloed state vector, pressure (and viscocities). Used to
        save "state" between aeroProblems

        """
        return self.adflow.nksolver.getinfo(self.adflow.nksolver.getinfosize())

    def _setInfo(self, info):
        """Set the haloed state vector, pressure (and viscocities). Used to
        restore "state" between aeroProblems
        """
        self.adflow.nksolver.setinfo(info)

    def setAdjoint(self, adjoint, objective=None):
        """Sets the adjoint vector externally. Used in coupled solver"""
        if objective is not None:
            self.curAP.adflowData.adjoints[objective] = adjoint.copy()

    def getAdjoint(self, objective):
        """Return the adjoint values for objective if they
        exist. Otherwise just return zeros"""

        if objective in self.curAP.adflowData.adjoints:
            return self.curAP.adflowData.adjoints[objective]
        else:
            return numpy.zeros(self.getAdjointStateSize(), self.dtype)

    def getResidual(self, aeroProblem, res=None, releaseAdjointMemory=True):
        """Return the residual on this processor. Used in aerostructural
        analysis"""
        self.setAeroProblem(aeroProblem, releaseAdjointMemory)
        if res is None:
            res = numpy.zeros(self.getStateSize())
        res = self.adflow.nksolver.getres(res)

        return res

    def solveErrorEstimate(self, aeroProblem, funcError, evalFuncs=None):
        r"""
        Evaluate the desired function errors given in iterable object
        'evalFuncs', and add them to the dictionary 'funcError'. The keys
        in the funcError dictionary will be have an ``<ap.name>_`` prepended to them.

        Parameters
        ----------
        aeroProblem : :class:`~baseclasses:baseclasses.problems.pyAero_problem.AeroProblem`
            The aerodynamic problem to get the error for

        funcError : dict
            Dictionary into which the function errors are saved.
            We define error to be :math:`\epsilon = f^\ast - f`, where
            :math:`f^\ast` is the converged solution and :math:`f` is the unconverged solution.

        evalFuncs : iterable object containing strings
            If not None, use these functions to evaluate.

        Examples
        --------
        >>> CFDsolver(ap)
        >>> funcsSens = {}
        >>> CFDSolver.evalFunctionsSens(ap, funcsSens)
        >>> funcError = {}
        >>> CFDsolver.solveErrorEstimate(ap, funcError)
        >>> # Result will look like (if aeroProblem, ap, has name of 'wing'):
        >>> print(funcError)
        >>> # {'wing_cl':0.00085, 'wing_cd':0.000021}
        """
        self.setAeroProblem(aeroProblem)
        res = self.getResidual(aeroProblem)
        if evalFuncs is None:
            evalFuncs = sorted(self.curAP.evalFuncs)

        # Make sure we have a list that has only lower-cased entries
        tmp = []
        for f in evalFuncs:
            tmp.append(f.lower())
        evalFuncs = tmp

        # Do the functions one at a time:
        for f in evalFuncs:
            key = f"{self.curAP.name}_{f}"

            # Set dict structure for this derivative
            psi = self.getAdjoint(f)
            error = numpy.dot(psi, res)
            errorTot = self.comm.allreduce(error)
            errorTot = -1 * errorTot  # multiply by -1 so that estimate is added to current output
            funcError[key] = errorTot

    def getFreeStreamResidual(self, aeroProblem):
        self.setAeroProblem(aeroProblem)
        rhoRes, totalRRes = self.adflow.nksolver.getfreestreamresidual()
        return totalRRes

    def _getSurfaceSize(self, groupName, includeZipper=True):
        """Internal routine to return the size of a particular surface. This
        does *NOT* set the actual family group"""
        if groupName is None:
            groupName = self.allFamilies
        if groupName not in self.families:
            raise Error(
                "'%s' is not a family in the CGNS file or has not been added"
                " as a combination of families" % groupName
            )

        [nPts, nCells] = self.adflow.surfaceutils.getsurfacesize(self.families[groupName], includeZipper)
        return nPts, nCells

    def rootSetOption(self, name, value, reset=False):
        """
        Set ADflow options from the root proc.
        The option will be broadcast to all procs when __call__ is invoked.

        Parameters
        ----------
        name : str
            The name of the option
        value : Any
            The value of the option
        reset : Bool
            If True, we reset all previously-set rootChangedOptions.

        See Also
        --------
        :func:`setOption`
        """
        if reset:
            self.rootChangedOptions = {}
        self.rootChangedOptions[name] = value

    def setOption(self, name, value):
        """
        Set Solver Option Value
        """
        super().setOption(name, value)
        name = name.lower()

        # If the option is only used in Python, we just return
        if name in self.pythonOptions:
            return

        # Do special Options individually
        if name in self.specialOptions:
            if name in ["monitorvariables", "surfacevariables", "volumevariables", "isovariables"]:
                varStr = ""
                for i in range(len(value)):
                    varStr = varStr + value[i] + "_"
                # end if
                varStr = varStr[0:-1]  # Get rid of last '_'
                if name == "monitorvariables":
                    self.adflow.inputparamroutines.monitorvariables(varStr)
                if name == "surfacevariables":
                    self.adflow.inputparamroutines.surfacevariables(varStr)
                if name == "volumevariables":
                    self.adflow.inputparamroutines.volumevariables(varStr)
                if name == "isovariables":
                    self.adflow.inputparamroutines.isovariables(varStr)

            elif name == "restartfile":
                # If value is None no value has been specified by the
                # user. None is the default value.
                if value is not None:
                    # Check its type, if a string its a single value,
                    # but if list multiple
                    if type(value) is str:
                        # Check empty string
                        if value:
                            # Allocate only one slot since we have
                            # only one filename
                            self.adflow.initializeflow.allocrestartfiles(1)
                            self.adflow.initializeflow.setrestartfiles(value, 1)
                        else:
                            # Empty string. Raise error
                            raise Error(
                                "Option 'restartfile' string cannot be empty. "
                                "If not performing a restart, 'restartFile' "
                                "option must be set to None"
                            )
                    elif type(value) is list:
                        # Check input
                        nFiles = len(value)
                        if nFiles > 0:
                            if type(value[0]) is str:
                                # Allocate for the entire list
                                self.adflow.initializeflow.allocrestartfiles(nFiles)
                                # Populate the array
                                for i, val in enumerate(value):
                                    # The +1 is to match fortran indexing
                                    self.adflow.initializeflow.setrestartfiles(val, i + 1)
                            else:
                                raise Error(
                                    "Datatype for Option %-35s was not "
                                    "valid. Expected list of <type 'str'>. "
                                    "Received data type is "
                                    "%-47s" % (name, type(value[0]))
                                )
                        else:
                            raise Error(
                                "Option %-35s of %-35s contains %-35s "
                                "elements. Must contain at least 1 "
                                "restart file of type <type "
                                "'str'>" % (name, type(value), len(value))
                            )
                    else:
                        raise Error(
                            "Datatype for Option %-35s not valid. "
                            "Expected data type is <type 'str'> or <type "
                            "'list'>. Received data type is %-47s" % (name, type(value))
                        )

            elif name == "isosurface":
                # We have a bit of work to do...extract out the
                # names, and there can be more than 1 value per variables
                var = []
                val = []
                isoDict = value
                for key in isoDict.keys():
                    isoVals = numpy.atleast_1d(isoDict[key])
                    for i in range(len(isoVals)):
                        var.append(key)
                        val.append(isoVals[i])

                val = numpy.array(val)

                self.adflow.inputparamroutines.initializeisosurfacevariables(val)
                for i in range(len(val)):
                    self.adflow.inputparamroutines.setisosurfacevariable(var[i], i + 1)

            elif name == "turbresscale":
                # If value is None no value has been specified by the
                # user. None is the default value.  Do nothing as it
                # will be updated with _updateTurbResScale from
                # __init__
                if value is not None:
                    tmp_turbresscalar = [0.0, 0.0, 0.0, 0.0]
                    # Check type to handle the insert properly
                    if type(value) is float:
                        tmp_turbresscalar[0] = value
                    elif type(value) is list:
                        if 1 <= len(value) and len(value) <= 4:
                            if type(value[0]) is float:
                                tmp_turbresscalar[0 : len(value)] = value[:]
                            else:
                                raise Error(
                                    "Datatype for Option %-35s was not "
                                    "valid. Expected list of "
                                    "<type 'float'>. Received data type "
                                    "is %-47s" % (name, type(value[0]))
                                )
                        else:
                            raise Error(
                                "Option %-35s of %-35s contains %-35s "
                                "elements. Min and max number of "
                                "elements are 1 and 4 "
                                "respectively" % (name, type(value), len(value))
                            )
                    else:
                        raise Error(
                            "Datatype for Option %-35s not valid. Expected "
                            "data type is <type 'float'> or <type "
                            "'list'>. Received data type is %-47s" % (name, type(value))
                        )

                    module = self.moduleMap[self.optionMap[name][0]]
                    variable = self.optionMap[name][1]
                    setattr(module, variable, tmp_turbresscalar)
            elif name == "oversetpriority":
                # Loop over each of the block names and call the fortran setter:
                for blkName in value:
                    setValue = self.adflow.oversetapi.setblockpriority(blkName, value[blkName])
                    if not setValue and self.myid == 0:
                        ADFLOWWarning(
                            "The block name %s was not found in the CGNS file "
                            "and could not set it's priority" % blkName
                        )

            # Special option has been set so return from function
            return

        # All other options do genericaly by setting value in module:
        # Check if there is an additional mapping to what actually
        # has to be set in the solver

        if isinstance(self.optionMap[name], dict):
            module = self.moduleMap[self.optionMap[name]["location"][0]]
            variable = self.optionMap[name]["location"][1]
            value = self.optionMap[name][value.lower()]
        else:
            module = self.moduleMap[self.optionMap[name][0]]
            variable = self.optionMap[name][1]

        # If the value is a string, pads additional spaces
        if isinstance(value, str):
            spacesToAdd = self.adflow.constants.maxstringlen - len(value)
            value = "".join([value, " " * spacesToAdd])

        # Set in the correct module
        setattr(module, variable, value)

    @staticmethod
    def _getDefaultOptions():
        defOpts = {
            # Input file parameters
            "gridFile": [str, "default.cgns"],
            "restartFile": [(str, list, type(None)), None],
            # Surface definition parameters
            "meshSurfaceFamily": [(str, type(None)), None],
            "designSurfaceFamily": [(str, type(None)), None],
            "closedSurfaceFamilies": [(list, type(None)), None],
            # Output Parameters
            "storeRindLayer": [bool, True],
            "outputDirectory": [str, "./"],
            "outputSurfaceFamily": [str, "allSurfaces"],
            "writeSurfaceSolution": [bool, True],
            "writeVolumeSolution": [bool, True],
            "writeSolutionEachIter": [bool, False],
            "writeTecplotSurfaceSolution": [bool, False],
            "nSaveVolume": [int, 1],
            "nSaveSurface": [int, 1],
            "solutionPrecision": [str, ["single", "double"]],
            "gridPrecision": [str, ["double", "single"]],
            "solutionPrecisionSurface": [str, ["single", "double"]],
            "gridPrecisionSurface": [str, ["single", "double"]],
            "isosurface": [dict, {}],
            "isoVariables": [list, []],
            "viscousSurfaceVelocities": [bool, True],
            # Physics Parameters
            "discretization": [str, ["central plus scalar dissipation", "central plus matrix dissipation", "upwind"]],
            "coarseDiscretization": [
                str,
                ["central plus scalar dissipation", "central plus matrix dissipation", "upwind"],
            ],
            "limiter": [str, ["van Albada", "minmod", "no limiter"]],
            "smoother": [str, ["DADI", "Runge-Kutta"]],
            "equationType": [str, ["RANS", "Euler", "laminar NS"]],
            "equationMode": [str, ["steady", "unsteady", "time spectral"]],
            "flowType": [str, ["external", "internal"]],
            "turbulenceModel": [
                str,
                ["SA", "SA-Edwards", "k-omega Wilcox", "k-omega modified", "k-tau", "Menter SST", "v2f"],
            ],
            "turbulenceOrder": [str, ["first order", "second order"]],
            "turbResScale": [(float, list, type(None)), None],
            "meshMaxSkewness": [float, 1.0],
            "useSkewnessCheck": [bool, False],
            "turbulenceProduction": [str, ["strain", "vorticity", "Kato-Launder"]],
            "useQCR": [bool, False],
            "useRotationSA": [bool, False],
            "useft2SA": [bool, True],
            "eddyVisInfRatio": [float, 0.009],
            "useWallFunctions": [bool, False],
            "useApproxWallDistance": [bool, True],
            "eulerWallTreatment": [
                str,
                [
                    "linear pressure extrapolation",
                    "constant pressure extrapolation",
                    "quadratic pressure extrapolation",
                    "normal momentum",
                ],
            ],
            "viscWallTreatment": [str, ["constant pressure extrapolation", "linear pressure extrapolation"]],
            "dissipationScalingExponent": [float, 0.67],
            "acousticScaleFactor": [float, 1.0],
            "vis4": [float, 0.0156],
            "vis2": [float, 0.25],
            "vis2Coarse": [float, 0.5],
            "restrictionRelaxation": [float, 0.80],
            "liftIndex": [int, [2, 3]],
            "lowSpeedPreconditioner": [bool, False],
            "wallDistCutoff": [float, 1e20],
            "infChangeCorrection": [bool, True],
            "infChangeCorrectionTol": [float, 1e-12],
            "infChangeCorrectionType": [str, ["offset", "rotate"]],
            "cavitationNumber": [float, 1.4],
            "cpMinRho": [float, 100.0],
            # Common Parameters
            "nCycles": [int, 2000],
            "timeLimit": [float, -1.0],
            "nCyclesCoarse": [int, 500],
            "nSubiterTurb": [int, 3],
            "nSubiter": [int, 1],
            "CFL": [float, 1.7],
            "CFLCoarse": [float, 1.0],
            "MGCycle": [str, "sg"],
            "MGStartLevel": [int, -1],
            "resAveraging": [str, ["alternate", "never", "always"]],
            "smoothParameter": [float, 1.5],
            "CFLLimit": [float, 1.5],
            "useBlockettes": [bool, True],
            "useLinResMonitor": [bool, False],
            "useDissContinuation": [bool, False],
            "dissContMagnitude": [float, 1.0],
            "dissContMidpoint": [float, 3.0],
            "dissContSharpness": [float, 3.0],
            # Overset Parameters
            "nearWallDist": [float, 0.1],
            "backgroundVolScale": [float, 1.0],
            "oversetProjTol": [float, 1e-12],
            "overlapFactor": [float, 0.9],
            "oversetLoadBalance": [bool, True],
            "debugZipper": [bool, False],
            "zipperSurfaceFamily": [(str, type(None)), None],
            "cutCallback": [(types.FunctionType, type(None)), None],
            "explicitSurfaceCallback": [(types.FunctionType, type(None)), None],
            "oversetUpdateMode": [str, ["frozen", "fast", "full"]],
            "nRefine": [int, 10],
            "nFloodIter": [int, -1],
            "useZipperMesh": [bool, True],
            "useOversetWallScaling": [bool, False],
            "selfZipCutoff": [float, 120.0],
            "oversetPriority": [dict, {}],
            "oversetDebugPrint": [bool, False],
            # Unsteady Parameters
            "timeIntegrationScheme": [str, ["BDF", "explicit RK", "implicit RK"]],
            "timeAccuracy": [int, [2, 1, 3]],
            "nTimeStepsCoarse": [int, 48],
            "nTimeStepsFine": [int, 400],
            "deltaT": [float, 0.010],
            "useALE": [bool, True],
            "useGridMotion": [bool, False],
            "coupledSolution": [bool, False],
            # Time Spectral Parameters
            "timeIntervals": [int, 1],
            "alphaMode": [bool, False],
            "betaMode": [bool, False],
            "machMode": [bool, False],
            "pMode": [bool, False],
            "qMode": [bool, False],
            "rMode": [bool, False],
            "altitudeMode": [bool, False],
            "windAxis": [bool, False],
            "alphaFollowing": [bool, True],
            "TSStability": [bool, False],
            "useTSInterpolatedGridVelocity": [bool, False],
            "useExternalDynamicMesh": [bool, False],
            # Convergence Parameters
            "L2Convergence": [float, 1e-8],
            "L2ConvergenceRel": [float, 1e-16],
            "L2ConvergenceCoarse": [float, 1e-2],
            "maxL2DeviationFactor": [float, 1.0],
            # Newton-Krylov Parameters
            "useNKSolver": [bool, False],
            "NKSwitchTol": [float, 1e-5],
            "NKSubspaceSize": [int, 60],
            "NKLinearSolveTol": [float, 0.3],
            "NKUseEW": [bool, True],
            "NKADPC": [bool, False],
            "NKViscPC": [bool, False],
            "NKGlobalPreconditioner": [str, ["additive Schwarz", "multigrid"]],
            "NKASMOverlap": [int, 1],
            "NKPCILUFill": [int, 2],
            "NKJacobianLag": [int, 20],
            "applyPCSubspaceSize": [int, 10],
            "NKInnerPreconIts": [int, 1],
            "NKOuterPreconIts": [int, 1],
            "NKAMGLevels": [int, 2],
            "NKAMGNSmooth": [int, 1],
            "NKLS": [str, ["cubic", "none", "non-monotone"]],
            "NKFixedStep": [float, 0.25],
            "RKReset": [bool, False],
            "nRKReset": [int, 5],
            # Approximate Newton-Krylov Parameters
            "useANKSolver": [bool, True],
            "ANKUseTurbDADI": [bool, True],
            "ANKUseApproxSA": [bool, False],
            "ANKSwitchTol": [float, 1e3],
            "ANKSubspaceSize": [int, -1],
            "ANKMaxIter": [int, 40],
            "ANKLinearSolveTol": [float, 0.05],
            "ANKLinearSolveBuffer": [float, 0.01],
            "ANKLinResMax": [float, 0.1],
            "ANKGlobalPreconditioner": [str, ["additive Schwarz", "multigrid"]],
            "ANKASMOverlap": [int, 1],
            "ANKPCILUFill": [int, 2],
            "ANKJacobianLag": [int, 10],
            "ANKInnerPreconIts": [int, 1],
            "ANKOuterPreconIts": [int, 1],
            "ANKAMGLevels": [int, 2],
            "ANKAMGNSmooth": [int, 1],
            "ANKCFL0": [float, 5.0],
            "ANKCFLMin": [float, 1.0],
            "ANKCFLLimit": [float, 1e5],
            "ANKCFLFactor": [float, 10.0],
            "ANKCFLExponent": [float, 0.5],
            "ANKCFLCutback": [float, 0.5],
            "ANKCFLReset": [bool, True],
            "ANKStepFactor": [float, 1.0],
            "ANKStepMin": [float, 0.01],
            "ANKConstCFLStep": [float, 0.4],
            "ANKPhysicalLSTol": [float, 0.2],
            "ANKPhysicalLSTolTurb": [float, 0.99],
            "ANKUnsteadyLSTol": [float, 1.0],
            "ANKSecondOrdSwitchTol": [float, 1e-16],
            "ANKCoupledSwitchTol": [float, 1e-16],
            "ANKTurbCFLScale": [float, 1.0],
            "ANKUseFullVisc": [bool, True],
            "ANKPCUpdateTol": [float, 0.5],
            "ANKPCUpdateCutoff": [float, 1e-16],
            "ANKPCUpdateTolAfterCutoff": [float, 1e-4],
            "ANKADPC": [bool, False],
            "ANKNSubiterTurb": [int, 1],
            "ANKTurbKSPDebug": [bool, False],
            "ANKUseMatrixFree": [bool, True],
            "ANKCharTimeStepType": [str, ["None", "Turkel", "VLR"]],
            # Load Balance/partitioning parameters
            "blockSplitting": [bool, True],
            "loadImbalance": [float, 0.1],
            "loadBalanceIter": [int, 10],
            "partitionOnly": [bool, False],
            "partitionLikeNProc": [int, -1],
            # Misc Parameters
            "numberSolutions": [bool, True],
            "writeSolutionDigits": [int, 3],
            "printIterations": [bool, True],
            "printTiming": [bool, True],
            "printIntro": [bool, True],
            "printAllOptions": [bool, True],
            "setMonitor": [bool, True],
            "printWarnings": [bool, True],
            "printNegativeVolumes": [bool, False],
            "printBadlySkewedCells": [bool, False],
            "monitorVariables": [list, ["cpu", "resrho", "resturb", "cl", "cd"]],
            "surfaceVariables": [list, ["cp", "vx", "vy", "vz", "mach"]],
            "volumeVariables": [list, ["resrho"]],
            "storeConvHist": [bool, True],
            # Multidisciplinary Coupling Parameters
            "forcesAsTractions": [bool, True],
            # Adjoint Parameters
            "adjointL2Convergence": [float, 1e-6],
            "adjointL2ConvergenceRel": [float, 1e-16],
            "adjointL2ConvergenceAbs": [float, 1e-16],
            "adjointDivTol": [float, 1e5],
            "adjointMaxL2DeviationFactor": [float, 1.0],
            "approxPC": [bool, True],
            "ADPC": [bool, False],
            "viscPC": [bool, False],
            "useDiagTSPC": [bool, True],
            "restartAdjoint": [bool, True],
            "adjointSolver": [str, ["GMRES", "TFQMR", "Richardson", "BCGS", "IBCGS"]],
            "adjointMaxIter": [int, 500],
            "adjointSubspaceSize": [int, 100],
            "GMRESOrthogonalizationType": [
                str,
                ["modified Gram-Schmidt", "CGS never refine", "CGS refine if needed", "CGS always refine"],
            ],
            "adjointMonitorStep": [int, 10],
            "dissipationLumpingParameter": [float, 6.0],
            "preconditionerSide": [str, ["right", "left"]],
            "matrixOrdering": [
                str,
                ["RCM", "natural", "nested dissection", "one way dissection", "quotient minimum degree"],
            ],
            "globalPreconditioner": [str, ["additive Schwarz", "multigrid"]],
            "localPreconditioner": [str, ["ILU"]],
            "ILUFill": [int, 2],
            "ASMOverlap": [int, 1],
            "innerPreconIts": [int, 1],
            "outerPreconIts": [int, 3],
            "adjointAMGLevels": [int, 2],
            "adjointAMGNSmooth": [int, 1],
            "applyAdjointPCSubspaceSize": [int, 20],
            "frozenTurbulence": [bool, False],
            "useMatrixFreedrdw": [bool, True],
            "skipAfterFailedAdjoint": [bool, True],
            # ADjoint debugger
            "firstRun": [bool, True],
            "verifyState": [bool, True],
            "verifySpatial": [bool, True],
            "verifyExtra": [bool, True],
            # Function parmeters
            "sepSensorOffset": [float, 0.0],
            "sepSensorSharpness": [float, 10.0],
            "cavSensorOffset": [float, 0.0],
            "cavSensorSharpness": [float, 10.0],
            "cavExponent": [int, 0],
            "computeCavitation": [bool, False],
        }

        return defOpts

    def _getImmutableOptions(self):
        """We define the list of options that *cannot* be changed after the
        object is created. ADflow will raise an error if a user tries to
        change these. The strings for these options are placed in a set"""

        return (
            "gridfile",
            "equationtype",
            "equationmode",
            "flowtype",
            "useapproxwalldistance",
            "liftindex",
            "mgcycle",
            "mgstartlevel",
            "timeintegrationscheme",
            "timeaccuracy",
            "useale",
            "timeintervals",
            "blocksplitting",
            "loadimbalance",
            "loadbalanceiter",
            "partitiononly",
            "meshsurfacefamily",
            "designsurfacefamily",
            "closedsurfacefamilies",
            "zippersurfacefamily",
            "cutcallback",
            "explicitsurfacecallback",
            "useskewnesscheck",
        )

    def _getOptionMap(self):
        """The ADflow option map and module mapping"""

        moduleMap = {
            "io": self.adflow.inputio,
            "discr": self.adflow.inputdiscretization,
            "iter": self.adflow.inputiteration,
            "physics": self.adflow.inputphysics,
            "stab": self.adflow.inputtsstabderiv,
            "nk": self.adflow.nksolver,
            "ank": self.adflow.anksolver,
            "amg": self.adflow.amg,
            "adjoint": self.adflow.inputadjoint,
            "cost": self.adflow.inputcostfunctions,
            "unsteady": self.adflow.inputunsteady,
            "motion": self.adflow.inputmotion,
            "parallel": self.adflow.inputparallel,
            "ts": self.adflow.inputtimespectral,
            "overset": self.adflow.inputoverset,
            "monitor": self.adflow.monitor,
        }

        # In the option map, we first list the "module" defined in
        # module map, and "variable" the variable to set in that module.

        optionMap = {
            # Common Parameters
            "gridfile": ["io", "gridfile"],
            "storerindlayer": ["io", "storerindlayer"],
            "nsavevolume": ["io", "nsavevolume"],
            "nsavesurface": ["iter", "nsavesurface"],
            "viscoussurfacevelocities": ["io", "viscoussurfacevelocities"],
            "solutionprecision": {
                "single": self.adflow.constants.precisionsingle,
                "double": self.adflow.constants.precisiondouble,
                "location": ["io", "precisionsol"],
            },
            "gridprecision": {
                "single": self.adflow.constants.precisionsingle,
                "double": self.adflow.constants.precisiondouble,
                "location": ["io", "precisiongrid"],
            },
            "solutionprecisionsurface": {
                "single": self.adflow.constants.precisionsingle,
                "double": self.adflow.constants.precisiondouble,
                "location": ["io", "precisionsurfsol"],
            },
            "gridprecisionsurface": {
                "single": self.adflow.constants.precisionsingle,
                "double": self.adflow.constants.precisiondouble,
                "location": ["io", "precisionsurfgrid"],
            },
            # Physics Parameters
            "discretization": {
                "central plus scalar dissipation": self.adflow.constants.dissscalar,
                "central plus matrix dissipation": self.adflow.constants.dissmatrix,
                "central plus cusp dissipation": self.adflow.constants.disscusp,
                "upwind": self.adflow.constants.upwind,
                "location": ["discr", "spacediscr"],
            },
            "coarsediscretization": {
                "central plus scalar dissipation": self.adflow.constants.dissscalar,
                "central plus matrix dissipation": self.adflow.constants.dissmatrix,
                "central plus cusp dissipation": self.adflow.constants.disscusp,
                "upwind": self.adflow.constants.upwind,
                "location": ["discr", "spacediscrcoarse"],
            },
            "limiter": {
                "van albada": self.adflow.constants.vanalbeda,
                "minmod": self.adflow.constants.minmod,
                "no limiter": self.adflow.constants.nolimiter,
                "location": ["discr", "limiter"],
            },
            "smoother": {
                "runge-kutta": self.adflow.constants.rungekutta,
                "lu sgs": self.adflow.constants.nllusgs,
                "lu sgs line": self.adflow.constants.nllusgsline,
                "dadi": self.adflow.constants.dadi,
                "location": ["iter", "smoother"],
            },
            "equationtype": {
                "euler": self.adflow.constants.eulerequations,
                "laminar ns": self.adflow.constants.nsequations,
                "rans": self.adflow.constants.ransequations,
                "location": ["physics", "equations"],
            },
            "equationmode": {
                "steady": self.adflow.constants.steady,
                "unsteady": self.adflow.constants.unsteady,
                "time spectral": self.adflow.constants.timespectral,
                "location": ["physics", "equationmode"],
            },
            "flowtype": {
                "internal": self.adflow.constants.internalflow,
                "external": self.adflow.constants.externalflow,
                "location": ["physics", "flowtype"],
            },
            "turbulencemodel": {
                "sa": self.adflow.constants.spalartallmaras,
                "sa-edwards": self.adflow.constants.spalartallmarasedwards,
                "k-omega wilcox": self.adflow.constants.komegawilcox,
                "k-omega modified": self.adflow.constants.komegamodified,
                "k-tau": self.adflow.constants.ktau,
                "menter sst": self.adflow.constants.mentersst,
                "v2f": self.adflow.constants.v2f,
                "location": ["physics", "turbmodel"],
            },
            "turbulenceorder": {"first order": 1, "second order": 2, "location": ["discr", "orderturb"]},
            "turbresscale": ["iter", "turbresscale"],
            "meshmaxskewness": ["iter", "meshmaxskewness"],
            "useskewnesscheck": ["iter", "useskewnesscheck"],
            "turbulenceproduction": {
                "strain": self.adflow.constants.strain,
                "vorticity": self.adflow.constants.vorticity,
                "kato-launder": self.adflow.constants.katolaunder,
                "location": ["physics", "turbprod"],
            },
            "useqcr": ["physics", "useqcr"],
            "userotationsa": ["physics", "userotationsa"],
            "useft2sa": ["physics", "useft2sa"],
            "eddyvisinfratio": ["physics", "eddyvisinfratio"],
            "usewallfunctions": ["physics", "wallfunctions"],
            "walldistcutoff": ["physics", "walldistcutoff"],
            "useapproxwalldistance": ["discr", "useapproxwalldistance"],
            "eulerwalltreatment": {
                "linear pressure extrapolation": self.adflow.constants.linextrapolpressure,
                "constant pressure extrapolation": self.adflow.constants.constantpressure,
                "quadratic pressure extrapolation": self.adflow.constants.quadextrapolpressure,
                "normal momentum": self.adflow.constants.normalmomentum,
                "location": ["discr", "eulerwallbctreatment"],
            },
            "viscwalltreatment": {
                "linear pressure extrapolation": self.adflow.constants.linextrapolpressure,
                "constant pressure extrapolation": self.adflow.constants.constantpressure,
                "location": ["discr", "viscwallbctreatment"],
            },
            "dissipationscalingexponent": ["discr", "adis"],
            "acousticscalefactor": ["discr", "acousticscalefactor"],
            "vis4": ["discr", "vis4"],
            "vis2": ["discr", "vis2"],
            "vis2coarse": ["discr", "vis2coarse"],
            "restrictionrelaxation": ["iter", "fcoll"],
            "forcesastractions": ["physics", "forcesastractions"],
            "lowspeedpreconditioner": ["discr", "lowspeedpreconditioner"],
            "cavitationnumber": ["physics", "cavitationnumber"],
            "cpminrho": ["physics", "cpmin_rho"],
            # Common Parameters
            "ncycles": ["iter", "ncycles"],
            "timelimit": ["iter", "timelimit"],
            "ncyclescoarse": ["iter", "ncyclescoarse"],
            "nsubiterturb": ["iter", "nsubiterturb"],
            "nsubiter": ["iter", "nsubiterations"],
            "cfl": ["iter", "cfl"],
            "cflcoarse": ["iter", "cflcoarse"],
            "mgcycle": ["iter", "mgdescription"],
            "mgstartlevel": ["iter", "mgstartlevel"],
            "resaveraging": {
                "never": self.adflow.constants.noresaveraging,
                "always": self.adflow.constants.alwaysresaveraging,
                "alternate": self.adflow.constants.alternateresaveraging,
                "location": ["iter", "resaveraging"],
            },
            "smoothparameter": ["iter", "smoop"],
            "cfllimit": ["iter", "cfllimit"],
            "useblockettes": ["discr", "useblockettes"],
            "uselinresmonitor": ["iter", "uselinresmonitor"],
            "usedisscontinuation": ["iter", "usedisscontinuation"],
            "disscontmagnitude": ["iter", "disscontmagnitude"],
            "disscontmidpoint": ["iter", "disscontmidpoint"],
            "disscontsharpness": ["iter", "disscontsharpness"],
            # Overset Parameters
            "nearwalldist": ["overset", "nearwalldist"],
            "backgroundvolscale": ["overset", "backgroundvolscale"],
            "oversetprojtol": ["overset", "oversetprojtol"],
            "overlapfactor": ["overset", "overlapfactor"],
            "oversetloadbalance": ["overset", "useoversetloadbalance"],
            "debugzipper": ["overset", "debugzipper"],
            "oversetupdatemode": {
                "frozen": self.adflow.constants.updatefrozen,
                "fast": self.adflow.constants.updatefast,
                "full": self.adflow.constants.updatefull,
                "location": ["overset", "oversetupdatemode"],
            },
            "nrefine": ["overset", "nrefine"],
            "nflooditer": ["overset", "nflooditer"],
            "usezippermesh": ["overset", "usezippermesh"],
            "useoversetwallscaling": ["overset", "useoversetwallscaling"],
            "selfzipcutoff": ["overset", "selfzipcutoff"],
            "oversetdebugprint": ["overset", "oversetdebugprint"],
            # Unsteady Params
            "timeintegrationscheme": {
                "bdf": self.adflow.constants.bdf,
                "explicit rk": self.adflow.constants.explicitrk,
                "implicit rk": self.adflow.constants.implicitrk,
                "location": ["unsteady", "timeintegrationscheme"],
            },
            "timeaccuracy": ["unsteady", "timeaccuracy"],
            "ntimestepscoarse": ["unsteady", "ntimestepscoarse"],
            "ntimestepsfine": ["unsteady", "ntimestepsfine"],
            "deltat": ["unsteady", "deltat"],
            "useale": ["unsteady", "useale"],
            # Grid motion Params
            "usegridmotion": ["motion", "gridmotionspecified"],
            # Time Spectral Parameters
            "timeintervals": ["ts", "ntimeintervalsspectral"],
            "alphamode": ["stab", "tsalphamode"],
            "betamode": ["stab", "tsbetamode"],
            "machmode": ["stab", "tsmachmode"],
            "pmode": ["stab", "tspmode"],
            "qmode": ["stab", "tsqmode"],
            "rmode": ["stab", "tsrmode"],
            "altitudemode": ["stab", "tsaltitudemode"],
            "windaxis": ["stab", "usewindaxis"],
            "alphafollowing": ["stab", "tsalphafollowing"],
            "tsstability": ["stab", "tsstability"],
            "usetsinterpolatedgridvelocity": ["ts", "usetsinterpolatedgridvelocity"],
            # Convergence Parameters
            "l2convergence": ["iter", "l2conv"],
            "l2convergencerel": ["iter", "l2convrel"],
            "l2convergencecoarse": ["iter", "l2convcoarse"],
            "maxl2deviationfactor": ["iter", "maxl2deviationfactor"],
            # Newton-Krylov Parameters
            "usenksolver": ["nk", "usenksolver"],
            "nkuseew": ["nk", "nk_useew"],
            "nkswitchtol": ["nk", "nk_switchtol"],
            "nksubspacesize": ["nk", "nk_subspace"],
            "nklinearsolvetol": ["nk", "nk_rtolinit"],
            "nkglobalpreconditioner": {
                "additive schwarz": "asm",
                "multigrid": "mg",
                "location": ["nk", "nk_precondtype"],
            },
            "nkasmoverlap": ["nk", "nk_asmoverlap"],
            "nkpcilufill": ["nk", "nk_ilufill"],
            "nkjacobianlag": ["nk", "nk_jacobianlag"],
            "nkadpc": ["nk", "nk_adpc"],
            "nkviscpc": ["nk", "nk_viscpc"],
            "applypcsubspacesize": ["nk", "applypcsubspacesize"],
            "nkinnerpreconits": ["nk", "nk_innerpreconits"],
            "nkouterpreconits": ["nk", "nk_outerpreconits"],
            "nkamglevels": ["nk", "nk_amglevels"],
            "nkamgnsmooth": ["nk", "nk_amgnsmooth"],
            "nkls": {
                "none": self.adflow.constants.nolinesearch,
                "cubic": self.adflow.constants.cubiclinesearch,
                "non-monotone": self.adflow.constants.nonmonotonelinesearch,
                "location": ["nk", "nk_ls"],
            },
            "nkfixedstep": ["nk", "nk_fixedstep"],
            "rkreset": ["iter", "rkreset"],
            "nrkreset": ["iter", "miniternum"],
            # Approximate Newton-Krylov Parameters
            "useanksolver": ["ank", "useanksolver"],
            "ankuseturbdadi": ["ank", "ank_useturbdadi"],
            "ankuseapproxsa": ["ank", "ank_useapproxsa"],
            "ankswitchtol": ["ank", "ank_switchtol"],
            "anksubspacesize": ["ank", "ank_subspace"],
            "ankmaxiter": ["ank", "ank_maxiter"],
            "anklinearsolvetol": ["ank", "ank_rtol"],
            "anklinearsolvebuffer": ["ank", "ank_atol_buffer"],
            "anklinresmax": ["ank", "ank_linresmax"],
            "ankglobalpreconditioner": {
                "additive schwarz": "asm",
                "multigrid": "mg",
                "location": ["ank", "ank_precondtype"],
            },
            "ankasmoverlap": ["ank", "ank_asmoverlap"],
            "ankpcilufill": ["ank", "ank_ilufill"],
            "ankjacobianlag": ["ank", "ank_jacobianlag"],
            "ankinnerpreconits": ["ank", "ank_innerpreconits"],
            "ankouterpreconits": ["ank", "ank_outerpreconits"],
            "ankamglevels": ["ank", "ank_amglevels"],
            "ankamgnsmooth": ["ank", "ank_amgnsmooth"],
            "ankcfl0": ["ank", "ank_cfl0"],
            "ankcflmin": ["ank", "ank_cflmin0"],
            "ankcfllimit": ["ank", "ank_cfllimit"],
            "ankcflfactor": ["ank", "ank_cflfactor"],
            "ankcflexponent": ["ank", "ank_cflexponent"],
            "ankcflcutback": ["ank", "ank_cflcutback"],
            "ankcflreset": ["ank", "ank_cflreset"],
            "ankstepfactor": ["ank", "ank_stepfactor"],
            "ankstepmin": ["ank", "ank_stepmin"],
            "ankconstcflstep": ["ank", "ank_constcflstep"],
            "ankphysicallstol": ["ank", "ank_physlstol"],
            "ankphysicallstolturb": ["ank", "ank_physlstolturb"],
            "ankunsteadylstol": ["ank", "ank_unstdylstol"],
            "anksecondordswitchtol": ["ank", "ank_secondordswitchtol"],
            "ankcoupledswitchtol": ["ank", "ank_coupledswitchtol"],
            "ankturbcflscale": ["ank", "ank_turbcflscale"],
            "ankusefullvisc": ["ank", "ank_usefullvisc"],
            "ankpcupdatetol": ["ank", "ank_pcupdatetol"],
            "ankpcupdatecutoff": ["ank", "ank_pcupdatecutoff"],
            "ankpcupdatetolaftercutoff": ["ank", "ank_pcupdatetol2"],
            "ankadpc": ["ank", "ank_adpc"],
            "anknsubiterturb": ["ank", "ank_nsubiterturb"],
            "ankturbkspdebug": ["ank", "ank_turbdebug"],
            "ankusematrixfree": ["ank", "ank_usematrixfree"],
            "ankchartimesteptype": {
                "none": "None",
                "turkel": "Turkel",
                "vlr": "VLR",
                "location": ["ank", "ank_chartimesteptype"],
            },
            # Load Balance Parameters
            "blocksplitting": ["parallel", "splitblocks"],
            "loadimbalance": ["parallel", "loadimbalance"],
            "loadbalanceiter": ["parallel", "loadbalanceiter"],
            "partitionlikenproc": ["parallel", "partitionlikenproc"],
            # Misc Parameters
            "printiterations": ["iter", "printiterations"],
            "printwarnings": ["iter", "printwarnings"],
            "printnegativevolumes": ["iter", "printnegativevolumes"],
            "printbadlyskewedcells": ["iter", "printbadlyskewedcells"],
            "printtiming": ["adjoint", "printtiming"],
            "setmonitor": ["adjoint", "setmonitor"],
            "storeconvhist": ["io", "storeconvinneriter"],
            # Adjoint Params
            "adjointl2convergence": ["adjoint", "adjreltol"],
            "adjointl2convergencerel": ["adjoint", "adjreltolrel"],
            "adjointl2convergenceabs": ["adjoint", "adjabstol"],
            "adjointdivtol": ["adjoint", "adjdivtol"],
            "adjointmaxl2deviationfactor": ["adjoint", "adjmaxl2dev"],
            "approxpc": ["adjoint", "approxpc"],
            "adpc": ["adjoint", "adpc"],
            "viscpc": ["adjoint", "viscpc"],
            "frozenturbulence": ["adjoint", "frozenturbulence"],
            "usediagtspc": ["adjoint", "usediagtspc"],
            "adjointsolver": {
                "gmres": "gmres",
                "tfqmr": "tfqmr",
                "richardson": "richardson",
                "bcgs": "bcgs",
                "ibcgs": "ibcgs",
                "location": ["adjoint", "adjointsolvertype"],
            },
            "gmresorthogonalizationtype": {
                "modified gram-schmidt": "modified_gram_schmidt",
                "cgs never refine": "cgs_never_refine",
                "cgs refine if needed": "cgs_refine_if_needed",
                "cgs always refine": "cgs_always_refine",
                "location": ["adjoint", "gmresorthogtype"],
            },
            "adjointmaxiter": ["adjoint", "adjmaxiter"],
            "adjointsubspacesize": ["adjoint", "adjrestart"],
            "adjointmonitorstep": ["adjoint", "adjmonstep"],
            "dissipationlumpingparameter": ["discr", "sigma"],
            "preconditionerside": {"left": "left", "right": "right", "location": ["adjoint", "adjointpcside"]},
            "matrixordering": {
                "natural": "natural",
                "rcm": "rcm",
                "nested dissection": "nd",
                "one way dissection": "1wd",
                "quotient minimum degree": "qmd",
                "location": ["adjoint", "matrixordering"],
            },
            "globalpreconditioner": {
                "additive schwarz": "asm",
                "multigrid": "mg",
                "location": ["adjoint", "precondtype"],
            },
            "localpreconditioner": {"ilu": "ilu", "location": ["adjoint", "localpctype"]},
            "ilufill": ["adjoint", "filllevel"],
            "applyadjointpcsubspacesize": ["adjoint", "applyadjointpcsubspacesize"],
            "asmoverlap": ["adjoint", "overlap"],
            "innerpreconits": ["adjoint", "innerpreconits"],
            "outerpreconits": ["adjoint", "outerpreconits"],
            "adjointamglevels": ["adjoint", "adjamglevels"],
            "adjointamgnsmooth": ["adjoint", "adjamgnsmooth"],
            "firstrun": ["adjoint", "firstrun"],
            "verifystate": ["adjoint", "verifystate"],
            "verifyspatial": ["adjoint", "verifyspatial"],
            "verifyextra": ["adjoint", "verifyextra"],
            "usematrixfreedrdw": ["adjoint", "usematrixfreedrdw"],
            # Parameters for functions
            "sepsensoroffset": ["cost", "sepsensoroffset"],
            "sepsensorsharpness": ["cost", "sepsensorsharpness"],
            "cavsensoroffset": ["cost", "cavsensoroffset"],
            "cavsensorsharpness": ["cost", "cavsensorsharpness"],
            "cavexponent": ["cost", "cavexponent"],
            "computecavitation": ["cost", "computecavitation"],
            "writesolutioneachiter": ["monitor", "writesoleachiter"],
            "writesurfacesolution": ["monitor", "writesurface"],
            "writevolumesolution": ["monitor", "writevolume"],
        }
        return optionMap, moduleMap

    def _getSpecialOptionLists(self):
        """
        Lists of special options
        """
        # pythonOptions do not get set in the Fortran code.
        # They are used strictly in Python.

        pythonOptions = {
            "numbersolutions",
            "writesolutiondigits",
            "writetecplotsurfacesolution",
            "coupledsolution",
            "partitiononly",
            "liftindex",
            "meshsurfacefamily",
            "designsurfacefamily",
            "closedsurfacefamilies",
            "zippersurfacefamily",
            "outputsurfacefamily",
            "cutcallback",
            "explicitsurfacecallback",
            "infchangecorrection",
            "infchangecorrectiontol",
            "infchangecorrectiontype",
            "restartadjoint",
            "skipafterfailedadjoint",
            "useexternaldynamicmesh",
            "printalloptions",
            "printintro",
        }

        # Deprecated options that may be in old scripts and should not be used.

        deprecatedOptions = {
            "finitedifferencepc": "Use the ADPC option.",
            "writesolution": "Use writeSurfaceSolution and writeVolumeSolution options instead.",
            "autosolveretry": "This feature is not implemented.",
            "autoadjointretry": "This feature is not implemented.",
            "nkcfl0": "The NK solver does not use a CFL value anymore. "
            + "The CFL is set to infinity and the true Newton method is used.",
        }

        specialOptions = {
            "surfacevariables",
            "volumevariables",
            "monitorvariables",
            "outputdirectory",
            "isovariables",
            "isosurface",
            "turbresscale",
            "restartfile",
            "oversetpriority",
        }

        return pythonOptions, deprecatedOptions, specialOptions

    def _getObjectivesAndDVs(self):
        iDV = OrderedDict()
        iDV["alpha"] = self.adflow.adjointvars.ialpha
        iDV["beta"] = self.adflow.adjointvars.ibeta
        iDV["mach"] = self.adflow.adjointvars.imach
        iDV["machgrid"] = self.adflow.adjointvars.imachgrid
        iDV["p"] = self.adflow.adjointvars.ipressure
        iDV["rho"] = self.adflow.adjointvars.idensity
        iDV["t"] = self.adflow.adjointvars.itemperature
        iDV["rotx"] = self.adflow.adjointvars.irotx
        iDV["roty"] = self.adflow.adjointvars.iroty
        iDV["rotz"] = self.adflow.adjointvars.irotz
        iDV["rotcenx"] = self.adflow.adjointvars.irotcenx
        iDV["rotceny"] = self.adflow.adjointvars.irotceny
        iDV["rotcenz"] = self.adflow.adjointvars.irotcenz
        iDV["xref"] = self.adflow.adjointvars.ipointrefx
        iDV["yref"] = self.adflow.adjointvars.ipointrefy
        iDV["zref"] = self.adflow.adjointvars.ipointrefz

        # Convert to python indexing
        for key in iDV:
            iDV[key] = iDV[key] - 1

        # Extra DVs for the boundary condition variables
        BCDV = ["pressure", "pressurestagnation", "temperaturestagnation", "thrust", "heat"]

        # This is ADflow's internal mapping for cost functions
        adflowCostFunctions = {
            "lift": self.adflow.constants.costfunclift,
            "drag": self.adflow.constants.costfuncdrag,
            "cl": self.adflow.constants.costfuncliftcoef,
            "cd": self.adflow.constants.costfuncdragcoef,
            "clp": self.adflow.constants.costfuncliftcoefpressure,
            "clv": self.adflow.constants.costfuncliftcoefviscous,
            "clm": self.adflow.constants.costfuncliftcoefmomentum,
            "cdp": self.adflow.constants.costfuncdragcoefpressure,
            "cdv": self.adflow.constants.costfuncdragcoefviscous,
            "cdm": self.adflow.constants.costfuncdragcoefmomentum,
            "fx": self.adflow.constants.costfuncforcex,
            "fy": self.adflow.constants.costfuncforcey,
            "fz": self.adflow.constants.costfuncforcez,
            "cfx": self.adflow.constants.costfuncforcexcoef,
            "cfxp": self.adflow.constants.costfuncforcexcoefpressure,
            "cfxv": self.adflow.constants.costfuncforcexcoefviscous,
            "cfxm": self.adflow.constants.costfuncforcexcoefmomentum,
            "cfy": self.adflow.constants.costfuncforceycoef,
            "cfyp": self.adflow.constants.costfuncforceycoefpressure,
            "cfyv": self.adflow.constants.costfuncforceycoefviscous,
            "cfym": self.adflow.constants.costfuncforceycoefmomentum,
            "cfz": self.adflow.constants.costfuncforcezcoef,
            "cfzp": self.adflow.constants.costfuncforcezcoefpressure,
            "cfzv": self.adflow.constants.costfuncforcezcoefviscous,
            "cfzm": self.adflow.constants.costfuncforcezcoefmomentum,
            "mx": self.adflow.constants.costfuncmomx,
            "my": self.adflow.constants.costfuncmomy,
            "mz": self.adflow.constants.costfuncmomz,
            "cmx": self.adflow.constants.costfuncmomxcoef,
            "cmy": self.adflow.constants.costfuncmomycoef,
            "cmz": self.adflow.constants.costfuncmomzcoef,
            "cm0": self.adflow.constants.costfunccm0,
            "cmzalpha": self.adflow.constants.costfunccmzalpha,
            "cmzalphadot": self.adflow.constants.costfunccmzalphadot,
            "cl0": self.adflow.constants.costfunccl0,
            "clalpha": self.adflow.constants.costfuncclalpha,
            "clalphadot": self.adflow.constants.costfuncclalphadot,
            "cfy0": self.adflow.constants.costfunccfy0,
            "cfyalpha": self.adflow.constants.costfunccfyalpha,
            "cfyalphadot": self.adflow.constants.costfunccfyalphadot,
            "cd0": self.adflow.constants.costfunccd0,
            "cdalpha": self.adflow.constants.costfunccdalpha,
            "cdalphadot": self.adflow.constants.costfunccdalphadot,
            "cmzq": self.adflow.constants.costfunccmzq,
            "cmzqdot": self.adflow.constants.costfunccmzqdot,
            "clq": self.adflow.constants.costfuncclq,
            "clqdot": self.adflow.constants.costfuncclqdot,
            "cbend": self.adflow.constants.costfuncbendingcoef,
            "sepsensor": self.adflow.constants.costfuncsepsensor,
            "sepsensoravgx": self.adflow.constants.costfuncsepsensoravgx,
            "sepsensoravgy": self.adflow.constants.costfuncsepsensoravgy,
            "sepsensoravgz": self.adflow.constants.costfuncsepsensoravgz,
            "cavitation": self.adflow.constants.costfunccavitation,
            "cpmin": self.adflow.constants.costfunccpmin,
            "mdot": self.adflow.constants.costfuncmdot,
            "mavgptot": self.adflow.constants.costfuncmavgptot,
            "aavgptot": self.adflow.constants.costfuncaavgptot,
            "aavgps": self.adflow.constants.costfuncaavgps,
            "mavgttot": self.adflow.constants.costfuncmavgttot,
            "mavgps": self.adflow.constants.costfuncmavgps,
            "mavgmn": self.adflow.constants.costfuncmavgmn,
            "area": self.adflow.constants.costfuncarea,
            "axismoment": self.adflow.constants.costfuncaxismoment,
            "flowpower": self.adflow.constants.costfuncflowpower,
            "forcexpressure": self.adflow.constants.costfuncforcexpressure,
            "forceypressure": self.adflow.constants.costfuncforceypressure,
            "forcezpressure": self.adflow.constants.costfuncforcezpressure,
            "forcexviscous": self.adflow.constants.costfuncforcexviscous,
            "forceyviscous": self.adflow.constants.costfuncforceyviscous,
            "forcezviscous": self.adflow.constants.costfuncforcezviscous,
            "forcexmomentum": self.adflow.constants.costfuncforcexmomentum,
            "forceymomentum": self.adflow.constants.costfuncforceymomentum,
            "forcezmomentum": self.adflow.constants.costfuncforcezmomentum,
            "dragpressure": self.adflow.constants.costfuncdragpressure,
            "dragviscous": self.adflow.constants.costfuncdragviscous,
            "dragmomentum": self.adflow.constants.costfuncdragmomentum,
            "liftpressure": self.adflow.constants.costfuncliftpressure,
            "liftviscous": self.adflow.constants.costfuncliftviscous,
            "liftmomentum": self.adflow.constants.costfuncliftmomentum,
            "mavgvx": self.adflow.constants.costfuncmavgvx,
            "mavgvy": self.adflow.constants.costfuncmavgvy,
            "mavgvz": self.adflow.constants.costfuncmavgvz,
            "mavgvi": self.adflow.constants.costfuncmavgvi,
            "cperror2": self.adflow.constants.costfunccperror2,
            "cofxx": self.adflow.constants.costfunccoforcexx,
            "cofxy": self.adflow.constants.costfunccoforcexy,
            "cofxz": self.adflow.constants.costfunccoforcexz,
            "cofyx": self.adflow.constants.costfunccoforceyx,
            "cofyy": self.adflow.constants.costfunccoforceyy,
            "cofyz": self.adflow.constants.costfunccoforceyz,
            "cofzx": self.adflow.constants.costfunccoforcezx,
            "cofzy": self.adflow.constants.costfunccoforcezy,
            "cofzz": self.adflow.constants.costfunccoforcezz,
            "colx": self.adflow.constants.costfunccofliftx,
            "coly": self.adflow.constants.costfunccoflifty,
            "colz": self.adflow.constants.costfunccofliftz,
        }

        return iDV, BCDV, adflowCostFunctions

    def _updateTurbResScale(self):
        # If turbresscale is None it has not been set by the user in script
        # thus set default values depending on turbulence model;
        # else do nothing since it already contains value specified by user

        if self.getOption("turbresscale") is None:
            turbModel = self.getOption("turbulencemodel")
            if turbModel == "SA":
                self.setOption("turbresscale", 10000.0)
            elif turbModel == "Menter SST":
                self.setOption("turbresscale", [1e3, 1e-6])
            else:
                raise Error(
                    "Turbulence model %-35s does not have default values specified for turbresscale. Specify turbresscale manually or update the python interface"
                    % (turbModel)
                )

    def _setUnsteadyFileParameters(self):
        """
        This is a temporary function that sets the appropriate filenames
        when using the unsteady equations and the bdf time integration scheme.
        This function should drop out when unsteady solver is stepped through python
        TEMPORARY
        """
        # THIS DOES NOT CURRENTLY WORK DUE TO INTERNAL LOGIC. REFACTOR FORTRAN
        # Set parameters for outputing data
        # if self.getOption('writevolumesolution'):
        #    self.adflow.monitor.writevolume = True
        #    self.adflow.monitor.writegrid = True
        # else:
        #    self.adflow.monitor.writevolume = False
        #    self.adflow.monitor.writegrid = False

        # if self.getOption('writesurfacesolution'):
        #    self.adflow.monitor.writesurface = True
        # else:
        #    self.adflow.monitor.writesurface = False

        outputDir = self.getOption("outputDirectory")
        baseName = self.curAP.name

        # Join to get the actual filename root
        volFileName = os.path.join(outputDir, baseName + "_vol.cgns")
        surfFileName = os.path.join(outputDir, baseName + "_surf.cgns")
        sliceFileName = os.path.join(outputDir, baseName + "_slices")
        liftDistributionFileName = os.path.join(outputDir, baseName + "_lift")

        # Set fileName in adflow
        self.adflow.inputio.solfile = self._expandString(volFileName)

        # Set the grid file to the same name so the grids will be written
        # to the volume files
        self.adflow.inputio.newgridfile = self._expandString(volFileName)
        self.adflow.inputio.surfacesolfile = self._expandString(surfFileName)
        self.adflow.inputio.slicesolfile = self._expandString(sliceFileName)
        self.adflow.inputio.liftdistributionfile = self._expandString(liftDistributionFileName)

    def _setForcedFileNames(self):
        # Set the filenames that will be used if the user forces a
        # write during a solution.

        outputDir = self.getOption("outputDirectory")
        baseName = self.curAP.name
        base = os.path.join(outputDir, baseName)

        volFileName = base + "_forced_vol.cgns"
        surfFileName = base + "_forced_surf.cgns"
        liftFileName = base + "_forced_lift.dat"
        sliceFileName = base + "_forced_slices.dat"

        convSolFileBaseName = base + "_intermediate_sol"

        self.adflow.inputio.convsolfilebasename = self._expandString(convSolFileBaseName)
        self.adflow.inputio.forcedvolumefile = self._expandString(volFileName)
        self.adflow.inputio.forcedsurfacefile = self._expandString(surfFileName)
        self.adflow.inputio.forcedliftfile = self._expandString(liftFileName)
        self.adflow.inputio.forcedslicefile = self._expandString(sliceFileName)

    def _createZipperMesh(self):
        """Internal routine for generating the zipper mesh. This operation is
        postposted as long as possible and now it cannot wait any longer."""

        # Verify if we already have previous failures, such as negative volumes
        self.adflow.killsignals.routinefailed = self.comm.allreduce(
            bool(self.adflow.killsignals.routinefailed), op=MPI.LOR
        )

        # Avoid zipper if we have failures or if it already exists
        if self.zipperCreated or self.adflow.killsignals.routinefailed:
            return

        zipFam = self.getOption("zipperSurfaceFamily")

        if zipFam is None:
            # The user didn't tell us anything. So we will use all
            # walls. Remind the user what those are.
            zipperFamList = self.families[self.allWallsGroup]
            if self.myid == 0:
                ADFLOWWarning(
                    "'zipperSurfaceFamily' option was not given. Using all "
                    "wall boundary conditions for the zipper mesh."
                )
        else:
            if zipFam not in self.families:
                raise Error(
                    "Trying to create the zipper mesh, but '%s' is not a "
                    "family in the CGNS file or has not been added"
                    " as a combination of families" % zipFam
                )
            zipperFamList = self.families[zipFam]

        self.adflow.zippermesh.createzippermesh(zipperFamList)
        self.zipperCreated = True
        self.coords0 = self.getSurfaceCoordinates(self.allFamilies)

    def finalizeUserIntegrationSurfaces(self):
        if self.hasIntegrationSurfaces and (not self.userSurfaceIntegrationsFinalized):
            # We can also do the surface plane integrations here if necessary
            self.adflow.usersurfaceintegrations.interpolateintegrationsurfaces()
            self.userSurfaceIntegrationsFinalized = True

    def _expandString(self, s):
        """Expand a supplied string 's' to be of the constants.maxstring
        length so we can set them in fortran"""
        return s + " " * (256 - len(s))

    def _createFortranStringArray(self, strList):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takesa list of strings and returns the array"""

        arr = numpy.zeros((len(strList), self.adflow.constants.maxcgnsnamelen), dtype="str")
        arr[:] = " "
        for i, s in enumerate(strList):
            for j in range(len(s)):
                arr[i, j] = s[j]

        return arr

    def _convertFortranStringArrayToList(self, fortArray):
        """Undoes the _createFotranStringArray"""
        strList = []
        for ii in range(len(fortArray)):
            if fortArray[ii].dtype == numpy.dtype("|S1"):
                strList.append(b"".join(fortArray[ii]).decode().strip())
            elif fortArray[ii].dtype == numpy.dtype("<U1"):
                strList.append("".join(fortArray[ii]).strip())
            else:
                raise TypeError(f"Unable to convert {fortArray[ii]} of type {fortArray[ii].dtype} to string")

        return strList

    def _readPlot3DSurfFile(self, fileName, convertToTris=True, coordXfer=None):
        """Read a plot3d file and return the points and connectivity in
        an unstructured mesh format"""

        pts = None
        conn = None

        if self.comm.rank == 0:
            f = open(fileName, "r")
            nSurf = numpy.fromfile(f, "int", count=1, sep=" ")[0]
            sizes = numpy.fromfile(f, "int", count=3 * nSurf, sep=" ").reshape((nSurf, 3))
            nPts = 0
            nElem = 0
            for i in range(nSurf):
                curSize = sizes[i, 0] * sizes[i, 1]
                nPts += curSize
                nElem += (sizes[i, 0] - 1) * (sizes[i, 1] - 1)

            # Generate the uncompacted point and connectivity list:
            pts = numpy.zeros((nPts, 3), dtype=self.dtype)
            conn = numpy.zeros((nElem, 4), dtype=int)

            nodeCount = 0
            elemCount = 0
            for iSurf in range(nSurf):
                curSize = sizes[iSurf, 0] * sizes[iSurf, 1]
                for idim in range(3):
                    pts[nodeCount : nodeCount + curSize, idim] = numpy.fromfile(f, "float", curSize, sep=" ")
                # Add in the connectivity.
                iSize = sizes[iSurf, 0]
                for j in range(sizes[iSurf, 1] - 1):
                    for i in range(sizes[iSurf, 0] - 1):
                        conn[elemCount, 0] = nodeCount + j * iSize + i
                        conn[elemCount, 1] = nodeCount + j * iSize + i + 1
                        conn[elemCount, 2] = nodeCount + (j + 1) * iSize + i + 1
                        conn[elemCount, 3] = nodeCount + (j + 1) * iSize + i
                        elemCount += 1
                nodeCount += curSize

            # Next reduce to eliminate the duplicates and form a
            # FE-like structure. This isn't strictly necessary but
            # cheap anyway
            uniquePts, link, nUnique = self.adflow.utils.pointreduce(pts.T, 1e-12)
            pts = uniquePts[:, 0:nUnique].T

            # coordinate transformation for the coordinates if we have any.
            # this is where the user can rotate, translate, or generally manipulate
            # the coordinates coming from the plot3d file before they are used
            if coordXfer is not None:
                # this callback function has the same call signature with the
                # method described in the DVGeometry addPointSet method:
                # https://github.com/mdolab/pygeo/blob/main/pygeo/parameterization/DVGeo.py
                pts = coordXfer(pts, mode="fwd", applyDisplacement=True)

            # Update the conn
            newConn = numpy.zeros_like(conn)
            for i in range(elemCount):
                for j in range(4):
                    newConn[i, j] = link[int(conn[i, j])]
            conn = newConn

            # Now convert the connectivity to triangles if requested
            # since we actually will be doing the integration on a
            # triangles
            if convertToTris:
                newConn = numpy.zeros((len(conn) * 2, 3))
                for i in range(len(conn)):
                    newConn[2 * i, :] = [conn[i, 0], conn[i, 1], conn[i, 2]]
                    newConn[2 * i + 1, :] = [conn[i, 0], conn[i, 2], conn[i, 3]]

                conn = newConn

        # Now bcast to everyone
        pts = self.comm.bcast(pts)
        conn = self.comm.bcast(conn)

        return pts, conn


class adflowFlowCase(object):
    """
    Class containing the data that ADflow requires to be saved to an
    aeroProblem to permit the analysis of multiple flow cases
    """

    def __init__(self):
        self.stateInfo = None
        self.adjoints = {}
        self.adjointRHS = {}
        self.coords = None
        self.callCounter = -1
        self.disp = None
        self.oldWinf = None
        self.ank_cfl = None


class adflowUserFunc(object):
    """Class containing the user-supplied function information"""

    def __init__(self, funcName, functions, callBack):
        self.funcName = funcName.lower()
        self.functions = functions
        self.callBack = callBack
        self.funcs = None

    def evalFunctions(self, funcs):
        # Cache the input funcs for reference
        self.funcs = funcs
        self.callBack(funcs)
        # Make sure the funcName was actually added:
        if self.funcName not in funcs:
            raise Error(
                "The func '%s' (must be lower-case) was " "not supplied from user-supplied function." % self.funcName
            )

    def evalFunctionsSens(self):
        # We need to get the derivative of 'self.funcName' as a
        # function of the 'self.functions'
        refFuncs = copy.deepcopy(self.funcs)
        deriv = OrderedDict()

        for f in self.functions:
            funcs = refFuncs.copy()
            funcs[f] += 1e-40j
            self.callBack(funcs)
            deriv[f] = numpy.imag(funcs[self.funcName]) / 1e-40

        return deriv
