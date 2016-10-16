#!/usr/bin/python
from __future__ import print_function
from __future__ import division
from six import iteritems
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
import numpy
from mpi4py import MPI
from petsc4py import PETSc
from baseclasses import AeroSolver, AeroProblem
from . import MExt
from pprint import pprint as pp

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyADFLOW Error: '
        i = 17
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
        print(msg)
        Exception.__init__(self)

class ADFLOWWarning(object):
    """
    Format a warning message
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pyADFLOW Warning: '
        i = 19
        for word in message.split():
            if len(word) + i + 1 > 78: # Finish line and start new one
                msg += ' '*(78-i)+'|\n| ' + word + ' '
                i = 1 + len(word)+1
            else:
                msg += word + ' '
                i += len(word)+1
        msg += ' '*(78-i) + '|\n' + '+'+'-'*78+'+'+'\n'
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
    options : dictionary
        The list of options to use with ADflow. This keyword arguement
        is NOT OPTIONAL. It must always be provided. It must contain, at least
        the 'gridFile' entry for the filename of the grid to load
    debug : bool
        Set this flag to true when debugging with a symbolic
        debugger. The MExt module deletes the copied .so file when not
        required which causes issues debugging.
    dtype : str
        String type for float: 'd' or 'D'. Not needed to be uset by user.
        """
    def __init__(self, comm=None, options=None, debug=False, dtype='d'):

        # Load the compiled module using MExt, allowing multiple
        # imports
        try:
            self.adflow
        except:
            curDir = os.path.dirname(os.path.realpath(__file__))
            self.adflow = MExt.MExt('libadflow', [curDir], debug=debug)._module

        # Information for base class:
        name = 'ADFLOW'
        category = 'Three Dimensional CFD'
        informs = {}

        # Load all the option/objective/DV information:
        defOpts = self._getDefOptions()
        self.optionMap, self.moduleMap = self._getOptionMap()
        self.ignoreOptions, self.deprecatedOptions, self.specialOptions = \
                           self._getSpecialOptionLists()
        self.imOptions = self._getImmutableOptions()

        self.possibleAeroDVs, self.basicCostFunctions = (
            self._getObjectivesAndDVs())

        # Now add the group for each of the "basic" cost functions:
        self.adflowCostFunctions = {}
        for key in self.basicCostFunctions:
            self.adflowCostFunctions[key] = [None, key]

        # This is the real solver so dtype is 'd'
        self.dtype = dtype

        # Next set the MPI Communicators and associated info
        if comm is None:
            comm = MPI.COMM_WORLD

        self.comm = comm
        self.adflow.communication.adflow_comm_world = self.comm.py2f()
        self.adflow.communication.adflow_comm_self = MPI.COMM_SELF.py2f()
        self.adflow.communication.sendrequests = numpy.zeros(self.comm.size)
        self.adflow.communication.recvrequests = numpy.zeros(self.comm.size)
        self.myid = self.adflow.communication.myid = self.comm.rank
        self.adflow.communication.nproc = self.comm.size

        # Initialize the inherited aerosolver.
        if options is None:
            raise Error("The 'options' keyword argument must be passed "
                        "adflow. The options dictionary must contain (at least) "
                        "the gridFile entry for the grid.")

        # Set all internal adflow default options before we set anything from python
        self.adflow.inputparamroutines.setdefaultvalues()

        AeroSolver.__init__(self, name, category, defOpts, informs,
                            options=options)

        # Update turbresscale depending on the turbulence model specified
        self._updateTurbResScale()

        # Initialize petec in case the user has not already
        self.adflow.adjointutils.initializepetsc()

        # Set the stand-alone adflow flag to false...this changes how
        # terminate calls are handled.
        self.adflow.iteration.standalonemode = False

        # Set the frompython flag to true... this also changes how
        # terminate calls are handled
        self.adflow.killsignals.frompython = True

        # Dictionary of design varibales and their index
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
        # Matrix Setup Flag
        self.adjointSetup = False

        # Write the intro message
        self.adflow.utils.writeintromessage()

        # Remind the user of all the adflow options:
        self.printCurrentOptions()

        # Do the remainder of the operations that would have been done
        # had we read in a param file
        self.adflow.iteration.deforming_grid = True

        # In order to properly initialize we need to have mach number
        # and a few other things set. Just create a dummy aeroproblem,
        # use it, and then it will be deleted.

        dummyAP = AeroProblem(name='dummy', mach=0.5, altitude=10000.0,
                              areaRef=1.0, chordRef=1.0, alpha=0.0, degreePol=0,
                              coefPol=[0, 0], degreeFourier=1,
                              omegaFourier=6.28, sinCoefFourier=[0, 0],
                              cosCoefFourier=[0, 0])

        self._setAeroProblemData(dummyAP, firstCall=True)
        # Now set it back to None so the user is none the wiser
        self.curAP = None

        # Finally complete loading
        self.adflow.inputparamroutines.dummyreadparamfile()
        if self.getOption('partitionOnly'):
            self.adflow.partitioning.partitionandreadgrid(True)
            return

        self.adflow.partitioning.partitionandreadgrid(False)
        self.adflow.preprocessingapi.preprocessing()

        # Do explict hole cutting.
        ncells = self.adflow.adjointvars.ncellslocal[0]
        ntime  = self.adflow.inputtimespectral.ntimeintervalsspectral
        n  = ncells*ntime

        # Call the user supplied callback if necessary

        cutCallBack = self.getOption('cutCallBack')
        flag = numpy.zeros(n)
        if cutCallBack is not None:
            xCen = self.adflow.utils.getcellcenters(1, n).T
            cutCallBack(xCen, flag)
        self.adflow.preprocessingapi.preprocessingoverset(flag)

        self.adflow.tecplotio.initializeliftdistributiondata()
        self.adflow.initializeflow.updatebcdataalllevels()
        self.adflow.initializeflow.initflow()

        famList = self._processFortranStringArray(
            self.adflow.surfacefamilies.famnames)

        # Add the initial families that already exist in the CGNS
        # file.
        for i in range(len(famList)):
            self.families[famList[i]] = [i+1]

        # Add the special "all surfaces" family.
        self.allFamilies = 'allSurfaces'
        self.addFamilyGroup(self.allFamilies, famList)

        # We also need to know about which surfaces are walls. Pull
        # out that information from fortran and set the special
        # familyGroup. 
        wallList, nWalls = self.adflow.surfaceutils.getwalllist(len(famList))
        wallList = wallList[0:nWalls]
        wallListFam = []
        for i in range(len(wallList)):
            wallListFam.append(famList[wallList[i]-1])
        self.allWallsGroup = 'allWalls'
        self.addFamilyGroup(self.allWallsGroup, wallListFam)

        # Set the design families if given, otherwise default to all
        # walls
        self.designFamilyGroup = self.getOption('designSurfaceFamily')
        if self.designFamilyGroup is None:
            self.designFamilyGroup = self.allWallsGroup

        # Set the mesh families if given, otherwise default to all
        # walls
        self.meshFamilyGroup = self.getOption('meshSurfaceFamily')
        if self.meshFamilyGroup is None:
            self.meshFamilyGroup = self.allWallsGroup

        self.coords0 = self.getSurfaceCoordinates(self.allFamilies)


    def getSolverMeshIndices(self):
        '''
        Get the list of indices to pass to the mesh object for the
        volume mesh mapping
        '''
        ndof_1_instance = self.adflow.adjointvars.nnodeslocal[0]*3
        meshInd = self.adflow.warping.getcgnsmeshindices(ndof_1_instance)

        return meshInd

    def setDisplacements(self, aeroProblem, dispFile):
        """
        This function allows the user to perform aerodyanmic
        analysis/optimization while using a fixed set of displacements
        computed from a previous structural analysis. Essentially this
        allows the jig shape to designed, but performing anlysis on
        the flying shape. Note that the fixed set of displacements do
        not affect the sensitivities.
        Parameters
        ----------
        aeroProblem : aeroProblem class
           The AP object that the displacements should be applied to.
        dispFile : str
           The file contaning the displacments. This file should have
           been obtained from TACS
        """
        self.setAeroProblem(aeroProblem)

        # All processors read the file
        f = open(dispFile,'r')
        n = int(f.readline())
        X = []
        D = []
        for i in range(n):
            aux = f.readline().split()
            X.append([float(aux[0]), float(aux[1]), float(aux[2])])
            D.append([float(aux[3]), float(aux[4]), float(aux[5])])
        f.close()

        localX = self.getSurfaceCoordinates()
        localD = numpy.zeros_like(localX)
        # Now we need to search each localX in X to find the corresponding D
        try:
            from scipy.spatial import KDTree
        except:
            raise Error('scip.spatial must be available to use setDisplacements')
        tree = KDTree(numpy.array(X))
        d, index = tree.query(localX)
        for j in range(len(localX)):
            localD[j] = D[index[j]]

        # Fianlly we have the local displacements for this processor we save
        self.curAP.adflowData.disp = localD

    def addLiftDistribution(self, nSegments, direction,
                            groupName=None):
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

        groupTag = '%s: '% groupName
        famList = self._getFamilyList(groupName)
        direction=direction.lower()
        if direction not in ['x','y','z']:
            if direction not in ['x','y','z']:
                raise Error("direction must be one of 'x', 'y', or 'z'")

        if direction == 'x':
            dirVec = [1.0, 0.0, 0.0]
            dirInd = 1
        elif direction == 'y':
            dirVec = [0.0, 1.0, 0.0]
            dirInd = 2
        else:
            dirVec = [0.0, 0.0, 1.0]
            dirInd = 3

        distName = 'LiftDist_%2.2d %s: %s normal'% (
            self.nLiftDist + 1, groupTag, direction)
        self.nLiftDist += 1

        self.adflow.tecplotio.addliftdistribution(
            nSegments, dirVec, dirInd, distName, famList)

    def addSlices(self, direction, positions, sliceType='relative',
                  groupName=None):
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
            Relative slices are 'sliced' at the beginning and then parametricly
            move as the geometry deforms. As a result, the slice through the
            geometry may not remain planar. An abolute slice is re-sliced for
            every out put so is always exactly planar and always at the initial
            position the user indicated.
        groupName : str
             The family to use for the slices. Default is None corresponding to all
             wall groups.
            """

        # Create the zipper mesh if not done so
        self._createZipperMesh()

        # Determine the families we want to use
        if groupName is None:
            groupName = self.allWallsGroup
        groupTag = '%s: '% groupName

        famList = self._getFamilyList(groupName)
        direction = direction.lower()
        if direction not in ['x', 'y', 'z']:
            raise Error("'direction' must be one of 'x', 'y', or 'z'")

        sliceType = sliceType.lower()
        if sliceType not in ['relative', 'absolute']:
            raise Error("'sliceType' must be 'relative' or 'absolute'.")

        positions = numpy.atleast_1d(positions)
        N = len(positions)
        tmp = numpy.zeros((N, 3),self.dtype)
        if direction == 'x':
            tmp[:, 0] = positions
            dirVec = [1.0, 0.0, 0.0]
        elif direction == 'y':
            tmp[:, 1] = positions
            dirVec = [0.0, 1.0, 0.0]
        elif direction == 'z':
            tmp[:, 2] = positions
            dirVec = [0.0, 0.0, 1.0]

        for i in xrange(len(positions)):
            # It is important to ensure each slice get a unique
            # name...so we will number sequentially from pythhon
            j = self.nSlice + i + 1
            if sliceType == 'relative':
                sliceName = 'Slice_%4.4d %s Para Init %s=%7.3f'% (
                    j, groupTag, direction, positions[i])
                self.adflow.tecplotio.addparaslice(sliceName, tmp[i], dirVec, famList)
            else:
                sliceName = 'Slice_%4.4d %s Absolute %s=%7.3f'% (
                    j, groupTag, direction, positions[i])
                self.adflow.tecplotio.addabsslice(sliceName, tmp[i], dirVec, famList)

        self.nSlice += N

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
        the addFunction() routine. See that routine for more documentation. """

        if names is None:
            names = [None]*len(funcNames)

        if len(funcNames) != len(groupNames) or len(funcNames) != len(names):
            raise Error("funcNames, groupNames, and names all have to be "
                        "lists of the same length")

        newFuncNames = []
        for i in range(len(funcNames)):
            funcName = funcNames[i]
            groupName = groupNames[i]
            name = names[i]

            # First make sure the supplied function is already known to adflow
            if funcName.lower() not in self.basicCostFunctions:
                raise Error('Supplied function name is not known to ADflow')

            # Check if the goupName has been added to the mesh.
            if groupName is not None:
                if groupName not in self.families:
                    raise Error("'%s' is not a family in the CGNS file or has not been added"
                                " as a combination of families"%(groupName))

            if name is None:
                adflowFuncName = '%s_%s'%(funcName, groupName)
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
            a  = numpy.sqrt(self.adflow.flowvarrefstate.gammainf*\
                                self.adflow.flowvarrefstate.pinfdim/ \
                                self.adflow.flowvarrefstate.rhoinfdim)
            V = (self.adflow.inputphysics.machgrid+self.adflow.inputphysics.mach)*a

            p = aeroProblem.phat*V/aeroProblem.spanRef
            q = aeroProblem.qhat*V/aeroProblem.chordRef
            r = aeroProblem.rhat*V/aeroProblem.spanRef
            if aeroProblem.liftIndex == 2:
                rotations = [p,r,q]
            elif aeroProblem.liftIndex == 3:
                rotations = [p,q,r]
            else:
                raise Error("Invalid lift direction. Must be 2 or 3 for "
                            "steady rotations and specifying an aeroProblem")

        else:
            rotations = rotRate

        if cgnsBlocks is None:
            # By Default set all the blocks:
            cgnsBlocks = numpy.arange(1, self.adflow.cgnsgrid.cgnsndom+1)

        self.adflow.updaterotationrate(rotCenter, rotations, cgnsBlocks)
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
        aeroProblem : pyAero_problem class
            The complete description of the problem to solve
        mdCallBack : python function
            Call back for doing unsteady aeroelastic. Currently
            not supported.
        writeSolution : bool
            Flag to override any solution writing parameters. This
            is used in a multidisciplinary enviornment when the outer
            solver can suppress all I/O during intermediate solves.
            """

        # Make sure the user isn't trying to solve a slave
        # aeroproblem. Cannot do that
        if hasattr(aeroProblem, 'isSlave'):
            if aeroProblem.isSlave:
                raise Error('Cannot solve an aeroProblem created as a slave')

        # Get option about adjoint memory
        releaseAdjointMemory = kwargs.pop('relaseAdjointMemory', True)

        # Set the aeroProblem
        self.setAeroProblem(aeroProblem, releaseAdjointMemory)

        # Set filenames for possible forced write during solution
        self._setForcedFileNames()

        # Remind the users of the modified options:
        if self.getOption('printIterations'):
            self.printModifiedOptions()

        # Possibly release adjoint memory if not already done
        # so. Don't remove this. If switching problems, this is done
        # in setAeroProblem, but for when not switching it must be
        # done here regardless. 
        if releaseAdjointMemory:
            self.releaseAdjointMemory()

        # Clear out any saved adjoint RHS since they are now out of
        # data. Also increment the counter for this case.
        self.curAP.adflowData.adjointRHS = {}
        self.curAP.adflowData.callCounter += 1

        # --------------------------------------------------------------
        # Setup interation arrays ---- don't touch this unless you
        # REALLY REALLY know what you're doing!

        # iterTot is non-zero which means we already did a solution,
        # so we can run on the finest level.
        if self.adflow.iteration.itertot != 0:
            self.adflow.inputiteration.mgstartlevel = 1

        self.adflow.iteration.itertot = 0
        desiredSize = self.adflow.inputiteration.nsgstartup + \
                      self.adflow.inputiteration.ncycles
        self.adflow.utils.allocconvarrays(desiredSize)
        # --------------------------------------------------------------

        if self.getOption('equationMode').lower() == 'unsteady':
            self.adflow.utils.alloctimearrays(self.getOption('nTimeStepsFine'))
            self._setUnsteadyFileParameters()

        # Mesh warp may have already failed:
        if not self.adflow.killsignals.fatalfail:

            if (self.getOption('equationMode').lower() == 'steady' or
                self.getOption('equationMode').lower() == 'time spectral'):
                self.updateGeometryInfo()

                # Check to see if the above update routines failed.
                self.adflow.killsignals.fatalfail = \
                    self.comm.allreduce(
                    bool(self.adflow.killsignals.fatalfail), op=MPI.LOR)

        if self.adflow.killsignals.fatalfail:
            print("Fatal failure during mesh warp! Bad mesh is "
                  "written in output directory as failed_mesh.cgns")
            fileName = os.path.join(self.getOption('outputDirectory'),
                                    'failed_mesh.cgns')
            self.writeMeshFile(fileName)
            self.curAP.fatalFail = True
            self.curAP.solveFailed = True
            return

        # We now now the mesh warping was ok so reset the flags:
        self.adflow.killsignals.routinefailed =  False
        self.adflow.killsignals.fatalfail = False

        t1 = time.time()
        
        # Solve the equations in appropriate mode
        mode = self.getOption('equationMode').lower()
        if mode in ['steady', 'time spectral']:
            # In steady mode, the wall temperature has to be set after the solver
            # switches to current AP, otherwise the data will be overwritten.
            wallTemperature = kwargs.pop('wallTemperature', None)
            groupName       = kwargs.pop('groupName', None)
            if wallTemperature is not None:
                self.setWallTemperature(wallTemperature, groupName=groupName)
            
            self.adflow.solvers.solver()
        elif mode == 'unsteady':
            # Initialization specific for unsteady simulation
            self.adflow.solvers.solverunsteadyinit()

            # If ADflow is to be coupled with other solvers, we are done
            if self.getOption('coupledsolution'):
                return
            # Otherwise, the problem is solved within this function
            else:
                # Get the callback functions
                surfaceMeshCallback = kwargs.pop('surfaceMeshCallback', None)
                volumeMeshCallback  = kwargs.pop('volumeMeshCallback', None)
                
                # Call solver wrapper for unsteady simulations
                self._solverunsteady(surfaceMeshCallback = surfaceMeshCallback,
                                     volumeMeshCallback  = volumeMeshCallback)
        else:
            raise Error("Mode {0} not supported".format(mode))

        # Save the states into the aeroProblem
        self.curAP.adflowData.stateInfo = self._getInfo()

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

        if self.getOption('printTiming') and self.comm.rank == 0:
            print('Solution Time: %10.3f sec'% solTime)

        # Post-Processing
        # --------------------------------------------------------------
        # Solution is written if writeSolution argument true or does not exist
        if kwargs.pop('writeSolution', True):
            self.writeSolution()

        if self.getOption('TSStability'):
            self.computeStabilityParameters()

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
        nTimeStepsFine = self.getOption('ntimestepsfine')
        for tdx in xrange(1, nTimeStepsFine+1):
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

    def advanceTimeStepCounter(self):
        """
        Advance one unit of timestep and physical time.
        """
        self.adflow.monitor.timestepunsteady = self.adflow.monitor.timestepunsteady + 1
        self.adflow.monitor.timeunsteady     = self.adflow.monitor.timeunsteady     + \
            self.adflow.inputunsteady.deltat

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
        if self.getOption('printTiming') and self.comm.rank == 0:
            print('Timestep solution time: {0:10.3f} sec'.format(t2-t1))

    def evalFunctions(self, aeroProblem, funcs, evalFuncs=None, 
                      ignoreMissing=False):
        """
        Evaluate the desired functions given in iterable object,
        'evalFuncs' and add them to the dictionary 'funcs'. The keys
        in the funcs dictioary will be have an _<ap.name> appended to
        them.

        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to to get the solution for

        funcs : dict
            Dictionary into which the functions are saved.

        evalFuncs : iterable object containing strings
          If not None, use these functions to evaluate.

        ignoreMissing : bool
            Flag to supress checking for a valid function. Please use
            this option with caution.

        Examples
        --------
        >>> funcs = {}
        >>> CFDsolver(ap)
        >>> CFDsolver.evalFunctions(ap1, funcs, ['cl', 'cd'])
        >>> funcs
        >>> # Result will look like (if aeroProblem, ap1, has name of 'c1'):
        >>> # {'cl_c1':0.501, 'cd_c1':0.02750}
        """

        # Set the AP
        self.setAeroProblem(aeroProblem)
        if evalFuncs is None:
            evalFuncs = self.curAP.evalFuncs

        # We need to determine how many different groupings we have,
        # since we will have to call getSolution for each *unique*
        # function grouping. We can also do the error checking
        groupMap = {}
        for f in evalFuncs:
            if f.lower() in self.adflowCostFunctions:
                group = self.adflowCostFunctions[f][0]
                basicFunc = self.adflowCostFunctions[f][1]
                if group not in groupMap:
                    groupMap[group] = [[basicFunc, f]]
                else:
                    groupMap[group].append([basicFunc, f])
            else:
                if not ignoreMissing:
                    raise Error('Supplied function %s is not known to ADflow.'%f)

        # Now we loop over the unique groups calling the required
        # getSolution, there may be just one. No need for error
        # checking here.
        for group in groupMap:

            # Call the lower level getSolution() function
            res = self.getSolution(groupName=group)

            # g contains the "basic function" (index 0) and the actual
            # function name (index 1)
            for g in groupMap[group]:
                key = self.curAP.name + '_%s'% g[1]
                self.curAP.funcNames[g[1]] = key
                funcs[key] = res[g[0]]

    def evalFunctionsSens(self, aeroProblem, funcsSens, evalFuncs=None):
        """
        Evaluate the sensitivity of the desired functions given in
        iterable object,'evalFuncs' and add them to the dictionary
        'funcSens'. The keys in the funcs dictioary will be have an
        _<ap.name> appended to them.

        Parameters
        ----------
        funcSens : dict
            Dictionary into which the function derivatives are saved.

        evalFuncs : iterable object containing strings
            The functions the user wants the derivatives of

        evalFuncs : iterable object containing strings
            The additaion functions the user wants returned that are
            not already defined in the aeroProblem

        Examples
        --------
        >>> funcSens = {}
        >>> CFDsolver.evalFunctionsSens(ap1, funcSens, ['cl', 'cd'])
        """

        # This is the one and only gateway to the getting derivatives
        # out of adflow. If you want a derivative, it should come from
        # here.

        self.setAeroProblem(aeroProblem)

        if evalFuncs is None:
            evalFuncs = self.curAP.evalFuncs
        else:
            evalFuncs = list(evalFuncs)

        # Do the functions one at a time:
        for f in evalFuncs:
            if f.lower() not in self.adflowCostFunctions:
                raise Error('Supplied %s function is not known to ADflow.'%f)

            if self.comm.rank == 0:
                print('Solving adjoint: %s'%f)

            key = self.curAP.name + '_%s'% f
            ptSetName = self.curAP.ptSetName

            # Set dict structure for this derivative
            funcsSens[key] = {}

            # Solve adjoint equation
            self.solveAdjoint(aeroProblem, f)

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
            funcsBar = {f.lower():1.0}

            # Compute everything and update into the dictionary
            funcsSens[key].update(self.computeJacobianVectorProductBwd(
                resBar=psi, funcsBar=funcsBar, xDvDeriv=True))

    def solveCL(self, aeroProblem, CLStar, alpha0=None,
                delta=0.5, tol=1e-3, autoReset=True, CLalphaGuess=None,
                maxIter = 20, nReset=25):

        """This is a simple secant method search for solving for a
        fixed CL. This really should only be used to determine the
        starting alpha for a lift constraint in an optimization.

        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to solve
        CLStar : float
            The desired target CL
        alpha0 : angle (deg)
            Initial guess for secant seach (deg). If None, use the
            value in the aeroProblem
        delta : angle (deg)
            Initial step direction for secant search
        tol : float
            Desired tolerance for CL
        autoReset : bool
            Flag to reset flow between solves. The Euler NK method has
            issues when only the alpha is changed (Martois effect we
            think). This will reset the flow after each solve which
            solves this problem. Not necessary (or desired) when using
            the RK solver.
        CLalphaGuess : float or None
            The user can provide an estimate for the lift curve slope
            in order to accelerate convergence. If the user supply a
            value to this option, it will not use the delta value anymore
            to select the angle of attack of the second run. The value
            should be in 1/deg.

        Returns
        -------
        None, but the correct alpha is stored in the aeroProblem

        """
        self.setAeroProblem(aeroProblem)
        if alpha0 is not None:
            aeroProblem.alpha = alpha0

        # Set the startign value
        anm2 = aeroProblem.alpha

        # Print alpha
        if self.comm.rank == 0:
            print ('Current alpha is: ', aeroProblem.alpha)

        self.__call__(aeroProblem, writeSolution=False)
        self.curAP.adflowData.callCounter -= 1
        sol = self.getSolution()
        fnm2 = sol['cl'] - CLStar

        if CLalphaGuess is None:
            # Use the delta option to define the next Aoa
            if fnm2 < 0:
                anm1 = alpha0 + abs(delta)
            else:
                anm1 = alpha0 - abs(delta)
        else:
            # Use CLalphaGuess option to define the next Aoa
            anm1 = alpha0 - fnm2/CLalphaGuess

        # Check for convergence if the starting point is already ok.
        if abs(fnm2) < tol:
            return

        # Overwrite the RK options momentarily.
        # We need to do this to avoid using NK right at the beggining
        # of the new AoA iteration. When we change the AoA, we only change
        # the residuals at the boundaries, and the Newton solver will not
        # allow these residuals to change the rest of the solution.
        # A few RK iterations allow the total residual to "go uphill",
        # so that we can converge to a new solution.
        minIterSave = self.getOption('nRKReset')
        rkresetSave = self.getOption('rkreset')
        self.setOption('nRKReset', nReset)
        self.setOption('rkreset', True)

        # Secant method iterations
        for iIter in range(maxIter):
            # We may need to reset the flow since changing strictly
            # alpha leads to problems with the NK solver
            if autoReset:
                self.resetFlow(aeroProblem)

            # Set current alpha
            aeroProblem.alpha = anm1

            # Print alpha
            if self.comm.rank == 0:
                print ('Current alpha is: ', aeroProblem.alpha)

            # Solve for n-1 value (anm1)
            self.__call__(aeroProblem, writeSolution=False)
            self.curAP.adflowData.callCounter -= 1
            sol = self.getSolution()
            fnm1 = sol['cl'] - CLStar

            # Check for convergence
            if abs(fnm1) < tol:
                break

            # Secant Update
            anew = anm1 - fnm1*(anm1-anm2)/(fnm1-fnm2)

            # Shift n-1 values to n-2 values
            fnm2 = fnm1
            anm2 = anm1

            # Se the n-1 alpha value from update
            anm1 = anew

        # Restore the min iter option given initially by user
        self.setOption('nRKReset', minIterSave)
        self.setOption('rkreset', rkresetSave)

    def solveTrimCL(self, aeroProblem, trimFunc, trimDV, dvIndex,
                    CLStar, trimStar=0.0, alpha0=None, trim0=None, da=1e-3,
                    deta=1e-2, tol=1e-4, nIter=10, Jac0=None, liftFunc='cl'):
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
            Initial alpha step for jacobian
        deta : float
            Initial  stet in the 'eta' or trim dv function
        tol : float
            Tolerance for trimCL solve solution
        nIter : int
            Maximum number of iterations.
        Jac0 : 2x2 numpy array
            Initial guess for the trim-cl jacobian. Usually obtained
            from a previous analysis and saves two function
            evaluations to produce the intial jacobian.
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
            F = numpy.array([funcs['%s_%s'%(self.curAP.name, liftFunc)]-CLstar,
                             funcs['%s_%s'%(self.curAP.name, trimFunc)]])
            return F

        # Generate initial point
        Xn = numpy.array([alpha0, trim0])
        Fn, sol = Func(Xn, CLStar)

        # Next we generate the initial jacobian if we haven't been
        # provided with one.
        if Jac0 is not None:
            J = Jac0.copy()
        else:
            J = numpy.zeros((2,2))

            # Perturb alpha
            Xn[0] += da
            Fpda = Func(Xn, CLStar)
            J[:, 0] = (Fpda - Fn)/da
            Xn[0] -= da

            # Perturb eta:
            Xn[1] += deta
            Fpdeta = Func(Xn, CLStar)
            J[:, 1] = (Fpdeta - Fn)/deta
            Xn[1] -= deta

        # Main iteration loop
        for jj in range(nIter):
            if self.comm.rank == 0:
                print ('Fn:', Fn)

            # We now have a point and a jacobian...Newton's Method!
            Xnp1 = Xn - numpy.linalg.solve(J, Fn)
            if self.comm.rank == 0:
                print ("Xnp1:", Xnp1)

            # Solve the new Xnp1
            Fnp1 = Func(Xnp1, CLStar)
            if self.comm.rank == 0:
                print ("Fnp1:", Fnp1)

            # Update the jacobian using Broyden's method
            dx = Xnp1 - Xn
            dF = Fnp1 - Fn
            J = J + numpy.outer((dF - numpy.dot(J, dx))/(numpy.linalg.norm(dx)**2), dx)

            if self.comm.rank == 0:
                print ("New J:", J)

            # Shuffle the Fn and Xn backwards
            Fn = Fnp1.copy()
            Xn = Xnp1.copy()

            # Check for convergence
            if numpy.linalg.norm(Fn) < tol:
                if self.comm.rank == 0:
                    print ('Converged!', jj, Fn, Xn)
                break
        return Xn

    def solveSep(self, aeroProblem, sepStar, nIter=10, alpha0=None,
                 delta=0.1, tol=1e-3, expansionRatio=1.2, sepName=None):
        """This is a safe-guarded secant search method to determine
        the alpha that yields a specified value of the sep
        sensor. Since this function is highly nonlinear we use a
        linear search to get the bounding range first.

        Parameters
        ----------
        aeroProblem : pyAero_problem class
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
            sepName = 'sepsensor'

        if alpha0 is None:
            alpha0 = ap.alpha

        # Name of function to use
        funcName = '%s_%s'%(ap.name, sepName)

        if not self.getOption('rkreset') and self.getOption('usenksolve'):
            ADFLOWWarning("nRKReset option is not set. It is usually necessary "
                        "for solveSep() when NK solver is used.")

        # Solve first problem
        ap.alpha = alpha0
        self.__call__(ap, writeSolution=False)
        funcs = {}
        self.evalFunctions(ap, funcs, [sepName])
        f = funcs[funcName] - sepStar
        da = delta

        if self.comm.rank == 0:
            print('+----------------------------------------------+')
            print('| Presolve ')
            print('| Alpha      = %s'%ap.alpha)
            print('| Sep value  = %s'%funcs[funcName])
            print('| F          = %s'%f)
            print('| Sep*       = %s'%sepStar)
            print('+----------------------------------------------+')

        if self.comm.rank == 0:
            print ('Searching for Correct Interval')

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
                print('+----------------------------------------------+')
                print('| Finished It= %d'%i)
                print('| Alpha      = %s'%ap.alpha)
                print('| Sep value  = %s'%funcs[funcName])
                print('| F          = %s'%fnew)
                print('| Sep*       = %s'%sepStar)
                print('+----------------------------------------------+')

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
            print ('Switching to Secant Search')

        for iIter in range(nIter):
            if iIter != 0:
                ap.alpha = anm1
                self.__call__(ap, writeSolution=False)
                self.evalFunctions(ap, funcs, [sepName])
                fnm1 = funcs[funcName] - sepStar

            if self.comm.rank == 0:
                print('+----------------------------------------------+')
                print('| Finished It= %d'%i)
                print('| Alpha      = %s'%ap.alpha)
                print('| Sep value  = %s'%funcs[funcName])
                print('| F          = %s'%fnm1)
                print('| Sep*       = %s'%sepStar)
                print('+----------------------------------------------+')

            # Secant update
            anew = anm1 - fnm1*(anm1 - anm2)/(fnm1-fnm2)

            # Make sure the anew update doesn't go outside (lowAlpha,
            # highAlpha). Just use bisection in this case.
            if anew < lowAlpha:
                anew = 0.5*(anm1 + lowAlpha)
            if anew > highAlpha:
                anew = 0.5*(anm1 + highAlpha)

            # Shift the n-1 values to n-2
            fnm2 = fnm1
            anm2 = anm1

            # Set the new n-1 alpha value from update
            anm1 = anew

            # Finally, convergence check
            if abs(fnm1) < tol:
                break

    def writeSolution(self, outputDir=None, baseName=None, number=None):
        """This is a generic shell function that potentially writes
        the various output files. The intent is that the user or
        calling program can call this file and ADflow write all the
        files that the user has defined. It is recommneded that this
        function is used along with the associated logical flags in
        the options to determine the desired writing procedure

        Optional arguments

        outputDir: Use the supplied output directory

        baseName: Use this supplied string for the base filename. Typically
                  only used from an external solver.
        number: Use the user supplied number to index solution. Again, only
                typically used from an external solver.
                """
        if outputDir is None:
            outputDir = self.getOption('outputDirectory')

        if baseName is None:
            baseName = self.curAP.name

        # If we are numbering solution, it saving the sequence of
        # calls, add the call number
        if number is not None:
            # We need number based on the provided number:
            baseName = baseName + '_%3.3d'% number
        else:
            # if number is none, i.e. standalone, but we need to
            # number solutions, use internal counter
            if self.getOption('numberSolutions'):
                baseName = baseName + '_%3.3d'% self.curAP.adflowData.callCounter

        # Join to get the actual filename root
        base = os.path.join(outputDir, baseName)

        # Get the timestep if this is time dependent solution, else we default
        # to the same values as the nsave values (which is 1) for steady
        # and time-spectral solution to be written,
        if self.getOption('equationMode').lower() == 'unsteady':
            ts = self.adflow.monitor.timestepunsteady
        else:
            ts = 1

        # Now call each of the 4 routines with the appropriate file name:
        if self.getOption('writevolumesolution') and \
          numpy.mod(ts, self.getOption('nsavevolume')) == 0 :
            self.writeVolumeSolutionFile(base + '_vol.cgns')

        if self.getOption('writesurfacesolution') and \
          numpy.mod(ts, self.getOption('nsavesurface')) == 0 :
            self.writeSurfaceSolutionFile(base + '_surf.cgns')

        # ADD NSAVE TEST AROUND THESE AS WELL BUT
        # THIS IS SMALL COMPARED TO OTHER. REFACTOR
        if self.getOption('equationMode').lower() == 'unsteady':
            liftName = base + '_lift_Timestep%4.4d.dat'% ts
            sliceName = base + '_slices_Timestep%4.4d.dat'% ts
            surfName = base + '_surf_Timestep%4.4d.dat' %ts
        else:
            liftName = base + '_lift.dat'
            sliceName = base + '_slices.dat'
            surfName = base + '_surf.dat'

        
        # Get the family list to write for the surface.
        famList = self._getFamilyList(self.getOption('outputSurfaceFamily'))

        # Flag to write the tecplot surface solution or not
        writeSurf = self.getOption('writeTecplotSurfaceSolution')

        # # Call fully compbined fortran routine. 
        self.adflow.tecplotio.writetecplot(sliceName, True, liftName, True,
                                           surfName, writeSurf, famList)
        
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
        fileName += '.cgns'

        # Set Flags for writing
        self.adflow.monitor.writegrid = True
        self.adflow.monitor.writevolume = False
        self.adflow.monitor.writesurface = False

        # Set filename in adflow
        self.adflow.inputio.solfile[:] = ''
        self.adflow.inputio.solfile[0:len(fileName)] = fileName

        self.adflow.inputio.newgridfile[:] = ''
        self.adflow.inputio.newgridfile[0:len(fileName)] = fileName

        # Actual fortran write call. Family list doesn't matter.
        famList = self._getFamilyList(self.allFamilies)
        self.adflow.writesol(famList)

    def writeVolumeSolutionFile(self, fileName, writeGrid=True):
        """Write the current state of the volume flow solution to a CGNS
        file. This is a lower level routine; Normally one should call
        writeSolution().

        Parameters
        ----------
        fileName : str
            Name of the file. Should have .cgns extension.
        writeGrid : bool
            Flag specifying whether the grid should be included or if
            links should be used. Always writing the grid is
            recommended even in cases when it is not strictly necessary.
            Note that if writeGrid == False the volume files do not contain
            any grid coordinates rendering the file useless if a separate
            grid file was written out and is linked to it.
        """
        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += '.cgns'

        # Set Flags for writing
        self.adflow.monitor.writegrid = writeGrid
        self.adflow.monitor.writevolume = True
        self.adflow.monitor.writesurface = False

        # Set fileName in adflow
        self.adflow.inputio.solfile[:] = ''
        self.adflow.inputio.solfile[0:len(fileName)] = fileName

        self.adflow.inputio.newgridfile[:] = ''
        self.adflow.inputio.newgridfile[0:len(fileName)] = fileName

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
        fileName += '.cgns'

        # Set Flags for writing
        self.adflow.monitor.writegrid=False
        self.adflow.monitor.writevolume=False
        self.adflow.monitor.writesurface=True

        # Set fileName in adflow
        self.adflow.inputio.surfacesolfile[:] = ''
        self.adflow.inputio.surfacesolfile[0:len(fileName)] = fileName

        # Actual fortran write call. Fam list matters. 
        famList = self._getFamilyList(self.getOption('outputSurfaceFamily'))
        self.adflow.writesol(famList)

    def writeLiftDistributionFile(self, fileName):
        """Evaluate and write the lift distibution to a tecplot file.

        Parameters
        ----------
        fileName : str
            File of lift distribution. Should have .dat extension.
        """
        # Ensure fileName is .dat even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += '.dat'

        # Actual write command
        sliceName = ""
        surfName = ""
        self.adflow.tecplotio.writetecplot(sliceName, False, fileName, True,
                                         surfName, False)

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
        fileName += '.dat'

        # Actual write command
        sliceName = ""
        surfName = ""
        self.adflow.tecplotio.writetecplot(sliceName, False, fileName, True,
                                         surfName, False)

    def writeForceFile(self, fileName, TS=0, groupName=None,
                       cfdForcePts=None):
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
            pts = self.comm.gather(self.getSurfacePoints(groupName, TS), root=0)
        else:
            pts = self.comm.gather(cfdForcePts)

        # Get the forces and connectivity
        forces = self.comm.gather(self.getForces(groupName, TS=TS), root=0)
        conn, faceSize = self.getSurfaceConnectivity(groupName)
        conn = self.comm.gather(conn, root=0)

        # Write out Data only on root proc:
        if self.myid == 0:
            # First sum up the total number of nodes and elements:
            nPt = 0
            nCell = 0
            for iProc in xrange(len(pts)):
                nPt += len(pts[iProc])
                nCell += len(conn[iProc])

            # Open output file
            f = open(fileName, 'w')

            # Write header with number of nodes and number of cells
            f.write("%d %d\n"% (nPt, nCell))

            # Now write out all the Nodes and Forces (or tractions)
            for iProc in xrange(len(pts)):
                for i in xrange(len(pts[iProc])):
                    f.write('%15.8g %15.8g %15.8g '% (
                            numpy.real(pts[iProc][i, 0]),
                            numpy.real(pts[iProc][i, 1]),
                            numpy.real(pts[iProc][i, 2])))
                    f.write('%15.8g %15.8g %15.8g\n'% (
                            numpy.real(forces[iProc][i, 0]),
                            numpy.real(forces[iProc][i, 1]),
                            numpy.real(forces[iProc][i, 2])))

            # Now write out the connectivity information. We have to
            # be a little careful, since the connectivitiy is given
            # locally per proc. As we loop over the data from each
            # proc, we need to increment the connectivity by the
            # number of nodes we've "used up" so far
            nodeOffset = 0
            for iProc in xrange(len(conn)):
                conn[iProc] += 1
                for i in xrange(len(conn[iProc])):
                    f.write('%d %d %d %d\n'%(
                        conn[iProc][i,0]+nodeOffset,
                        conn[iProc][i,1]+nodeOffset,
                        conn[iProc][i,2]+nodeOffset,
                        conn[iProc][i,3]+nodeOffset))

                nodeOffset += len(pts[iProc])
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
        aeroProblem : pyAero_problem object
            The aeroproblem with the flow information we would like
            to reset the flow to.
            """
        self.setAeroProblem(aeroProblem, releaseAdjointMemory)
        self._resetFlow()

    def _resetFlow(self):
        strLvl = self.getOption('MGStartLevel')
        nLevels = self.adflow.inputiteration.nmglevels
        if strLvl < 0 or strLvl > nLevels:
            strLvl = nLevels

        self.adflow.inputiteration.mgstartlevel = strLvl
        self.adflow.iteration.groundlevel = strLvl
        self.adflow.iteration.currentlevel = strLvl
        self.adflow.monitor.nitercur = 0
        self.adflow.iteration.itertot = 0
        self.adflow.initializeflow.setuniformflow()
        self.adflow.killsignals.routinefailed =  False
        self.adflow.killsignals.fatalfail = False
        self.adflow.nksolver.freestreamresset = False

    def getSolution(self, groupName=None):
        """ Retrieve the solution variables from the solver. Note this
        is a collective function and must be called on all processors
        """

        # Extract the familiy list we want to use for evaluation
        famList = self._getFamilyList(groupName)

        # We should return the list of results that is the same as the
        # possibleObjectives list
        self.adflow.surfaceintegrations.getsolution(1, famList)

        funcVals = self.adflow.costfunctions.funcvalues
        ADflowsolution = {
            'lift':funcVals[self.adflow.costfunctions.costfunclift-1],
            'drag':funcVals[self.adflow.costfunctions.costfuncdrag-1],
            'cl'  :funcVals[self.adflow.costfunctions.costfuncliftcoef-1],
            'cd'  :funcVals[self.adflow.costfunctions.costfuncdragcoef-1],
            'fx'  :funcVals[self.adflow.costfunctions.costfuncforcex-1],
            'fy'  :funcVals[self.adflow.costfunctions.costfuncforcey-1],
            'fz'  :funcVals[self.adflow.costfunctions.costfuncforcez-1],
            'cfx' :funcVals[self.adflow.costfunctions.costfuncforcexcoef-1],
            'cfy' :funcVals[self.adflow.costfunctions.costfuncforceycoef-1],
            'cfz' :funcVals[self.adflow.costfunctions.costfuncforcezcoef-1],
            'mx'  :funcVals[self.adflow.costfunctions.costfuncmomx-1],
            'my'  :funcVals[self.adflow.costfunctions.costfuncmomy-1],
            'mz'  :funcVals[self.adflow.costfunctions.costfuncmomz-1],
            'cmx' :funcVals[self.adflow.costfunctions.costfuncmomxcoef-1],
            'cmy' :funcVals[self.adflow.costfunctions.costfuncmomycoef-1],
            'cmz' :funcVals[self.adflow.costfunctions.costfuncmomzcoef-1],
            'cmzalphadot':funcVals[self.adflow.costfunctions.costfunccmzalphadot-1],
            'cmzalpha'   :funcVals[self.adflow.costfunctions.costfunccmzalpha-1],
            'cm0'        :funcVals[self.adflow.costfunctions.costfunccm0-1],
            'clalphadot' :funcVals[self.adflow.costfunctions.costfuncclalphadot-1],
            'clalpha'    :funcVals[self.adflow.costfunctions.costfuncclalpha-1],
            'cl0'        :funcVals[self.adflow.costfunctions.costfunccl0-1],
            'cfyalphadot':funcVals[self.adflow.costfunctions.costfunccfyalphadot-1],
            'cfyalpha'   :funcVals[self.adflow.costfunctions.costfunccfyalpha-1],
            'cfy0'       :funcVals[self.adflow.costfunctions.costfunccfy0-1],
            'cdalphadot' :funcVals[self.adflow.costfunctions.costfunccdalphadot-1],
            'cdalpha'    :funcVals[self.adflow.costfunctions.costfunccdalpha-1],
            'cd0'        :funcVals[self.adflow.costfunctions.costfunccd0-1],
            'cmzqdot'    :funcVals[self.adflow.costfunctions.costfunccmzqdot-1],
            'cmzq'       :funcVals[self.adflow.costfunctions.costfunccmzq-1],
            'clqdot'     :funcVals[self.adflow.costfunctions.costfuncclqdot-1],
            'clq'        :funcVals[self.adflow.costfunctions.costfuncclq-1],
            'cbend'      :funcVals[self.adflow.costfunctions.costfuncbendingcoef-1],
            'sepsensor'  :funcVals[self.adflow.costfunctions.costfuncsepsensor-1],
            'sepsensoravgx'  :funcVals[self.adflow.costfunctions.costfuncsepsensoravgx-1],
            'sepsensoravgy'  :funcVals[self.adflow.costfunctions.costfuncsepsensoravgy-1],
            'sepsensoravgz'  :funcVals[self.adflow.costfunctions.costfuncsepsensoravgz-1],
            'cavitation' :funcVals[self.adflow.costfunctions.costfunccavitation-1],
            'mdot' :funcVals[self.adflow.costfunctions.costfuncmdot-1],
            'mavgptot':funcVals[self.adflow.costfunctions.costfuncmavgptot-1],
            'mavgttot':funcVals[self.adflow.costfunctions.costfuncmavgttot-1],
            'mavgps':funcVals[self.adflow.costfunctions.costfuncmavgps-1]
            }

        return ADflowsolution

    # =========================================================================
    #   The following routines are public functions however, they should
    #   not need to be used by a user using this class directly. They are
    #   primarly used internally OR by a solver that uses this class
    #   i.e. an Aerostructural solver
    # =========================================================================

    def getSurfaceCoordinates(self, groupName=None):
        # This is an alias for getSurfacePoints
        return self.getSurfacePoints(groupName)

    def getPointSetName(self,apName):
        """
        Take the apName and return the mangled point set name.

        """
        return 'adflow_%s_coords'% apName


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
            raise Error("Cannot set new surface coordinate locations without a mesh"
                        "warping object present.")

        # First get the surface coordinates of the meshFamily in case
        # the groupName is a subset, those values will remain unchanged.
        meshSurfCoords = self.getSurfaceCoordinates(self.meshFamilyGroup)
        meshSurfCoords = self.mapVector(coordinates, groupName,
                                        self.meshFamilyGroup, meshSurfCoords)
        self.mesh.setSurfaceCoordinates(meshSurfCoords)

    def setAeroProblem(self, aeroProblem, releaseAdjointMemory=True):
        """Set the supplied aeroProblem to be used in ADflow"""

        ptSetName = 'adflow_%s_coords'% aeroProblem.name

        newAP = False
        # Tell the user if we are switching aeroProblems
        if self.curAP != aeroProblem:
            newAP = True
            if self.comm.rank == 0:
                print('+'+'-'*70+'+')
                print('|  Switching to Aero Problem: %-41s|'% aeroProblem.name)
                print('+'+'-'*70+'+')

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

            for key in self.curAP.savedOptions['adflow']:
                # Saved Val: This is the main option value when the
                # aeroProblem set it's own option
                # setVal: This is the value the aeroProblem set itself
                # curVal: This is the actual value currently set

                savedVal = self.curAP.savedOptions['adflow'][key]
                setVal = self.curAP.solverOptions['adflow'][key]
                curVal = self.getOption(key)

                if curVal == setVal:
                    # Restore the saved value, if it still what the
                    # aeroProblem had set
                    self.setOption(key, savedVal)

        # Now check if we have an DVGeo object to deal with:
        if self.DVGeo is not None:
            # DVGeo appeared and we have not embedded points!
            if not ptSetName in self.DVGeo.points:
                coords0 = self.mapVector(self.coords0, self.allFamilies,
                                         self.designFamilyGroup)
                self.DVGeo.addPointSet(coords0, ptSetName)

            # Check if our point-set is up to date:
            if not self.DVGeo.pointSetUpToDate(ptSetName):
                coords = self.DVGeo.update(ptSetName, config=aeroProblem.name)

                # Potentially add a fixed set of displacements to it.
                if aeroProblem.adflowData.disp is not None:
                    coords += self.curAP.adflowData.disp
                self.setSurfaceCoordinates(coords, self.designFamilyGroup)

        self._setAeroProblemData(aeroProblem)

        # Note that it is safe to call updateGeometryInfo since it
        # only updates the mesh if necessary
        self.updateGeometryInfo()

        # Create the zipper mesh if not done so
        self._createZipperMesh()

        # If the state info none, initialize to the supplied
        # restart file or the free-stream values by calling
        # resetFlow()
        if aeroProblem.adflowData.stateInfo is None:
            if self.getOption('restartFile') is not None:
                self.adflow.inputiteration.mgstartlevel = 1
                self.adflow.initializeflow.initflowrestart()
                # Save this state information
                aeroProblem.adflowData.stateInfo = self._getInfo()
            else:
                self._resetFlow()
     
        # Now set the data from the incomming aeroProblem:
        stateInfo = aeroProblem.adflowData.stateInfo
        if stateInfo is not None and newAP:
            self._setInfo(stateInfo)

        self.adflow.killsignals.routinefailed = False
        self.adflow.killsignals.fatalFail = False

        # We are now ready to associate self.curAP with the supplied AP
        self.curAP = aeroProblem
        self.curAP.adjointRHS = None

        # Destroy the ANK/NK solver and the adjoint memory
        if newAP:
            self.adflow.nksolver.destroynksolver()
            self.adflow.anksolver.destroyanksolver()
            if releaseAdjointMemory:
                self.releaseAdjointMemory()

    def _setAeroProblemData(self, aeroProblem, firstCall=False):
        """
        After an aeroProblem has been associated with self.cuAP, set
        all the updated information in ADflow."""

        # Set any additional adflow options that may be defined in the
        # aeroproblem. While we do it we save the options that we've
        # modified if they are different than the current option.
        AP = aeroProblem
        try:
            AP.savedOptions
        except:
            AP.savedOptions = {'adflow':{}}


        if 'adflow' in AP.solverOptions:
            for key in AP.solverOptions['adflow']:
                curVal = self.getOption(key)
                overwriteVal =  AP.solverOptions['adflow'][key]
                if overwriteVal != curVal:
                    self.setOption(key, overwriteVal)
                    AP.savedOptions['adflow'][key] = curVal

        alpha = AP.alpha
        beta = AP.beta
        beta = AP.beta
        mach = AP.mach
        machRef = AP.machRef
        machGrid = AP.machGrid

        xRef = AP.xRef; yRef = AP.yRef; zRef = AP.zRef
        xRot = AP.xRot; yRot = AP.yRot; zRot = AP.zRot
        areaRef = AP.areaRef
        chordRef = AP.chordRef
        liftIndex = self.getOption('liftIndex')

        if (AP.T is None or AP.P is None or AP.rho is None or
            AP.V is None or AP.mu is None):
            raise Error("Insufficient information is given in the "
                        "aeroProblem to determine physical state. "
                        "See AeroProblem documentation for how to "
                        "specify complete aerodynamic states.")

        if self.dtype == 'd':
            mach = numpy.real(mach)

        if self.dtype == 'd':
            T = numpy.real(AP.T)
            P = numpy.real(AP.P)
            rho = numpy.real(AP.rho)
            V = numpy.real(AP.V)
            mu = numpy.real(AP.mu)

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
            V = AP.V
            mu = AP.mu

            SSuthDim = AP.SSuthDim
            muSuthDim = AP.muSuthDim
            TSuthDim = AP.TSuthDim
            RGasDim = AP.R
            gammaConstant = AP.gamma
            Pr = AP.Pr

        # Do some checking here for things that MUST be specified:
        if AP.mach is None:
            raise Error("'mach' number must be specified in the aeroProblem"
                        " for ADflow.")
        if areaRef is None:
            raise Error("'areaRef' must be specified in aeroProblem"
                        " for ADflow.")
        if chordRef is None:
            raise Error("'chordRef' must be specified in aeroProblem"
                        " for ADflow.")

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

        # Set mach defaults if user did not specified any machRef or machGrid values

        # If the user is running time spectral but did not specify
        # machGrid then default it to be equal to mach.

        if machRef is None:
            machRef = mach

        if machGrid is None:
            if self.getOption('equationMode').lower()=='time spectral':
                machGrid = mach
            else:
                # Steady, unsteady
                machGrid = 0.0

        # 1. Angle of attack:
        dToR = numpy.pi/180.0
        self.adflow.inputphysics.alpha = alpha*dToR
        self.adflow.inputphysics.beta = beta*dToR
        self.adflow.inputphysics.liftindex = liftIndex
        self.adflow.flowutils.adjustinflowangle()

        if self.getOption('printIterations') and self.comm.rank == 0:
            print('-> Alpha... %f '% numpy.real(alpha))

        # 2. Reference Points:
        self.adflow.inputphysics.pointref = [xRef, yRef, zRef]
        self.adflow.inputmotion.rotpoint = [xRot, yRot, zRot]
        self.adflow.inputphysics.pointrefec = [0.0, 0.0, 0.0]

        # 3. Reference Areas
        self.adflow.inputphysics.surfaceref = areaRef
        self.adflow.inputphysics.lengthref = chordRef

        # 4. Set mach numbers
        # Mach number for time spectral needs to be set to set to zero
        # If time-spectral (TS) then mach = 0, machcoef = mach, machgrid = mach
        # If Steady-State (SS), time-accurate (TA) then mach = mach, machcoef = mach, machgrid = 0
        if self.getOption('equationMode').lower()=='time spectral':
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
            self.adflow.updategamma() # NOTE! It is absolutely necessary to call this function, otherwise gamma is not properly updated.

        # 4. Periodic Parameters --- These are not checked/verified
        # and come directly from aeroProblem. Make sure you specify
        # them there properly!!
        if  self.getOption('alphaMode'):
            self.adflow.inputmotion.degreepolalpha = int(AP.degreePol)
            self.adflow.inputmotion.coefpolalpha = AP.coefPol
            self.adflow.inputmotion.omegafouralpha   = AP.omegaFourier
            self.adflow.inputmotion.degreefouralpha  = AP.degreeFourier
            self.adflow.inputmotion.coscoeffouralpha = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffouralpha = AP.sinCoefFourier
        elif  self.getOption('betaMode'):
            self.adflow.inputmotion.degreepolmach = int(AP.degreePol)
            self.adflow.inputmotion.coefpolmach = AP.coefPol
            self.adflow.inputmotion.omegafourbeta   = AP.omegaFourier
            self.adflow.inputmotion.degreefourbeta  = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourbeta = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourbeta = AP.sinCoefFourier
        elif self.getOption('machMode'):
            self.adflow.inputmotion.degreepolmach = int(AP.degreePol)
            self.adflow.inputmotion.coefpolmach = AP.coefPol
            self.adflow.inputmotion.omegafourmach   = AP.omegaFourier
            self.adflow.inputmotion.degreefourmach  = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourmach = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourmach = AP.sinCoefFourier
        elif  self.getOption('pMode'):
            ### add in lift axis dependence
            self.adflow.inputmotion.degreepolxrot = int(AP.degreePol)
            self.adflow.inputmotion.coefpolxrot = AP.coefPol
            self.adflow.inputmotion.omegafourxrot = AP.omegaFourier
            self.adflow.inputmotion.degreefourxrot  = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourxrot = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourxrot = AP.sinCoefFourier
        elif self.getOption('qMode'):
            self.adflow.inputmotion.degreepolzrot = int(AP.degreePol)
            self.adflow.inputmotion.coefpolzrot = AP.coefPol
            self.adflow.inputmotion.omegafourzrot = AP.omegaFourier
            self.adflow.inputmotion.degreefourzrot  = AP.degreeFourier
            self.adflow.inputmotion.coscoeffourzrot = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffourzrot = AP.sinCoefFourier
        elif self.getOption('rMode'):
            self.adflow.inputmotion.degreepolyrot = int(AP.degreePol)
            self.adflow.inputmotion.coefpolyrot = AP.coefPol
            self.adflow.inputmotion.omegafouryrot = AP.omegaFourier
            self.adflow.inputmotion.degreefouryrot  = AP.degreeFourier
            self.adflow.inputmotion.coscoeffouryrot = AP.cosCoefFourier
            self.adflow.inputmotion.sincoeffouryrot = AP.sinCoefFourier

        if not firstCall:
            self.adflow.initializeflow.updatebcdataalllevels()
            self.adflow.preprocessingapi.updateperiodicinfoalllevels()
            self.adflow.preprocessingapi.updategridvelocitiesalllevels()

    def getForces(self, groupName=None, TS=0):
        """ Return the forces on this processor on the families defined by groupName.

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
        self.adflow.getforces(numpy.ravel(forces), TS+1)
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
        self.adflow.getheatflux(fluxes, TS+1)
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
        fullTemp = self.adflow.gettnswall(npts, TS+1)

        # Now map new values in and set.
        fullTemp = self.mapVector(temperature, groupName, self.allWallsGroup, fullTemp)
        self.adflow.settnswall(fullTemp, TS+1)
        
    def getSurfacePoints(self, groupName=None, TS=0):

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

        # Get the required size
        npts, ncell = self._getSurfaceSize(groupName)
        pts = numpy.zeros((npts, 3), self.dtype)
        famList = self._getFamilyList(groupName)

        if npts > 0:
            self.adflow.surfaceutils.getsurfacepoints(pts.T, TS+1, famList)

        return pts

    def getSurfaceConnectivity(self, groupName=None, useBlanking=True):
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
        """

        if groupName is None:
            groupName = self.allWallsGroup

        # Set the list of surfaces this family requires
        famList = self._getFamilyList(groupName)
        npts, ncell = self._getSurfaceSize(groupName)
        conn =  numpy.zeros((ncell, 4), dtype='intc')
        self.adflow.surfaceutils.getsurfaceconnectivity(numpy.ravel(conn), famList, useBlanking)

        faceSizes = 4*numpy.ones(len(conn), 'intc')

        # Conver to 0-based ordering becuase we are in python
        return conn-1, faceSizes

    def setBCData(self, data, groupName=None, sps=1):
        """
        Take in the data for a bc and set it into the solver

        Data is a dictionary of the form:

        {'variable1':value1, 'variable2':value2 , ...}
        """

        if groupName is None:
            raise Error("Cannot set BCdata without specifying a group")

        nameList = []
        valList = []
        for k,v in iteritems(data):
            nameList.append(k)
            valList.append(v)

        nameArray = self._createFortranStringArray(nameList)
        bcInArray = numpy.array(valList, dtype=self.dtype)
        self.adflow.bcdata.setbcdata(nameArray,bcInArray, self._getFamilyList(groupName), sps)


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
        """Add a single desgin variable that ADflow knows about.

        Parameters
        ----------
        dv : str
            dv name. Must be in the self.possibleAeroDVs list
            """
        dv = dv.lower()
        if dv not in self.possibleAeroDVs:
            raise Error("%s was not one of the possible AeroDVs. "
                        "The complete list of DVs for ADflow is %s. "%(
                            dv, repr(set(self.possibleAeroDVs.keys()))))

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

    def solveAdjoint(self, aeroProblem, objective, forcePoints=None,
                      structAdjoint=None, groupName=None):

        # May be switching aeroProblems here
        self.setAeroProblem(aeroProblem)

        # Possibly setup adjoint matrices/vector as required
        self._setupAdjoint()

        # # Check to see if the RHS Partials have been computed
        if objective not in self.curAP.adflowData.adjointRHS:
            RHS = self.computeJacobianVectorProductBwd(
                funcsBar={objective.lower():1.0}, wDeriv=True)
            self.curAP.adflowData.adjointRHS[objective] = RHS.copy()
        else:
            RHS = self.curAP.adflowData.adjointRHS[objective].copy()

        # Check to see if we need to agument the RHS with a structural
        # adjoint:
        if structAdjoint is not None and groupName is not None:
            phi = self.mapVector(structAdjoint, groupName, self.allWallsGroup)

            agument = self.computeJacobianVectorProductBwd(
                fBar=phi, wDeriv=True)
            RHS -= agument

        # Check if objective is python 'allocated':
        if objective not in self.curAP.adflowData.adjoints:
            self.curAP.adflowData.adjoints[objective] = (
                numpy.zeros(self.getAdjointStateSize(), float))

        # Extract the psi:
        psi = self.curAP.adflowData.adjoints[objective]

        # Actually Solve the adjoint system...psi is updated with the
        # new solution.
        self.adflow.adjointapi.solveadjoint(RHS, psi, True)

        # Now set the flags and possibly reset adjoint
        if self.adflow.killsignals.adjointfailed:
            self.adjointFailed = True
            # Reset stored adjoint
            self.curAP.adflowData.adjoints[objective][:] = 0.0
        else:
            self.curAP.adflowData.adjoints[objective] = psi
            self.adjointFailed = False

    def _processAeroDerivatives(self, dIda):
        """This internal furncion is used to convert the raw array ouput from
        the matrix-free product bwd routine into the required
        dictionary format."""

	funcsSens = {}

        DVsRequired = list(self.curAP.DVNames.keys())
        for dv in DVsRequired:
            dvl = dv.lower()
            tmp = {}
            if dvl in ['altitude']:
                # This design variable is special. It combines changes
                # in temperature, pressure and density into a single
                # variable. Since we have derivatives for T, P and
                # rho, we simply chain rule it back to the the
                # altitude variable.
	        self.curAP.evalFunctionsSens(tmp, ['P', 'T', 'rho'])

                # Extract the derivatives wrt the independent
                # parameters in ADflow
                dIdP = dIda[self.possibleAeroDVs['p']]
                dIdT = dIda[self.possibleAeroDVs['t']]
                dIdrho = dIda[self.possibleAeroDVs['rho']]

                # Chain-rule to get the final derivative:
                funcsSens[self.curAP.DVNames[dv]] = (
                    tmp[self.curAP['P']][self.curAP.DVNames[dv]]*dIdP +
                    tmp[self.curAP['T']][self.curAP.DVNames[dv]]*dIdT +
                    tmp[self.curAP['rho']][self.curAP.DVNames[dv]]*dIdrho)
            elif dv.lower() in ['mach']:
                self.curAP.evalFunctionsSens(tmp, ['P', 'rho'])
                # Simular story for Mach: It is technically possible
                # to use a mach number for a fixed RE simulation. For
                # the RE to stay fixed and change the mach number, the
                # 'P' and 'rho' must also change. We have to chain run
                # this dependence back through to the final mach
                # derivative. When Mach number is used with altitude
                # or P and T, this calc is unnecessary, but won't do
                # any harm.
                dIdP = dIda[self.possibleAeroDVs['p']]
                dIdrho = dIda[self.possibleAeroDVs['rho']]

                # Chain-rule to get the final derivative:
                funcsSens[self.curAP.DVNames[dv]] = (
                    tmp[self.curAP['P']][self.curAP.DVNames[dv]]*dIdP +
                    tmp[self.curAP['rho']][self.curAP.DVNames[dv]]*dIdrho +
                    dIda[self.possibleAeroDVs['mach']])

            elif dvl in self.possibleAeroDVs:
                funcsSens[self.curAP.DVNames[dv]] = (
                    dIda[self.possibleAeroDVs[dvl]])
                if dvl == 'alpha':
                    funcsSens[self.curAP.DVNames[dv]] *= numpy.pi/180.0

	return funcsSens

    def _setAeroDVs(self):

        """ Do everything that is required to deal with aerodynamic
        design variables in ADflow"""

        DVsRequired = list(self.curAP.DVNames.keys())
        DVMap = {}
        for dv in DVsRequired:
            dv = dv.lower()
            if dv in ['altitude']:
                # All these variables need to be compined
                self._addAeroDV('P')
                self._addAeroDV('T')
                self._addAeroDV('rho')
            elif dv in ['mach']:
                self._addAeroDV('mach')
                self._addAeroDV('P')
                self._addAeroDV('rho')

            elif dv in self.possibleAeroDVs:
                self._addAeroDV(dv)
            else:
                raise Error("The design variable '%s' as specified in the"
                            " aeroProblem cannot be used with ADflow."% dv)

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
            relTol = self.getOption('adjointl2convergence')
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
            relTol = self.getOption('adjointl2convergence')
        outVec = self.adflow.solvedirectforrhs(inVec, relTol)

        return outVec

    def saveAdjointMatrix(self, baseFileName):
        """ Save the adjoint matrix to a binary petsc file for
        possible future external testing

        Parameters
        ----------
        basefileName : str
            Filename to use. The Adjoint matrix, PC matrix(if it exists)
            and RHS  will be written
        """
        adjointMatrixName = baseFileNmae + '_drdw.bin'
        pcMatrixName = baseFileNmae + '_drdwPre.bin'
        rhsName = baseFileName + '_rsh.bin'
        cellCenterName = baseFileName + '_cellCen.bin'
        self.adflow.adjointapi.saveadjointmatrix(adjointMatrixName)
        self.adflow.adjointapi.saveadjointpc(pcMatrixName)
        self.adflow.saveadjointrhs(rhsName)

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
            if self.getOption('equationMode').lower() == 'unsteady':
                self.adflow.preprocessingapi.shiftcoorandvolumes()
                self.adflow.aleutils.shiftlevelale()
            # Warp the mesh if surface coordinates are modified
            if warpMesh:
                timeA = time.time()
                self.mesh.warpMesh()
                newGrid = self.mesh.getSolverGrid()
                self.adflow.killsignals.routinefailed = False
                self.adflow.killsignals.fatalFail = False
                self.updateTime = time.time()-timeA
                if newGrid is not None:
                    self.adflow.warping.setgrid(newGrid)
            # Update geometric data, depending on the type of simulation
            if self.getOption('equationMode').lower() == 'unsteady':
                self.adflow.solvers.updateunsteadygeometry()
            else:
                self.adflow.preprocessingapi.updatecoordinatesalllevels()
                self.adflow.walldistance.updatewalldistancealllevels()
                self.adflow.preprocessingapi.updatemetricsalllevels()
                self.adflow.preprocessingapi.updategridvelocitiesalllevels()
            # Update flags
            self._updateGeomInfo = False
            self.adflow.killsignals.routinefailed = \
                self.comm.allreduce(
                bool(self.adflow.killsignals.routinefailed), op=MPI.LOR)
            self.adflow.killsignals.fatalfail = self.adflow.killsignals.routinefailed

    def getAdjointResNorms(self):
        '''
        Return the following adjoint residual norms:
        initRes Norm: Norm the adjoint RHS
        startRes Norm: Norm at the start of adjoint call (with possible non-zero restart)
        finalCFD Norm: Norm at the end of adjoint solve
        '''
        initRes  = self.adflow.adjointpetsc.adjresinit
        startRes = self.adflow.adjointpetsc.adjresstart
        finalRes = self.adflow.adjointpetsc.adjresfinal
        fail = self.adflow.killsignals.adjointfailed

        return initRes, startRes, finalRes, fail

    def getResNorms(self):
        """Return the initial, starting and final Res Norms. Typically
        used by an external solver."""
        return (numpy.real(self.adflow.iteration.totalr0),
                numpy.real(self.adflow.iteration.totalrstart),
                numpy.real(self.adflow.iteration.totalrfinal))

    def setResNorms(self, initNorm=None, startNorm=None, finalNorm=None):
        """ Set one of these norms if not None. Typlically used by an
        external solver"""
        if initNorm is not None:
            self.adflow.iteration.totalr0 = initNorm
        if startNorm is not None:
            self.adflow.iteration.totalrstart = startNorm
        if finalNorm is not None:
            self.adflow.iteration.finalNorm = finalNorm

    def _prescribedTSMotion(self):
        """Determine if we have prescribed motion timespectral analysis"""

        if (self.getOption('alphamode') or self.getOption('betamode') or
            self.getOption('machmode') or self.getOption('pmode') or
            self.getOption('qmode') or self.getOption('rmode') or
            self.getOption('altitudemode')):
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
        if objective.lower() in self.adflowCostFunctions.keys():
            return True
        else:
            return False

    # =========================================================================
    #   The following routines two routines are the workhorse of the
    #   forward and reverse mode routines. These can compute *ANY*
    #   product that is possible from the solver.
    #   =========================================================================

    def computeJacobianVectorProductFwd(self, xDvDot=None, xSDot=None, xVDot=None, wDot=None,
                                        residualDeriv=False, funcDeriv=False, fDeriv=False, 
                                        groupName=None):
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
        Fderiv : bool
            Flag specifiying if the derviative of the surface forces (tractions)
            should be returned
        groupName : str
            Optional group name to use for evaluating functions. Defaults to all 
            surfaces. 

        Returns
        -------
        dwdot, funcsdot, fDot : array, dict, array
            One or more of the these are return depending ont he *Deriv flags
        """

        if xDvDot is None and xSDot is None and xVDot is None and wDot is None:
            raise Error('computeJacobianVectorProductFwd: xDvDot, xSDot, xVDot and wDot cannot '
                        'all be None')

        self._setAeroDVs()
        nTime  = self.adflow.inputtimespectral.ntimeintervalsspectral

        # Default flags
        useState = False
        useSpatial = False

        # Process the Xs perturbation
        if xSDot is None:
            xsdot = numpy.zeros_like(self.coords0)
            xsdot = self.mapVector(xsdot, self.allFamilies,
                                   self.designFamilyGroup)
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
        if xDvDot is not None:
            useSpatial = True
            for key in xDvDot:
                val = xDvDot[key]
                if key.lower() == 'alpha':
                    val *= numpy.pi/180
                extradot[self.possibleAeroDVs[key.lower()]] = val
                    
        # For the geometric xDvDot perturbation we accumulate into the
        # already existing (and possibly nonzero) xsdot and xvdot
        if xDvDot is not None or xSDot is not None:
            if xDvDot is not None and self.DVGeo is not None:
                xsdot += self.DVGeo.totalSensitivityProd(xDvDot, self.curAP.ptSetName).reshape(xsdot.shape)
            if self.mesh is not None:
                xvdot += self.mesh.warpDerivFwd(xsdot)
            useSpatial = True

        # Sizes for output arrays
        costSize = self.adflow.costfunctions.ncostfunction
        fSize, nCell = self._getSurfaceSize(self.allWallsGroup)

        # Get the famList from the groupName
        famList = self._getFamilyList(groupName)
        
        dwdot,tmp,fdot = self.adflow.adjointapi.computematrixfreeproductfwd(
            xvdot, extradot, wdot, useSpatial, useState, costSize,  max(1, fSize), nTime, famList)

        # Explictly put fdot to nothing if size is zero
        if fSize==0:
            fdot = numpy.zeros((0, 3))

        # Process the derivative of the functions
        funcsDot = {}
        for f in self.curAP.evalFuncs:
            basicFunc = self.adflowCostFunctions[f.lower()][1]
            mapping = self.basicCostFunctions[basicFunc]
            funcsDot[f] = tmp[mapping - 1]

        # Assemble the returns
        returns = []
        if residualDeriv:
            returns.append(dwdot)
        if funcDeriv:
            returns.append(funcsDot)
        if fDeriv:
            returns.append(fdot.T)

        return tuple(returns) if len(returns) > 1 else returns[0]

    def computeJacobianVectorProductBwd(self, resBar=None, funcsBar=None, fBar=None,
                                    wDeriv=None, xVDeriv=None, xSDeriv=None,
                                    xDvDeriv=None, xDvDerivAero=None):
        """This the main python gateway for producing reverse mode jacobian
        vector products. It is not generally called by the user by
        rather internally or from another solver. A mesh object must
        be present for the xSDeriv=True flag and a mesh and DVGeo
        object must be present for xDvDeriv=True flag. Note that more
        than one of the specified return flags may be spcified. If
        more than one return is specified, the order of return is :
        (wDeriv, xVDeriv, XsDeriv, xDvDeriv, dXdvDerivAero).

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
            One or more of these are returned depending on the *Deriv flags provided.

        """
        # Error Checking
        if resBar is None and funcsBar is None and fBar is None:
            raise Error("computeJacobianVectorProductBwd: One of resBar, funcsBar and fBar"
                        " must be given. resBar=%s, funcsBar=%s, fBar=%s"% (
                            resBar, funcsBar, fBar))
        if (wDeriv is None and xVDeriv is None and xDvDeriv is None and
            xSDeriv is None and xDvDerivAero is None):
            raise Error("computeJacobianVectorProductBwd: One of wDeriv, xVDeriv, "
                        "xDvDeriv and xDvDerivAero must be given as True. "
                        "wDeriv=%s, xVDeriv=%s, xDvDeriv=%s, xSDeriv=%s xDvDerivAero=%s."%(
                            wDeriv, xVDeriv, xDvDeriv, xSDeriv, xDvDerivAero))
        self._setAeroDVs()

        # ---------------------
        #  Check for resBar
        # ---------------------
        if resBar is None:
            if self.getOption('frozenTurbulence'):
                resBar = numpy.zeros(self.getAdjointStateSize())
            else:
                resBar = numpy.zeros(self.getStateSize())

        # -------------------------
        #  Check for fBar (forces)
        # ------------------------
        nTime  = self.adflow.inputtimespectral.ntimeintervalsspectral
        nPts, nCell = self._getSurfaceSize(self.allWallsGroup)

        if fBar is None:
            fBar = numpy.zeros((nTime, nPts, 3))
        else:
            # Expand out to the sps direction in case there were only
            # 2 dimensions.
            fBar= fBar.reshape((nTime, nPts, 3))

        # ---------------------
        #  Check for funcsBar
        # ---------------------

        if funcsBar is None:
            funcsBar = numpy.zeros(self.adflow.costfunctions.ncostfunction)
            famList = self._getFamilyList(self.allWallsGroup)
        else:
            tmp = numpy.zeros(self.adflow.costfunctions.ncostfunction)

            # We have to make sure that the user has supplied
            # functions that have the same group, otherwise, we can't
            # do it
            groups = set()

            for f in funcsBar:
                if f.lower() in self.adflowCostFunctions:

                    groupName = self.adflowCostFunctions[f][0]
                    basicFunc = self.adflowCostFunctions[f][1]
                    groups.add(groupName)

                    mapping = self.basicCostFunctions[basicFunc]
                    tmp[mapping-1] = funcsBar[f]

            if len(groups) == 1:
                # We're ok..there was only one group from the
                # functions..
                famList = self._getFamilyList(list(groups)[0])
            elif len(groups)==0:
                # this is probably an aerostructural function
                # take the self.allFamilies group
                famList = self._getFamilyList(self.allFamilies)

            elif len(groups) > 1:
                raise Error("Attemping to compute a jacobian vector product "
                            "with multiple functions that have different "
                            "groups. This is not allowed.")

            # If there wasn't any actual funcsBar, then tmp is still
            # just zeros
            funcsBar = tmp

        # Determine if we can do a cheaper state-variable only
        # computation or if we need to include the spatial terms:
        useSpatial = False
        useState = False
        if wDeriv:
            useState = True
        if xDvDeriv or xVDeriv or xSDeriv or xDvDerivAero:
            useSpatial = True

        # Do actual Fortran call.
        xvbar, extrabar, wbar = self.adflow.adjointapi.computematrixfreeproductbwd(
            resBar, funcsBar, fBar.T, useSpatial, useState, self.getSpatialSize(),
            self.adflow.adjointvars.ndesignextra, famList)

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
                xsbar = self.mapVector(xsbar, self.meshFamilyGroup,
                                       self.designFamilyGroup)

                if xSDeriv:
                    returns.append(xsbar)
            else:
                # Only an error if xSDeriv is requested...since we
                # can't do it. xDVDeriv may be specified even when no
                # mesh is present.
                if xSDeriv:
                    raise Error("Could not complete requested xSDeriv "
                                "derivatives since no mesh is present")

            # Process all the way back to the DVs:
            if xDvDeriv:
                xdvbar = {}
                if self.mesh is not None: # Include geometric
                                          # derivatives if mesh is
                                          # present
                    if self.DVGeo is not None and self.DVGeo.getNDV() > 0:
                        xdvbar.update(self.DVGeo.totalSensitivity(
                            xsbar, self.curAP.ptSetName, self.comm, config=self.curAP.name))
                    else:
                        if self.comm.rank == 0:
                            ADFLOWWarning("No DVGeo object is present or it has no "
                                        "design variables specified. No geometric "
                                        "derivatives computed.")
                else:
                    if self.comm.rank == 0:
                        ADFLOWWarning("No mesh object is present. No geometric "
                                    "derivatives computed.")

                # Include aero derivatives here:
                xdvbar.update(self._processAeroDerivatives(extrabar))
                returns.append(xdvbar)


        # Include the aerodynamic variables if requested to do so.
        if xDvDerivAero:
            xdvaerobar = {}
            xdvaerobar.update(self._processAeroDerivatives(extrabar))
            returns.append(xdvaerobar)

        # Single return (most frequent) is 'clean', otherwise a tuple.
        return tuple(returns) if len(returns) > 1 else returns[0]

    def mapVector(self, vec1, groupName1, groupName2, vec2=None):
        """This is the main workhorse routine of everything that deals with
        families in ADflow. The purpose of this routine is to convert a
        vector 'vec1' (of size Nx3) that was evaluated with
        'groupName1' and expand or contract it (and adjust the
        ordering) to produce 'vec2' evaluated on groupName2.

        A little ascii art might help. Consider the following "mesh"
        . Family 'fam1' has 9 points, 'fam2' has 10 pts and 'fam3' has
        5 points.  Consider that we have also also added two
        additional groups: 'f12' containing 'fam1' and 'fma2' and a
        group 'f23' that contains families 'fam2' and 'fam3'. The vector
        we want to map is 'vec1'. It is length 9+10. All the 'x's are
        significant values.

        The call: mapVector(vec1, 'f12', 'f23')

        will produce the "returned vec" array, containing the
        significant values from 'fam2', where the two groups overlap,
        and the new values from 'fam3' set to zero. The values from
        fam1 are lost. The returned vec has size 15.

            fam1     fam2      fam3
        |---------+----------+------|

        |xxxxxxxxx xxxxxxxxxx|        <- vec1
                  |xxxxxxxxxx 000000| <- returned vec (vec2)

        It is also possible to pass in vec2 into this routine. For
        that case, the existing values in the array will not be
        kept. In the previous examples, the values cooresponding to
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
            raise Error("'%s' or '%s' is not a family in the CGNS file or has not been added"
                        " as a combination of families"%(groupName1, groupName2))

        # Shortcut:
        if groupName1 == groupName2:
            return vec1

        if vec2 is None:
            npts, ncell = self._getSurfaceSize(groupName2)
            vec2 = numpy.zeros((npts, 3), self.dtype)

        famList1 = self.families[groupName1]
        famList2 = self.families[groupName2]
        self.adflow.surfaceutils.mapvector(vec1.T, famList1, vec2.T, famList2)

        return vec2

    def getStateSize(self):
        """Return the number of degrees of freedom (states) that are on this
        processor. This is (number of states)*(number of
        cells)*(number of spectral instances)

        """

        nstate = self.adflow.flowvarrefstate.nw
        ncells = self.adflow.adjointvars.ncellslocal[0]
        ntime  = self.adflow.inputtimespectral.ntimeintervalsspectral

        return nstate*ncells*ntime

    def getAdjointStateSize(self):
        """Return the number of ADJOINT degrees of freedom (states)
        that are on this processor. The reason this is different from
        getStateSize() is that if frozenTurbulence is used for RANS,
        the nonlinear system has 5+neq turb states per cell, while the
        adjoint still has 5."""
        if self.getOption('frozenTurbulence'):
            nstate = self.adflow.flowvarrefstate.nwf
        else:
            nstate = self.adflow.flowvarrefstate.nw

        ncells = self.adflow.adjointvars.ncellslocal[0]
        ntime  = self.adflow.inputtimespectral.ntimeintervalsspectral

        return nstate*ncells*ntime

    def getSpatialSize(self):
        """Return the number of degrees of spatial degrees of freedom
        on this processor. This is (number of nodes)*(number of
        spectral instances)*3"""

        nnodes = self.adflow.adjointvars.nnodeslocal[0]
        ntime  = self.adflow.inputtimespectral.ntimeintervalsspectral

        return 3*nnodes*ntime

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
        nPts, nCell = self._getSurfaceSize(self.allWallsGroup)
        xRand = self.getSpatialPerturbation(seed)
        surfRand = numpy.zeros((nPts, 3))
        famList = self._getFamilyList(self.allWallsGroup)
        self.adflow.warping.getsurfaceperturbation(xRand, numpy.ravel(surfRand), famList)
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

        return self.adflow.warping.getstateperturbation(
            randVec, self.getStateSize())
        
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
        offsets = numpy.zeros(len(ptSizes), 'intc')
        offsets[1:] = numpy.cumsum(ptSizes)[:-1]
        
        # Now Re-assemble the global CGNS vector. 
        nDOFCGNS = numpy.max(cgnsIndices) +1 
        CGNSVec = numpy.zeros(nDOFCGNS)
        CGNSVec[cgnsIndices] = allPts
        CGNSVec = CGNSVec.reshape((len(CGNSVec)/3, 3))
        
        # Run the pointReduce on the CGNS nodes
        uniquePts, linkTmp, nUnique = self.adflow.utils.pointreduce(CGNSVec.T, 1e-12)
        
        # Expand link out to the 3x the size and convert to 1 based ordering
        linkTmp -= 1
        link = numpy.zeros(len(linkTmp)*3, 'intc')
        link[0::3] = 3*linkTmp
        link[1::3] = 3*linkTmp+1
        link[2::3] = 3*linkTmp+2

        # Set the seed and everyone produces the random vector for
        # nUnique pts.
        numpy.random.seed(seed)
        randVec = numpy.random.random(nUnique*3)

        # Finally just extract out our own part:
        iStart = offsets[self.comm.rank]
        iEnd   = iStart + self.getSpatialSize()
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
        dXv    = numpy.hstack(self.comm.allgather(dXv))
        norm = None
        if self.myid == 0:
            allPts = allPts.reshape((len(allPts)/3, 3))
            dXv = dXv.reshape((len(dXv)/3,3))
            # Run the pointReduce on all nodes
            uniquePts, link, nUnique = self.adflow.utils.pointreduce(allPts.T, 1e-12)
            uniquePtsBar = numpy.zeros((nUnique,3))
            link = link -1 # Convert to zero-based for python:
            
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
        """ Return the adjoint values for objective if they
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
        
    def getFreeStreamResidual(self, aeroProblem):
        self.setAeroProblem(aeroProblem)
        rhoRes, totalRRes = self.adflow.nksolver.getfreestreamresidual()
        return totalRRes

    def _getSurfaceSize(self, groupName, useBlanking=True):
        """Internal routine to return the size of a particular surface. This
        does *NOT* set the actual family group"""
        if groupName is None:
            groupName = self.allFamilies

        if groupName not in self.families:
            raise Error("'%s' is not a family in the CGNS file or has not been added"
                        " as a combination of families"%groupName)

        [nPts, nCells] = self.adflow.surfaceutils.getsurfacesize(
            self.families[groupName], useBlanking)
        return nPts, nCells

    def setOption(self, name, value):
        """
        Set Solver Option Value
        """
        name = name.lower()

        # Make sure we are not trying to change an immutable option if
        # we are not allowed to.
        if self.solverCreated and name in self.imOptions:
            raise Error("Option '%-35s' cannot be modified after the solver "
                        "is created."%name)

        # Check to see if we have a deprecated option. Print a useful
        # warning that this is deprecated.
        if name in self.deprecatedOptions:
            if self.comm.rank == 0:
                ADFLOWWarning("Option '%-29s\' is a deprecated ADflow Option |"% name)
            return

        # Try the option in the option dictionary to make sure we are setting a valid option
        if name not in self.defaultOptions:
            if self.comm.rank == 0:
                ADFLOWWarning("Option '%-30s' is not a valid ADflow Option |"%name)
            return

        # Now we know the option exists, lets check if the type is ok:
        if isinstance(value, self.defaultOptions[name][0]):
            self.options[name] = [type(value),value]
        else:
            raise Error("Datatype for Option %-35s was not valid \n "
                        "Expected data type is %-47s \n "
                        "Received data type is %-47s"% (
                            name, self.defaultOptions[name][0], type(value)))

        # If the option is in the ignoredOption list, we just return.
        if name in self.ignoreOptions:
            return

        # Do special Options individually
        if name in self.specialOptions:
            if name in ['monitorvariables',
                        'surfacevariables',
                        'volumevariables',
                        'isovariables']:
                varStr = ''
                for i in xrange(len(value)):
                    varStr = varStr + value[i] + '_'
                # end if
                varStr = varStr[0:-1] # Get rid of last '_'
                if name == 'monitorvariables':
                    self.adflow.inputparamroutines.monitorvariables(varStr)
                if name == 'surfacevariables':
                    self.adflow.inputparamroutines.surfacevariables(varStr)
                if name == 'volumevariables':
                    self.adflow.inputparamroutines.volumevariables(varStr)
                if name == 'isovariables':
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
                            raise Error("Option 'restartfile' string cannot be empty. "
                                        "If not performing a restart, 'restartFile' "
                                        "option must be set to None")
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
                                    self.adflow.initializeflow.setrestartfiles(val,i+1)
                            else:
                                raise Error("Datatype for Option %-35s was not "
                                            "valid. Expected list of <type 'str'>. "
                                            "Received data type is "
                                            "%-47s"% (name, type(value[0])))
                        else:
                            raise Error("Option %-35s of %-35s contains %-35s "
                                        "elements. Must contain at least 1 "
                                        "restart file of type <type "
                                        "'str'>"% (name, type(value), len(value)))
                    else:
                        raise Error("Datatype for Option %-35s not valid. "
                                    "Expected data type is <type 'str'> or <type "
                                    "'list'>. Received data type is %-47s"% (name, type(value)))

            elif name == 'isosurface':
                # We have a bit of work to do...extract out the
                # names, and there can be more than 1 value per variables
                var = []
                val = []
                isoDict = value
                for key in isoDict.keys():

                    isoVals = numpy.atleast_1d(isoDict[key])
                    for i in xrange(len(isoVals)):
                        var.append(key)
                        val.append(isoVals[i])

                val = numpy.array(val)

                self.adflow.inputparamroutines.initializeisosurfacevariables(val)
                for i in xrange(len(val)):
                    self.adflow.inputparamroutines.setisosurfacevariable(var[i], i+1)

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
                        if  1 <= len(value) and len(value) <= 4:
                            if type(value[0]) is float:
                                tmp_turbresscalar[0:len(value)] = value[:]
                            else:
                                raise Error("Datatype for Option %-35s was not "
                                            "valid. Expected list of "
                                            "<type 'float'>. Received data type "
                                            "is %-47s"% (name, type(value[0])))
                        else:
                            raise Error("Option %-35s of %-35s contains %-35s "
                                        "elements. Min and max number of "
                                        "elements are 1 and 4 "
                                        "respectively"% (name, type(value), len(value)))
                    else:
                        raise Error("Datatype for Option %-35s not valid. Expected "
                                    "data type is <type 'float'> or <type "
                                    "'list'>. Received data type is %-47s"% (name, type(value)))

                    module = self.moduleMap[self.optionMap[name][0]]
                    variable = self.optionMap[name][1]
                    setattr(module, variable, tmp_turbresscalar)

            # Special option has been set so return from function
            return

        # All other options do genericaly by setting value in module:
        # Check if there is an additional mapping to what actually
        # has to be set in the solver

        if isinstance(self.optionMap[name], dict):
            module = self.moduleMap[self.optionMap[name]['location'][0]]
            variable = self.optionMap[name]['location'][1]
            value = self.optionMap[name][value.lower()]
        else:
            module = self.moduleMap[self.optionMap[name][0]]
            variable = self.optionMap[name][1]

        # If the value is a string, pads additional spaces
        if isinstance(value, str):
            spacesToAdd = self.adflow.constants.maxstringlen - len(value)
            value = ''.join([value,' '*spacesToAdd])

        # Set in the correct module
        setattr(module, variable, value)

    def getOption(self, name):
        # Redefine the getOption def from the base class so we can
        # make sure the name is lowercase

        if name.lower() in self.defaultOptions:
            return self.options[name.lower()][1]
        else:
            raise Error('%s is not a valid option name.'% name)

    def _getDefOptions(self):
        """
        There are many options for ADflow. These technically belong in
        the __init__ function but it gets far too long so we split
        them out.
        """
        defOpts = {
            # Input file parameters
            'gridfile':[str, 'default.cgns'],
            'restartfile':[object, None],

            # Surface definition parameters:
            'meshsurfacefamily':[object, None],
            'designsurfacefamily':[object, None],

            # Output Parameters
            'storerindlayer':[bool, True],
            'outputdirectory':[str, './'],
            'outputsurfacefamily':[str, 'allSurfaces'],
            'writesurfacesolution':[bool,True],
            'writevolumesolution':[bool,True],
            'writetecplotsurfacesolution':[bool,False],
            'nsavevolume':[int,1],
            'nsavesurface':[int,1],
            'solutionprecision':[str,'single'],
            'gridprecision':[str,'double'],
            'solutionprecisionsurface':[str,'single'],
            'gridprecisionsurface':[str,'single'],
            'isosurface':[dict, {}],
            'isovariables':[list, []],
            'viscoussurfacevelocities':[bool, True],

            # Physics Paramters
            'discretization':[str, 'central plus scalar dissipation'],
            'coarsediscretization':[str, 'central plus scalar dissipation'],
            'limiter':[str, 'vanalbeda'],
            'smoother':[str, 'runge kutta'],
            'equationtype': [str, 'euler'],
            'equationmode': [str, 'steady'],
            'flowtype':[str, 'external'],
            'turbulencemodel':[str, 'sa'],
            'turbulenceorder':[str, 'first order'],
            'turbresscale':[object, None],
            'turbulenceproduction':[str, 'strain'],
            'useqcr':[bool, False],
            'userotationsa':[bool, False],
            'useft2sa':[bool, True],
            'eddyvisinfratio':[float, .009],
            'usewallfunctions':[bool, False],
            'useapproxwalldistance':[bool, True],
            'eulerwalltreatment':[str, 'linear pressure extrapolation'],
            'viscwalltreatment':[str, 'constant pressure extrapolation'],
            'dissipationscalingexponent':[float, 0.67],
            'vis4':[float, 0.0156],
            'vis2':[float, 0.25],
            'vis2coarse':[float, 0.5],
            'restrictionrelaxation':[float, .80],
            'liftindex':[int, 2],
            'lowspeedpreconditioner':[bool, False],
            'walldistcutoff':[float, 1e20],

            # Common Paramters
            'ncycles':[int, 500],
            'ncyclescoarse':[int, 500],
            'nsubiterturb':[int, 1],
            'nsubiter':[int, 1],
            'cfl':[float, 1.7],
            'cflcoarse':[float, 1.0],
            'mgcycle':[str, '3w'],
            'mgstartlevel':[int, -1],
            'resaveraging':[str,'alternateresaveraging'],
            'smoothparameter':[float, 1.5],
            'cfllimit':[float, 1.5],

            # Overset Parameters:
            'nearwalldist':[float, 0.1],
            'backgroundvolscale':[float, 1.0],
            'oversetprojtol':[float, 1e-12],
            'overlapfactor':[float, 0.9],
            'debugzipper':[bool, False],
            'zippersurfacefamily':[object, None],
            'cutcallback':[object, None],

            # Unsteady Paramters
            'timeintegrationscheme':[str, 'bdf'],
            'timeaccuracy':[int, 2],
            'ntimestepscoarse':[int, 48],
            'ntimestepsfine':[int, 400],
            'deltat':[float, .010],
            'useale':[bool, True],
            'usegridmotion':[bool, False],
            'coupledsolution':[bool, False],

            # Time Spectral Paramters
            'timeintervals': [int, 1],
            'alphamode':[bool, False],
            'betamode':[bool, False],
            'machmode':[bool, False],
            'pmode':[bool, False],
            'qmode':[bool, False],
            'rmode':[bool, False],
            'altitudemode':[bool, False],
            'windaxis':[bool, False],
            'alphafollowing':[bool,True],
            'tsstability': [bool, False],

            # Convergence Paramters
            'l2convergence':[float, 1e-6],
            'l2convergencerel':[float, 1e-16],
            'l2convergencecoarse':[float, 1e-2],
            'maxl2deviationfactor':[float, 1.0],

            # Newton-Krylov Parameters
            'usenksolver':[bool, False],
            'nkswitchtol':[float, 2.5e-4],
            'nksubspacesize':[int, 60],
            'nklinearsolvetol':[float, 0.3],
            'nkuseew':[bool, True],
            'nkadpc':[bool, False],
            'nkviscpc':[bool, False],
            'nkasmoverlap':[int, 1],
            'nkpcilufill':[int, 2],
            'nkjacobianlag':[int, 20],
            'applypcsubspacesize':[int, 10],
            'nkinnerpreconits':[int, 1],
            'nkouterpreconits':[int, 1],
            'nkcfl0':[float, 100.0],
            'nkls':[str, 'cubic'],
            'rkreset':[bool, False],
            'nrkreset':[int, 5],

            # Approximate Newton-Krylov Parameters
            'useanksolver':[bool, False],
            'ankuseturbdadi':[bool, True],
            'ankswitchtol':[float, 1e-2],
            'anksubspacesize':[int, 5],
            'anklinearsolvetol':[float, 0.5],
            'ankasmoverlap':[int, 1],
            'ankpcilufill':[int, 1],
            'ankjacobianlag':[int, 20],
            'ankinnerpreconits':[int, 1],
            'ankcfl0':[float, 1.0],

            # Load Balance/partitioning parameters
            'blocksplitting':[bool, True],
            'loadimbalance':[float, 0.1],
            'loadbalanceiter':[int, 10],
            'partitiononly':[bool, False],
            'partitionlikenproc':[int, -1],

            # Misc Paramters
            'autosolveretry':[bool, False],
            'autoadjointretry':[bool, False],
            'numbersolutions':[bool, True],
            'printiterations':[bool, True],
            'printtiming':[bool, True],
            'setmonitor':[bool, True],
            'printwarnings':[bool, True],
            'monitorvariables':[list, ['cpu','resrho', 'resturb', 'cl', 'cd']],
            'surfacevariables':[list, ['cp','vx', 'vy','vz', 'mach']],
            'volumevariables':[list, ['resrho']],

            # Multidisciplinary Coupling Parameters:
            'forcesastractions':[bool, True],

            # Adjoint Paramters
            'adjointl2convergence':[float, 1e-6],
            'adjointl2convergencerel':[float, 1e-16],
            'adjointl2convergenceabs':[float, 1e-16],
            'adjointdivtol':[float, 1e5],
            'approxpc': [bool, True],
            'adpc': [bool, False],
            'viscpc':[bool,False],
            'usediagtspc':[bool, True],
            'restartadjoint':[bool, True],
            'adjointsolver': [str, 'gmres'],
            'adjointmaxiter': [int, 500],
            'adjointsubspacesize' : [int, 100],
            'adjointmonitorstep': [int, 10],
            'dissipationlumpingparameter':[float, 6.0],
            'preconditionerside': [str, 'right'],
            'matrixordering': [str, 'rcm'],
            'globalpreconditioner': [str, 'additive schwartz'],
            'localpreconditioner' : [str, 'ilu'],
            'ilufill': [int, 2],
            'asmoverlap' : [int, 1],
            'innerpreconits':[int, 1],
            'outerpreconits':[int, 3],
            'applyadjointpcsubspacesize':[int, 20],
            'frozenturbulence':[bool, True],
            'usematrixfreedrdw':[bool, True],

            # ADjoint debugger
            'firstrun':[bool, True],
            'verifystate':[bool, True],
            'verifyspatial':[bool, True],
            'verifyextra':[bool, True],

            # Function parmeters
            'sepsensoroffset':[float, 0.0],
            'sepsensorsharpness':[float, 10.0],
            }

        return defOpts

    def _getImmutableOptions(self):
        """We define the list of options that *cannot* be changed after the
        object is created. ADflow will raise an error if a user tries to
        change these. The strings for these options are placed in a set"""

        return ('gridfile', 'equationtype', 'equationmode', 'flowtype',
                'useapproxwalldistance', 'liftindex', 'mgcycle',
                'mgstartlevel', 'timeintegrationscheme', 'timeaccuracy',
                'useale', 'timeintervals', 'blocksplitting',
                'loadimbalance', 'loadbalanceiter', 'partitiononly',
                'meshSurfaceFamily', 'designSurfaceFamily', 
                'zippersurfacefamily', 'cutcallback')

    def _getOptionMap(self):
        """ The ADflow option map and module mapping"""

        moduleMap = {'io': self.adflow.inputio,
                     'discr':self.adflow.inputdiscretization,
                     'iter':self.adflow.inputiteration,
                     'physics':self.adflow.inputphysics,
                     'stab': self.adflow.inputtsstabderiv,
                     'nk': self.adflow.nksolver,
                     'ank': self.adflow.anksolver,
                     'adjoint': self.adflow.inputadjoint,
                     'cost': self.adflow.costfunctions,
                     'unsteady':self.adflow.inputunsteady,
                     'motion':self.adflow.inputmotion,
                     'parallel':self.adflow.inputparallel,
                     'ts':self.adflow.inputtimespectral,
                     'overset':self.adflow.inputoverset,
                 }

        # In the option map, we first list the "module" defined in
        # module map, and "variable" the variable to set in that module.

        optionMap = {
            # Common Paramters
            'gridfile':['io', 'gridfile'],
            'storerindlayer':['io', 'storerindlayer'],
            'nsavevolume':['io', 'nsavevolume'],
            'nsavesurface':['iter', 'nsavesurface'],
            'viscoussurfacevelocities':['io', 'viscoussurfacevelocities'],
            'solutionprecision':{'single':self.adflow.constants.precisionsingle,
                                 'double':self.adflow.constants.precisiondouble,
                                 'location':['io', 'precisionsol']},
            'gridprecision':{'single':self.adflow.constants.precisionsingle,
                             'double':self.adflow.constants.precisiondouble,
                             'location':['io', 'precisiongrid']},
            'solutionprecisionsurface':{'single':self.adflow.constants.precisionsingle,
                                        'double':self.adflow.constants.precisiondouble,
                                        'location':['io', 'precisionsol']},
            'gridprecisionsurface':{'single':self.adflow.constants.precisionsingle,
                                    'double':self.adflow.constants.precisiondouble,
                                    'location':['io', 'precisiongrid']},
            # Physics Paramters
            'discretization':{'central plus scalar dissipation': self.adflow.constants.dissscalar,
                              'central plus matrix dissipation': self.adflow.constants.dissmatrix,
                              'central plus cusp dissipation':self.adflow.constants.disscusp,
                              'upwind':self.adflow.constants.upwind,
                              'location':['discr', 'spacediscr']},
            'coarsediscretization':{'central plus scalar dissipation': self.adflow.constants.dissscalar,
                                    'central plus matrix dissipation': self.adflow.constants.dissmatrix,
                                    'central plus cusp dissipation': self.adflow.constants.disscusp,
                                    'upwind': self.adflow.constants.upwind,
                                    'location':['discr', 'spacediscrcoarse']},
            'limiter':{'vanalbeda':self.adflow.constants.vanalbeda,
                       'minmod':self.adflow.constants.minmod,
                       'nolimiter':self.adflow.constants.nolimiter,
                       'location':['discr', 'limiter']},
            'smoother':{'runge kutta':self.adflow.constants.rungekutta,
                        'lu sgs':self.adflow.constants.nllusgs,
                        'lu sgs line':self.adflow.constants.nllusgsline,
                        'dadi':self.adflow.constants.dadi,
                        'location':['iter', 'smoother']},

            'equationtype':{'euler':self.adflow.constants.eulerequations,
                            'laminar ns':self.adflow.constants.nsequations,
                            'rans':self.adflow.constants.ransequations,
                            'location':['physics', 'equations']},
            'equationmode':{'steady':self.adflow.constants.steady,
                            'unsteady':self.adflow.constants.unsteady,
                            'time spectral':self.adflow.constants.timespectral,
                            'location':['physics', 'equationmode']},
            'flowtype':{'internal':self.adflow.constants.internalflow,
                        'external':self.adflow.constants.externalflow,
                        'location':['physics', 'flowtype']},
            'turbulencemodel':{'sa':self.adflow.constants.spalartallmaras,
                               'sae':self.adflow.constants.spalartallmarasedwards,
                               'k omega wilcox':self.adflow.constants.komegawilcox,
                               'k omega modified':self.adflow.constants.komegamodified,
                               'ktau':self.adflow.constants.ktau,
                               'menter sst':self.adflow.constants.mentersst,
                               'v2f':self.adflow.constants.v2f,
                               'location':['physics', 'turbmodel']},
            'turbulenceorder':{'first order':1,
                               'second order':2,
                               'location':['discr', 'orderturb']},
            'turbresscale':['iter', 'turbresscale'],
            'turbulenceproduction':{'strain':self.adflow.constants.strain,
                                    'vorticity':self.adflow.constants.vorticity,
                                    'katolaunder':self.adflow.constants.katolaunder,
                                    'location':['physics', 'turbprod']},
            'useqcr':['physics', 'useqcr'],
            'userotationsa':['physics', 'userotationsa'],
            'useft2sa':['physics', 'useft2sa'],
            'eddyvisinfratio':['physics', 'eddyvisinfratio'],
            'usewallfunctions':['physics', 'wallfunctions'],
            'walldistcutoff':['physics', 'walldistcutoff'],
            'useapproxwalldistance':['discr', 'useapproxwalldistance'],
            'eulerwalltreatment':{'linear pressure extrapolation':self.adflow.constants.linextrapolpressure,
                                  'constant pressure extrapolation':self.adflow.constants.constantpressure,
                                  'quadratic pressure extrapolation':self.adflow.constants.quadextrapolpressure,
                                  'normal momentum':self.adflow.constants.normalmomentum,
                                  'location':['discr', 'eulerwallbctreatment']},
            'viscwalltreatment':{'linear pressure extrapolation':self.adflow.constants.linextrapolpressure,
                                 'constant pressure extrapolation':self.adflow.constants.constantpressure,
                                 'location':['discr', 'viscwallbctreatment']},
            'dissipationscalingexponent':['discr', 'adis'],
            'vis4':['discr', 'vis4'],
            'vis2':['discr', 'vis2'],
            'vis2coarse':['discr', 'vis2coarse'],
            'restrictionrelaxation':['iter', 'fcoll'],
            'forcesastractions':['physics', 'forcesastractions'],
            'lowspeedpreconditioner':['discr', 'lowspeedpreconditioner'],

            # Common Paramters
            'ncycles':['iter', 'ncycles'],
            'ncyclescoarse':['iter', 'ncyclescoarse'],
            'nsubiterturb':['iter', 'nsubiterturb'],
            'nsubiter':['iter', 'nsubiterations'],
            'cfl':['iter', 'cfl'],
            'cflcoarse':['iter', 'cflcoarse'],
            'mgcycle':['iter', 'mgdescription'],
            'mgstartlevel':['iter', 'mgstartlevel'],
            'resaveraging':{'noresaveraging':self.adflow.constants.noresaveraging,
                            'alwaysresaveraging':self.adflow.constants.alwaysresaveraging,
                            'alternateresaveraging':self.adflow.constants.alternateresaveraging,
                            'location':['iter', 'resaveraging']},
            'smoothparameter':['iter', 'smoop'],
            'cfllimit':['iter', 'cfllimit'],

            # Overset Parameters
            'nearwalldist':['overset','nearwalldist'],
            'backgroundvolscale':['overset','backgroundvolscale'],
            'oversetprojtol':['overset','oversetprojtol'],
            'overlapfactor':['overset','overlapfactor'],
            'debugzipper':['overset','debugzipper'],

            # Unsteady Params
            'timeintegrationscheme':{'bdf':self.adflow.constants.bdf,
                                     'explicitrk':self.adflow.constants.explicitrk,
                                     'implicitrk':self.adflow.constants.implicitrk,
                                     'location':['unsteady', 'timeintegrationscheme']},
            'timeaccuracy':['unsteady', 'timeaccuracy'],
            'ntimestepscoarse':['unsteady', 'ntimestepscoarse'],
            'ntimestepsfine':['unsteady', 'ntimestepsfine'],
            'deltat':['unsteady', 'deltat'],
            'useale':['unsteady', 'useale'],

            # Grid motion Params
            'usegridmotion':['motion', 'gridmotionspecified'],

            # Time Spectral Paramters
            'timeintervals':['ts', 'ntimeintervalsspectral'],
            'alphamode':['stab', 'tsalphamode'],
            'betamode':['stab', 'tsbetamode'],
            'machmode':['stab', 'tsmachmode'],
            'pmode':['stab', 'tspmode'],
            'qmode':['stab', 'tsqmode'],
            'rmode':['stab', 'tsrmode'],
            'altitudemode':['stab', 'tsaltitudemode'],
            'windaxis':['stab', 'usewindaxis'],
            'alphafollowing':['stab', 'tsalphafollowing'],
            'tsstability':['stab', 'tsstability'],

            # Convergence Paramters
            'l2convergence':['iter', 'l2conv'],
            'l2convergencerel':['iter', 'l2convrel'],
            'l2convergencecoarse':['iter', 'l2convcoarse'],
            'maxl2deviationfactor':['iter', 'maxl2deviationfactor'],

            # Newton-Krylov Paramters
            'usenksolver':['nk', 'usenksolver'],
            'nkuseew':['nk', 'nk_useew'],
            'nkswitchtol':['nk', 'nk_switchtol'],
            'nksubspacesize':['nk', 'nk_subspace'],
            'nklinearsolvetol':['nk', 'nk_rtolinit'],
            'nkasmoverlap':['nk', 'nk_asmoverlap'],
            'nkpcilufill':['nk', 'nk_ilufill'],
            'nkjacobianlag':['nk', 'nk_jacobianlag'],
            'nkadpc':['nk', 'nk_adpc'],
            'nkviscpc':['nk', 'nk_viscpc'],
            'applypcsubspacesize':['nk', 'applypcsubspacesize'],
            'nkinnerpreconits':['nk', 'nk_innerpreconits'],
            'nkouterpreconits':['nk', 'nk_outerpreconits'],
            'nkcfl0':['nk', 'nk_cfl0'],
            'nkls':{'none':self.adflow.constants.nolinesearch,
                    'cubic':self.adflow.constants.cubiclinesearch,
                    'non monotone':self.adflow.constants.nonmonotonelinesearch,
                    'location':['nk', 'nk_ls']},
            'rkreset':['iter', 'rkreset'],
            'nrkreset':['iter', 'miniternum'],

            # Approximate Newton-Krylov Paramters
            'useanksolver':['ank', 'useanksolver'],
            'ankuseturbdadi':['ank', 'ank_useturbdadi'],
            'ankswitchtol':['ank', 'ank_switchtol'],
            'anksubspacesize':['ank', 'ank_subspace'],
            'anklinearsolvetol':['ank', 'ank_rtol'],
            'ankasmoverlap':['ank', 'ank_asmoverlap'],
            'ankpcilufill':['ank', 'ank_ilufill'],
            'ankjacobianlag':['ank', 'ank_jacobianlag'],
            'ankinnerpreconits':['ank', 'ank_innerpreconits'],
            'ankcfl0':['ank', 'ank_cfl0'],

            # Load Balance Paramters
            'blocksplitting':['parallel', 'splitblocks'],
            'loadimbalance':['parallel', 'loadimbalance'],
            'loadbalanceiter':['parallel', 'loadbalanceiter'],
            'partitionlikenproc':['parallel', 'partitionlikenproc'],

            # Misc Paramters
            'printiterations':['iter', 'printiterations'],
            'printwarnings':['iter', 'printwarnings'],
            'printtiming':['adjoint', 'printtiming'],
            'setmonitor':['adjoint', 'setmonitor'],

            # Adjoint Params
            'adjointl2convergence':['adjoint', 'adjreltol'],
            'adjointl2convergencerel':['adjoint', 'adjreltolrel'],
            'adjointl2convergenceabs':['adjoint', 'adjabstol'],
            'adjointdivtol':['adjoint', 'adjdivtol'],
            'approxpc':['adjoint', 'approxpc'],
            'adpc':['adjoint', 'adpc'],
            'viscpc':['adjoint', 'viscpc'],
            'frozenturbulence':['adjoint', 'frozenturbulence'],
            'usediagtspc':['adjoint', 'usediagtspc'],
            'restartadjoint':['adjoint', 'restartadjoint'],
            'adjointsolver':{'gmres':'gmres',
                             'tfqmr':'tfqmr',
                             'richardson':'richardson',
                             'bcgs':'bcgs',
                             'ibcgs':'ibcgs',
                             'location':['adjoint', 'adjointsolvertype']},

            'adjointmaxiter':['adjoint', 'adjmaxiter'],
            'adjointsubspacesize':['adjoint', 'adjrestart'],
            'adjointmonitorstep':['adjoint', 'adjmonstep'],
            'dissipationlumpingparameter':['discr', 'sigma'],
            'preconditionerside':{'left':'left',
                                  'right':'right',
                                  'location':['adjoint', 'adjointpcside']},
            'matrixordering':{'natural':'natural',
                              'rcm':'rcm',
                              'nested dissection':'nd',
                              'one way dissection':'1wd',
                              'quotient minimum degree':'qmd',
                              'location':['adjoint', 'matrixordering']},

            'globalpreconditioner':{'additive schwartz':'asm',
                                    'multigrid':'mg',
                                    'location':['adjoint', 'precondtype']},
            'localpreconditioner':{'ilu':'ilu',
                                   'location':['adjoint', 'localpctype']},

            'ilufill':['adjoint', 'filllevel'],
            'applyadjointpcsubspacesize':['adjoint', 'applyadjointpcsubspacesize'],
            'asmoverlap':['adjoint', 'overlap'],
            'innerpreconits':['adjoint', 'innerpreconits'],
            'outerpreconits':['adjoint', 'outerpreconits'],
            'firstrun':['adjoint', 'firstrun'],
            'verifystate':['adjoint', 'verifystate'],
            'verifyspatial':['adjoint', 'verifyspatial'],
            'verifyextra':['adjoint', 'verifyextra'],
            'usematrixfreedrdw':['adjoint', 'usematrixfreedrdw'],

            # Parameters for functions
            'sepsensoroffset':['cost', 'sepsensoroffset'],
            'sepsensorsharpness':['cost', 'sepsensorsharpness'],
        }

        return optionMap, moduleMap

    def _getSpecialOptionLists(self):
        """
        Lists of special options
        """
        # These "ignore_options" are NOT actually', ignored, rather,
        # they DO NOT GET SET IN THE FORTRAN CODE. Rather, they are
        # used strictly in Python

        ignoreOptions = set(('numbersolutions',
                             'writesurfacesolution',
                             'writevolumesolution',
                             'writetecplotsurfacesolution',
                             'coupledsolution',
                             'autosolveretry',
                             'autoadjointretry',
                             'partitiononly',
                             'liftindex',
                             'meshsurfacefamily',
                             'designsurfacefamily',
                             'zippersurfacefamily',
                             'outputsurfacefamily',
                             'cutcallback',
                         ))

        # Deprecated options. These should not be used, but old
        # scripts can continue to run
        deprecatedOptions = {'finitedifferencepc':'Use the ADPC option.',
                             'writesolution':'Use writeSurfaceSolution and writeVolumeSolution options instead.',
                             }

        specialOptions = set(('surfacevariables',
                              'volumevariables',
                              'monitorvariables',
                              'outputdirectory',
                              'isovariables',
                              'isosurface',
                              'turbresscale',
                              'restartfile'
                          ))

        return ignoreOptions, deprecatedOptions, specialOptions

    def _getObjectivesAndDVs(self):
        iDV = {}
        iDV['alpha'] = self.adflow.adjointvars.ialpha
        iDV['beta'] = self.adflow.adjointvars.ibeta
        iDV['mach'] = self.adflow.adjointvars.imach
        iDV['machgrid'] = self.adflow.adjointvars.imachgrid
        iDV['p'] = self.adflow.adjointvars.ipressure
        iDV['rho'] = self.adflow.adjointvars.idensity
        iDV['t'] = self.adflow.adjointvars.itemperature
        iDV['rotx'] = self.adflow.adjointvars.irotx
        iDV['roty'] = self.adflow.adjointvars.iroty
        iDV['rotz'] = self.adflow.adjointvars.irotz
        iDV['rotcenx'] = self.adflow.adjointvars.irotcenx
        iDV['rotceny'] = self.adflow.adjointvars.irotceny
        iDV['rotcenz'] = self.adflow.adjointvars.irotcenz
        iDV['xref'] = self.adflow.adjointvars.ipointrefx
        iDV['yref'] = self.adflow.adjointvars.ipointrefy
        iDV['zref'] = self.adflow.adjointvars.ipointrefz

        # Convert to python indexing
        for key in iDV:
            iDV[key] = iDV[key] - 1 

        # This is ADflow's internal mapping for cost functions
        adflowCostFunctions = {
            'lift':self.adflow.costfunctions.costfunclift,
            'drag':self.adflow.costfunctions.costfuncdrag,
            'cl'  :self.adflow.costfunctions.costfuncliftcoef,
            'cd'  :self.adflow.costfunctions.costfuncdragcoef,
            'fx'  :self.adflow.costfunctions.costfuncforcex,
            'fy'  :self.adflow.costfunctions.costfuncforcey,
            'fz'  :self.adflow.costfunctions.costfuncforcez,
            'cfx' :self.adflow.costfunctions.costfuncforcexcoef,
            'cfy' :self.adflow.costfunctions.costfuncforceycoef,
            'cfz' :self.adflow.costfunctions.costfuncforcezcoef,
            'mx'  :self.adflow.costfunctions.costfuncmomx,
            'my'  :self.adflow.costfunctions.costfuncmomy,
            'mz'  :self.adflow.costfunctions.costfuncmomz,
            'cmx':self.adflow.costfunctions.costfuncmomxcoef,
            'cmy':self.adflow.costfunctions.costfuncmomycoef,
            'cmz':self.adflow.costfunctions.costfuncmomzcoef,
            'cm0':self.adflow.costfunctions.costfunccm0,
            'cmzalpha':self.adflow.costfunctions.costfunccmzalpha,
            'cmzalphadot':self.adflow.costfunctions.costfunccmzalphadot,
            'cl0':self.adflow.costfunctions.costfunccl0,
            'clalpha':self.adflow.costfunctions.costfuncclalpha,
            'clalphadot':self.adflow.costfunctions.costfuncclalphadot,
            'cfy0':self.adflow.costfunctions.costfunccfy0,
            'cfyalpha':self.adflow.costfunctions.costfunccfyalpha,
            'cfyalphdDot':self.adflow.costfunctions.costfunccfyalphadot,
            'cd0':self.adflow.costfunctions.costfunccd0,
            'cdalpha':self.adflow.costfunctions.costfunccdalpha,
            'cdalphadot':self.adflow.costfunctions.costfunccdalphadot,
            'cmzq':self.adflow.costfunctions.costfunccmzq,
            'cmzqdot':self.adflow.costfunctions.costfunccmzqdot,
            'clq':self.adflow.costfunctions.costfuncclq,
            'clqdot':self.adflow.costfunctions.costfuncclqdot,
            'cbend':self.adflow.costfunctions.costfuncbendingcoef,
            'sepsensor':self.adflow.costfunctions.costfuncsepsensor,
            'sepsensoravgx':self.adflow.costfunctions.costfuncsepsensoravgx,
            'sepsensoravgy':self.adflow.costfunctions.costfuncsepsensoravgy,
            'sepsensoravgz':self.adflow.costfunctions.costfuncsepsensoravgz,
            'cavitation':self.adflow.costfunctions.costfunccavitation,
            'mdot':self.adflow.costfunctions.costfuncmdot,
            'mavgptot':self.adflow.costfunctions.costfuncmavgptot,
            'mavgttot':self.adflow.costfunctions.costfuncmavgttot,
            'mavgps':self.adflow.costfunctions.costfuncmavgps
            }

        return iDV, adflowCostFunctions

    def _updateTurbResScale(self):
        # If turbresscale is None it has not been set by the user in script
        # thus set default values depending on turbulence model;
        # else do nothing since it already contains value specified by user

        if self.getOption("turbresscale") is None:
            turbModel = self.getOption("turbulencemodel")
            if turbModel == "sa":
                self.setOption("turbresscale", 10000.0)
            elif turbModel == "menter sst":
                self.setOption("turbresscale", [1e3, 1e-6])
            else:
                raise Error("Turbulence model %-35s does not have default values specified for turbresscale. Specify turbresscale manually or update the python interface"%(turbModel))

    def _setUnsteadyFileParameters(self):
        """
        This is a temporary function that sets the appropriate filenames
        when using the unsteady equations and the bdf time integration scheme.
        This function should drop out when unsteady solver is stepped through python
        TEMPORARY
        """
        # THIS DOES NOT CURRENTLY WORK DUE TO INTERNAL LOGIC. REFACTOR FORTRAN
        # Set parameters for outputing data
        #if self.getOption('writevolumesolution'):
        #    self.adflow.monitor.writevolume = True
        #    self.adflow.monitor.writegrid = True
        #else:
        #    self.adflow.monitor.writevolume = False
        #    self.adflow.monitor.writegrid = False

        #if self.getOption('writesurfacesolution'):
        #    self.adflow.monitor.writesurface = True
        #else:
        #    self.adflow.monitor.writesurface = False

        outputDir = self.getOption('outputDirectory')
        baseName = self.curAP.name

        # Join to get the actual filename root
        volFileName = os.path.join(outputDir, baseName + "_vol.cgns")
        surfFileName = os.path.join(outputDir, baseName + "_surf.cgns")
        sliceFileName = os.path.join(outputDir, baseName + "_slices")
        liftDistributionFileName = os.path.join(outputDir, baseName + "_lift")

        # Set fileName in adflow
        self.adflow.inputio.solfile[:] = ''
        self.adflow.inputio.solfile[0:len(volFileName)] = volFileName

        # Set the grid file to the same name so the grids will be written
        # to the volume files
        self.adflow.inputio.newgridfile[:] = ''
        self.adflow.inputio.newgridfile[0:len(volFileName)] = volFileName

        self.adflow.inputio.surfacesolfile[:] = ''
        self.adflow.inputio.surfacesolfile[0:len(surfFileName)] = surfFileName

        self.adflow.inputio.slicesolfile[:] = ''
        self.adflow.inputio.slicesolfile[0:len(sliceFileName)] = sliceFileName

        self.adflow.inputio.liftdistributionfile[:] = ''
        self.adflow.inputio.liftdistributionfile[0:len(liftDistributionFileName)] = liftDistributionFileName

    def _setForcedFileNames(self):
        # Set the filenames that will be used if the user forces a
        # write during a solution.

        self.adflow.inputio.forcedvolumefile[:] = ''
        self.adflow.inputio.forcedsurfacefile[:] = ''
        self.adflow.inputio.forcedliftfile[:] = ''
        self.adflow.inputio.forcedslicefile[:] = ''

        outputDir = self.getOption('outputDirectory')
        baseName = self.curAP.name
        base = os.path.join(outputDir, baseName)

        volFileName = base + '_forced_vol.cgns'
        surfFileName = base + '_forced_surf.cgns'
        liftFileName = base + '_forced_lift.dat'
        sliceFileName = base + '_forced_slices.dat'

        self.adflow.inputio.forcedvolumefile[0:len(volFileName)] = volFileName
        self.adflow.inputio.forcedsurfacefile[0:len(surfFileName)] = surfFileName
        self.adflow.inputio.forcedliftfile[0:len(liftFileName)] = liftFileName
        self.adflow.inputio.forcedslicefile[0:len(sliceFileName)] = sliceFileName

    def _createZipperMesh(self):
        """Internal routine for generating the zipper mesh. This operation is
        postposted as long as possible and now it cannot wait any longer."""
        if self.zipperCreated:
            return 

        zipFam = self.getOption('zipperSurfaceFamily')
        
        if zipFam is None:
            # The user didn't tell us anything. So we will use all
            # walls. Remind the user what those are. 
            zipperFamList = self.families[self.allWallsGroup]
            if self.myid == 0:
                ADFLOWWarning("'zipperSurfaceFamily' option was not given. Using all "
                              "wall boundary conditions for the zipper mesh.")
        else:

            if zipFam not in self.families:
                raise Error("Trying to create the zipper mesh, but '%s' is not a "
                            "family in the CGNS file or has not been added"
                            " as a combination of families"%groupName)
            zipperFamList = self.families[zipFam]

        self.adflow.zippermesh.createzippermesh(zipperFamList)
        self.zipperCreated = True

    def _processFortranStringArray(self, strArray):
        """Getting arrays of strings out of Fortran can be kinda nasty. This
        takes the array and returns a nice python list of strings"""
        shp = strArray.shape
        arr = strArray.reshape((shp[1],shp[0]), order='F')
        tmp = []
        for i in range(arr.shape[1]):
            tmp.append(''.join(arr[:, i]).strip().lower())

        return tmp

    def _createFortranStringArray(self, strList):
        """Setting arrays of strings in Fortran can be kinda nasty. This
        takesa list of strings and returns the array"""

        arr = numpy.zeros((len(strList),self.adflow.constants.maxcgnsnamelen), dtype="str")
        arr[:] = " "
        for i,s in enumerate(strList):
            for j in range(len(s)):
                arr[i,j] = s[j]
            
        return arr

    def createSlaveAeroProblem(self, master):
        """Create a slave aeroproblem"""

        # Make sure everything is created for the master
        self.setAeroProblem(master)

        slave = copy.deepcopy(master)
        slave.adflowData = master.adflowData
        slave.surfMesh = master.surfMesh
        slave.isSlave = True
        return slave

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
