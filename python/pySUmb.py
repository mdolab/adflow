#!/usr/bin/python
from __future__ import print_function
from __future__ import division
"""
pySUmb - A Python interface to SUmb.

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
from baseclasses import AeroSolver, AeroProblem
from . import MExt
from pprint import pprint as pp
    
class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pySUMB Error: '
        i = 15
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

class SUMBWarning(object):
    """
    Format a warning message
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| pySUMB Warning: '
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

# =============================================================================
# SUMB Class
# =============================================================================
class SUMB(AeroSolver):
    """
    Create the SUmb object.

    Parameters
    ----------
    comm : MPI intra comm
        The communicator on which to create SUmb. If not given, defaults
        to MPI.COMM_WORLD.
    options : dictionary
        The list of options to use with SUmb. This keyword arguement
        is NOT OPTIONAL. It must always be provided. It must contain, at least
        the 'gridFile' entry for the filename of the grid to load
    debug : bool
        Set this flag to true when debugging with a symbolic
        debugger. The MExt module deletes the copied .so file when not
        required which causes issues debugging.
        """
    def __init__(self, comm=None, options=None, debug=False):

        # Load the compiled module using MExt, allowing multiple
        # imports
        try:
            self.sumb
        except:
            curDir = os.path.dirname(os.path.realpath(__file__))
            self.sumb = MExt.MExt('libsumb', [curDir], debug=debug)._module

        # Information for base class:
        name = 'SUMB'
        category = 'Three Dimensional CFD'
        informs = {}

        # If 'options' is not None, go through and make sure all keys
        # are lower case:
        if options is not None:
            for key in options.keys():
                options[key.lower()] = options.pop(key)
        else:
            raise Error("The 'options' keyword argument must be passed "
                        "sumb. The options dictionary must contain (at least) "
                        "the gridFile entry for the grid")

        # Load all the option/objective/DV information:
        defOpts = self._getDefOptions()
        self.optionMap = self._getOptionMap()
        self.ignoreOptions, self.deprecatedOptions, self.specialOptions = \
                           self._getSpecialOptionLists()

        self.possibleAeroDVs, self.basicCostFunctions = (
            self._getObjectivesAndDVs())

        # Now add the group for each of the "basic" cost functions:
        self.sumbCostFunctions = {}
        for key in self.basicCostFunctions:
            self.sumbCostFunctions[key] = [None, key]

        # This is the real solver so dtype is 'd'
        self.dtype = 'd'

        # Next set the MPI Communicators and associated info
        if comm is None:
            comm = MPI.COMM_WORLD

        self.comm = comm
        self.sumb.communication.sumb_comm_world = self.comm.py2f()
        self.sumb.communication.sumb_comm_self = MPI.COMM_SELF.py2f()
        self.sumb.communication.sendrequests = numpy.zeros(self.comm.size)
        self.sumb.communication.recvrequests = numpy.zeros(self.comm.size)
        self.myid = self.sumb.communication.myid = self.comm.rank
        self.sumb.communication.nproc = self.comm.size

        # Initialize the inherited aerosolver
        AeroSolver.__init__(self, name, category, defOpts, informs,
                            options=options)

        # Setup the log file if necessary:
        self.logFile = None
        if self.getOption('logfile') != '' and self.comm.rank == 0:
            self.logFile = self.getOption('logfile')
            self.sumb.openlog(self.logFile)

        # Update turbresscale depending on the turbulence model specified
        self._updateTurbResScale()
        
        # Initialize petec in case the user has not already
        self.sumb.initializepetsc()

        # Set the stand-alone sumb flag to false...this changes how
        # terminate calls are handled.
        self.sumb.iteration.standalonemode = False

        # Set the frompython flag to true... this also changes how
        # terminate calls are handled
        self.sumb.killsignals.frompython = True

        # Dictionary of design varibales and their index
        self.aeroDVs = {}

        # Default counters
        self.updateTime = 0.0
        self.nSlice = 0
        self.nLiftDist = 0

        # Set default values
        self.sumb.setdefaultvalues()
        self.sumb.inputio.autoparameterupdate = False
        self._updateGeomInfo = True

        # By Default we don't have an external mesh object or a
        # geometric manipulation object, and the groupName is 'all'
        self.mesh = None
        self.DVGeo = None
        self.groupName = 'all'

        # Matrix Setup Flag
        self.adjointSetup = False

        # Write the intro message
        self.sumb.writeintromessage()

        # Make sure all the params are ok
        for option in self.options:
            if option != 'defaults':
                self.setOption(option.lower(), self.options[option][1])

        # Remind the user of all the sumb options:
        self.printCurrentOptions()

        # Do the remainder of the operations that would have been done
        # had we read in a param file
        self.sumb.iteration.deforming_grid = True

        # In order to properly initialize we need to have mach number
        # and a few other things set. Just create a dummy aeroproblem,
        # use it, and then it will be deleted.

        dummyAP = AeroProblem(name='dummy', mach=0.5, altitude=10000.0,
                              areaRef=1.0, chordRef=1.0, alpha=0.0, degreePol=0,
                              coefPol=[0, 0], degreeFourier=1,
                              omegaFourier=6.28, sinCoefFourier=[0, 0],
                              cosCoefFourier=[0, 0])

        self.curAP = dummyAP
        self._setAeroProblemData(firstCall=True)
        # Now set it back to None so the user is none the wiser
        self.curAP = None


        # Finally complete loading
        self.sumb.dummyreadparamfile()

        if self.getOption('partitionOnly'):
            self.sumb.partitionandreadgrid(True)
            return

        self.sumb.partitionandreadgrid(False)
        self.sumb.preprocessing()
        self.sumb.initflow()
        self.sumb.preprocessingpart2()
        self.sumb.initflowpart2()
        self.sumb.preprocessingadjoint()

    def setMesh(self, mesh, groupName='all'):
        """
        Set the mesh object to SUmb to do geometric deformations

        Parameters
        ----------
        mesh : multiBlockMesh object
            The pyWarp mesh object for doing the warping
        groupName : str
            The group which SUmb will consider to include all surfaces
            that will be manipulated by the geometry. 
        """

        self.mesh = mesh
        self.groupName = groupName

        # Setup External Warping
        ndof_1_instance = self.sumb.adjointvars.nnodeslocal[0]*3
        meshInd = self.sumb.getcgnsmeshindices(ndof_1_instance)
        self.mesh.setExternalMeshIndices(meshInd)

        # Setup Surface/Force info
        npatch = self.sumb.getnpatches()
        patchnames = []
        patchsizes = []
        for i in xrange(npatch):
            tmp = numpy.zeros(256, 'c')
            self.sumb.getpatchname(i+1, tmp)
            patchnames.append(
                ''.join([tmp[j] for j in range(256)]).lower().strip())
            patchsizes.append(self.sumb.getpatchsize(i+1))

        conn = self.getForceConnectivity()
        pts = self.getForcePoints()

        self.mesh.setExternalSurface(patchnames, patchsizes, conn, pts)
        # Get a inital copy of coordinates and save
        self.coords0 = self.getSurfaceCoordinates(self.groupName)

    def setDVGeo(self, DVGeo):
        """
        Set the DVGeometry object that will manipulate 'geometry' in
        this object. Note that SUmb doe not **strictly** need a
        DVGeometry object, but if optimization with geometric
        changes is desired, then it is required.

        Parameters
        ----------
        dvGeo : A DVGeometry object.
            Object responsible for manipulating the constraints that
            this object is responsible for.

        Examples
        --------
        >>> CFDsolver = SUMB(comm=comm, options=CFDoptions)
        >>> CFDsolver.setDVGeo(DVGeo)
        """

        self.DVGeo = DVGeo

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
        self.curAP.sumbData.disp = localD
                        


    def addLiftDistribution(self, nSegments, direction,
                            groupName=None, description=''):
        """
        Add a lift distribution to the surface output.

        Parameters
        ----------
        nSegments : int
            Number of slices to use for the distribution. Typically
            150-250 is sufficient
        direction : str
            One of 'x', 'y', or 'z'. The axis that is the 'lift' direction.

        groupName: str
            The family (as defined in pyWarp) to use for the lift
            distribution. Currently not coded.

        description : str
            An additional string that can be used to destingush
            between multiple lift distributions in the output.
        """
        if groupName is None:
            groupTag = ""
        else:
            groupTag = '%s: '% groupName

        direction=direction.lower()
        if direction not in ['x','y','z']:
            if direction not in ['x','y','z']:
                raise Error("direction must be one of 'x', 'y', or 'z'")

        # Determine the mask for the cells
        mask = self._getCellGroupMask(groupName)

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

        self.sumb.addliftdistribution(nSegments, dirVec, dirInd, distName, mask)

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
            The family (as defined in pyWarp) to use for the slices.
            """

        if groupName is None:
            groupTag = ""
        else:
            groupTag = '%s: '% groupName

        # Determine the mask for the cells
        mask = self._getCellGroupMask(groupName)

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
                self.sumb.addparaslice(sliceName, tmp[i], dirVec, mask)
            else:
                sliceName = 'Slice_%4.4d %s Absolute %s=%7.3f'% (
                    j, groupTag, direction, positions[i])
                self.sumb.addabsslice(sliceName, tmp[i], dirVec, mask)

        self.nSlice += N

    def addFunction(self, funcName, groupName, name=None):
        """Add a "new" function to SUmb by restricting the integration of an
        existing SUmb function by a section of the mesh defined by
        'groupName'. The function will be named 'funcName_groupName'
        provided that the 'name' keyword argument is not given. It is
        is, the function will be use that name. If necessary, this
        routine may be used to "change the name" of a function. For
        example,

        >>> addFunction('cd', None, 'super_cd')
    
        will add a function that is the same as 'cd', but call 'super_cd'
        supplied name 'name'.
        
        Parameters
        ----------
        funcName : str or list
            The name of the built-in sumb function
        groupName : str or list
            The family (as defined in warping module) to use for 
            the function.
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

        """Add a series of new functions to SUmb. This is a vector version of
        the addFunction() routine. See that routine for more documentation. """

        if len(funcNames) != len(groupNames) or len(funcNames) != len(names):
            raise Error("funcNames, groupNames, and names all have to be "
                        "lists of the same length")

        newFuncNames = []
        for i in range(len(funcNames)):
            funcName = funcNames[i]
            groupName = groupNames[i]
            name = names[i]

            # First make sure the supplied function is already known to sumb
            if funcName.lower() not in self.basicCostFunctions:
                raise Error('Supplied function name is not known to SUmb')

            # And make sure that the groupName is defined. 
            if self.mesh is None:
                raise Error("A mesh must be supplied to SUmb in order to use "
                            "the 'addFunction' functionality")

            # Check if the goupName has been added to the mesh.
            if groupName is not None:
                try:
                    famList = self.mesh.familyGroup[groupName]['families']
                except:
                    raise Error("The supplied groupName='%s' has not been "
                                "been added in the mesh"%groupName)
            if name is None:
                sumbFuncName = '%s_%s'%(funcName, groupName)
            else:
                sumbFuncName = name

            # Now register the function into sumbCostFunctions
            self.sumbCostFunctions[sumbFuncName] = [groupName, funcName]

            newFuncNames.append(name)

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
            a  = numpy.sqrt(self.sumb.flowvarrefstate.gammainf*\
                                self.sumb.flowvarrefstate.pinfdim/ \
                                self.sumb.flowvarrefstate.rhoinfdim)
            V = (self.sumb.inputphysics.machgrid+self.sumb.inputphysics.mach)*a

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
            cgnsBlocks = numpy.arange(1, self.sumb.cgnsgrid.cgnsndom+1)

        self.sumb.updaterotationrate(rotCenter, rotations, cgnsBlocks)
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
        loadInbalance, faceInbalance = self.sumb.checkpartitioning(nprocs)

        return loadInbalance, faceInbalance

    def getTriangulatedMeshSurface(self, groupName='all'):
        """
        This function returns a trianguled verision of the surface
        mesh on all processors. The intent is to use this for doing
        constraints in DVConstraints.

        Returns
        -------
        surf : list
           List of points and vectors describing the surface. This may
           be passed directly to DVConstraint setSurface() function.
        """

        # Use first spectral instance
        pts = self.comm.allgather(self.getForcePoints(0, groupName))
        conn = self.mesh.getSurfaceConnectivity(groupName)
        conn = self.comm.allgather(conn)

        # Triangle info...point and two vectors
        p0 = []
        v1 = []
        v2 = []

        for iProc in xrange(len(conn)):
            for i in xrange(len(conn[iProc])//4):
                i0 = conn[iProc][4*i+0]
                i1 = conn[iProc][4*i+1]
                i2 = conn[iProc][4*i+2]
                i3 = conn[iProc][4*i+3]

                p0.append(pts[iProc][i0])
                v1.append(pts[iProc][i1]-pts[iProc][i0])
                v2.append(pts[iProc][i3]-pts[iProc][i0])

                p0.append(pts[iProc][i2])
                v1.append(pts[iProc][i1]-pts[iProc][i2])
                v2.append(pts[iProc][i3]-pts[iProc][i2])

        return [p0, v1, v2]

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
        
        # Get option about adjoint memory
        releaseAdjointMemory = kwargs.pop('relaseAdjointMemory', True)

        # Set the aeroProblem
        self.setAeroProblem(aeroProblem, releaseAdjointMemory)

        # Set the full mask such that we get the full coefficients in
        # the printout.
        self.sumb.setfullmask()

        # If this problem has't been solved yet, reset flow to this
        # flight condition
        if self.curAP.sumbData.stateInfo is None:
            self.resetFlow(aeroProblem, releaseAdjointMemory)

        # Possibly release adjoint memory if not already done so.
        if releaseAdjointMemory:
            self.releaseAdjointMemory()

        # Clear out any saved adjoint RHS since they are now out of
        # data. Also increment the counter for this case.
        self.curAP.sumbData.adjointRHS = {}
        self.curAP.sumbData.callCounter += 1

        # --------------------------------------------------------------
        # Setup interation arrays ---- don't touch this unless you
        # REALLY REALLY know what you're doing!
        if self.sumb.monitor.niterold == 0 and \
            self.sumb.monitor.nitercur == 0 and \
            self.sumb.iteration.itertot == 0:
            if self.myid == 0:
                desiredSize = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles
                self.sumb.allocconvarrays(desiredSize)
        else:
            # More Time Steps / Iterations OR a restart
            # Reallocate convergence history array and time array
            # with new size,  storing old values from previous runs
            if self.getOption('storeHistory'):
                currentSize = len(self.sumb.monitor.convarray)
                desiredSize = currentSize + self.sumb.inputiteration.ncycles+1
                self.sumb.monitor.niterold  = self.sumb.monitor.nitercur+1
            else:
                self.sumb.monitor.nitercur  = 0
                self.sumb.monitor.niterold  = 1
                desiredSize = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles +1

            # Allocate Arrays
            if self.myid == 0:
                self.sumb.allocconvarrays(desiredSize)

            self.sumb.inputiteration.mgstartlevel = 1
            self.sumb.iteration.itertot = 0

        # --------------------------------------------------------------

        if self.getOption('equationMode') == 'unsteady':
            self.sumb.alloctimearrays(self.getOption('nTimeStepsFine'))

        # Mesh warp may have already failed:
        if not self.sumb.killsignals.fatalfail:

            if (self.getOption('equationMode').lower() == 'steady' or
                self.getOption('equationMode').lower() == 'time spectral'):
                self.updateGeometryInfo()

                # Check to see if the above update routines failed.
                self.sumb.killsignals.fatalfail = \
                    self.comm.allreduce(
                    bool(self.sumb.killsignals.fatalfail), op=MPI.LOR)

        if self.sumb.killsignals.fatalfail:
            print("Fatal failure during mesh warp! Bad mesh is "
                  "written in output directory as failed_mesh.cgns")
            fileName = os.path.join(self.getOption('outputDirectory'),
                                    'failed_mesh.cgns')
            self.writeMeshFile(fileName)
            self.curAP.fatalFail = True
            self.curAP.solveFailed = True
            return

        # We now now the mesh warping was ok so reset the flags:
        self.sumb.killsignals.routinefailed =  False
        self.sumb.killsignals.fatalfail = False

        t1 = time.time()

        # Call the Solver or the MD callback solver
        MDCallBack = kwargs.pop('MDCallBack', None)
        # When callback is not given, or time integration scheme is not md/ale
        # use normal solver instead
        if MDCallBack is None or not self.getOption('timeIntegrationScheme') == 'md' :
            self.sumb.solver()
        else:
            self.sumb.solverunsteady_ale(MDCallBack)

        # Save the states into the aeroProblem
        self.curAP.sumbData.stateInfo = self._getInfo()

        # Assign Fail Flags
        self.curAP.solveFailed = bool(self.sumb.killsignals.routinefailed)
        self.curAP.fatalFail = bool(self.sumb.killsignals.fatalfail)

        # Reset Flow if there's a fatal fail reset and return;
        # --> Do not write solution
        if self.curAP.fatalFail:
            self.resetFlow(aeroProblem)
            return

        t2 = time.time()
        solTime = t2 - t1

        if self.getOption('printTiming') and self.comm.rank == 0:
            print('Solution Time: %10.3f sec'% solTime)

        # Post-Processing -- Write Solutions is requested
        if kwargs.pop('writeSolution', True):
            self.writeSolution()

        if self.getOption('TSStability'):
            self.computeStabilityParameters()

    #======================================================
    #
    # Utilities for MD simulations
    #
    #======================================================

    def initSolver(self, aeroProblem, **kwargs):
        """
        This is essentially first half of __call__, which initializes the solver.
        Plus the routine that initializes zeroth unsteady time step.
        """
        
        # Get option about adjoint memory
        releaseAdjointMemory = kwargs.pop('relaseAdjointMemory', True)

        # Set the aeroProblem
        self.setAeroProblem(aeroProblem, releaseAdjointMemory)

        # If this problem has't been solved yet, reset flow to this
        # flight condition
        if self.curAP.sumbData.stateInfo is None:
            self.resetFlow(aeroProblem, releaseAdjointMemory)

        # Possibly release adjoint memory if not already done so.
        if releaseAdjointMemory:
            self.releaseAdjointMemory()

        # Clear out any saved adjoint RHS since they are now out of
        # data. Also increment the counter for this case.
        self.curAP.sumbData.adjointRHS = {}
        self.curAP.sumbData.callCounter += 1

        # --------------------------------------------------------------
        # Setup interation arrays ---- don't touch this unless you
        # REALLY REALLY know what you're doing!
        if self.sumb.monitor.niterold == 0 and \
            self.sumb.monitor.nitercur == 0 and \
            self.sumb.iteration.itertot == 0:
            if self.myid == 0:
                desiredSize = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles
                self.sumb.allocconvarrays(desiredSize)
        else:
            # More Time Steps / Iterations OR a restart
            # Reallocate convergence history array and time array
            # with new size,  storing old values from previous runs
            if self.getOption('storeHistory'):
                currentSize = len(self.sumb.monitor.convarray)
                desiredSize = currentSize + self.sumb.inputiteration.ncycles+1
                self.sumb.monitor.niterold  = self.sumb.monitor.nitercur+1
            else:
                self.sumb.monitor.nitercur  = 0
                self.sumb.monitor.niterold  = 1
                desiredSize = self.sumb.inputiteration.nsgstartup + \
                    self.sumb.inputiteration.ncycles +1

            # Allocate Arrays
            if self.myid == 0:
                self.sumb.allocconvarrays(desiredSize)

            self.sumb.inputiteration.mgstartlevel = 1
            self.sumb.iteration.itertot = 0

        # --------------------------------------------------------------

        if self.getOption('equationMode') == 'unsteady':
            self.sumb.alloctimearrays(self.getOption('nTimeStepsFine'))

        # Mesh warp may have already failed:
        if not self.sumb.killsignals.fatalfail:

            if (self.getOption('equationMode').lower() == 'steady' or
                self.getOption('equationMode').lower() == 'time spectral'):
                self.updateGeometryInfo()

                # Check to see if the above update routines failed.
                self.sumb.killsignals.fatalfail = \
                    self.comm.allreduce(
                    bool(self.sumb.killsignals.fatalfail), op=MPI.LOR)

        if self.sumb.killsignals.fatalfail:
            print("Fatal failure during mesh warp! Bad mesh is "
                  "written in output directory as failed_mesh.cgns")
            fileName = os.path.join(self.getOption('outputDirectory'),
                                    'failed_mesh.cgns')
            self.writeMeshFile(fileName)
            self.curAP.fatalFail = True
            self.curAP.solveFailed = True
            return

        # We now know the mesh warping was ok so reset the flags:
        self.sumb.killsignals.routinefailed =  False
        self.sumb.killsignals.fatalfail = False
        
        # Initialize for solverUnsteady_ALE
        self.sumb.solverunsteadywrapbegin()

    def stepCounter(self, adv=1):
        self.sumb.monitor.timestepunsteady = self.sumb.monitor.timestepunsteady + adv
        self.sumb.monitor.timeunsteady     = self.sumb.monitor.timeunsteady     + \
            adv*self.sumb.inputunsteady.deltat
        return self.sumb.monitor.timeunsteady, self.sumb.monitor.timestepunsteady
        
    def timeStepping(self):
        """
        This advances the solver by one physical time step
        """
        
        t1 = time.time()

        # Call the unsteady solver
        self.sumb.solverunsteadywrapinloop()

        # # Save the states into the aeroProblem
        # self.curAP.sumbData.stateInfo = self._getInfo()

        # # Assign Fail Flags
        # self.curAP.solveFailed = bool(self.sumb.killsignals.routinefailed)
        # self.curAP.fatalFail = bool(self.sumb.killsignals.fatalfail)

        # # Reset Flow if there's a fatal fail reset and return;
        # # --> Do not write solution
        # if self.curAP.fatalFail:
        #     self.resetFlow(aeroProblem)
        #     return

        t2 = time.time()
        solTime = t2 - t1

        ts     = self.sumb.monitor.timestepunsteady
        if self.getOption('writeVolumeSolution') and \
          numpy.mod(ts, self.getOption('nsavevolume')) == 0 :
            self.writeVolumeSolutionFile('V.cgns')
        if self.getOption('writeSurfaceSolution') and \
          numpy.mod(ts, self.getOption('nsavesurface')) == 0 :
            self.writeSurfaceSolutionFile('S.cgns')

        if self.getOption('printTiming') and self.comm.rank == 0:
            print('Solution Time: %10.3f sec'% solTime)

    def setDisplacement(self, coordinates, groupName='interface', ifShift=True):
        if ifShift and self.sumb.iteration.deforming_grid:
            self.sumb.shiftcoorandvolumes()
            self.sumb.shiftlevelale()
        self.setSurfaceCoordinates(coordinates, groupName)
        if self._updateGeomInfo and self.mesh is not None:
            timeA = time.time()
            self.mesh.warpMesh()
            newGrid = self.mesh.getSolverGrid()
            self.sumb.killsignals.routinefailed = False
            self.sumb.killsignals.fatalFail = False
            self.updateTime = time.time()-timeA
            if newGrid is not None:
                self.sumb.setgrid(newGrid)
            self._updateGeomInfo = False
            self.sumb.killsignals.routinefailed = \
                self.comm.allreduce(
                bool(self.sumb.killsignals.routinefailed), op=MPI.LOR)
            self.sumb.killsignals.fatalfail = self.sumb.killsignals.routinefailed

    def getCurForce(self):
    # Adapted from writeForceFile
    # Returns current tractions on the interface
        TS = 0
        groupName = 'interface'

        # Gather the data
        pts = self.comm.gather(self.getForcePoints(TS, groupName), root=0)

        # Forces are still evaluated on the displaced surface so do NOT pass in pts.
        forces = self.comm.gather(self.getForces(groupName, TS=TS), root=0)
        conn   = self.comm.gather(self.mesh.getSurfaceConnectivity(groupName), root=0)

        # Get data only on root proc
        f = [] # Tractions
        x = [] # Coordinates
        c = [] # Connectivity
        if self.myid == 0:
            # First sum up the total number of nodes and elements:
            nPt = 0
            nCell = 0
            for iProc in xrange(len(pts)):
                nPt += len(pts[iProc])
                nCell += len(conn[iProc])//4

            # Get the coordinates and surface tractions
            for iProc in xrange(len(pts)):
                for i in xrange(len(pts[iProc])):
                    f.append([ numpy.real(forces[iProc][i, 0]), \
                               numpy.real(forces[iProc][i, 1]), \
                               numpy.real(forces[iProc][i, 2]) ])
                    x.append([ numpy.real(pts[iProc][i, 0]), \
                               numpy.real(pts[iProc][i, 1]), \
                               numpy.real(pts[iProc][i, 2]) ])
            
            # Get the connectivity
            nodeOffset = 0
            for iProc in xrange(len(conn)):
                for i in xrange(len(conn[iProc])//4):
                    c.append([ int(conn[iProc][4*i+0]+nodeOffset), \
                               int(conn[iProc][4*i+1]+nodeOffset), \
                               int(conn[iProc][4*i+2]+nodeOffset), \
                               int(conn[iProc][4*i+3]+nodeOffset) ])
                nodeOffset += len(pts[iProc])
        return [f,x,c]
    
    #======================================================
    #
    # End of utilities for MD simulations
    #
    #======================================================
    
    def evalFunctions(self, aeroProblem, funcs, evalFuncs=None, sps=1,
                      ignoreMissing=False):
        """
        Evaluate the desired functions given in iterable object,
        'evalFuncs' and add them to the dictionary 'funcs'. The keys
        in the funcs dictioary will be have an _<ap.name> appended to
        them. Additionally, information regarding whether or not the
        last analysis with the aeroProblem was sucessful is
        included. This information is included as "funcs['fail']". If
        the 'fail' entry already exits in the dictionary the following
        operation is performed:

        funcs['fail'] = funcs['fail'] or <did this problem fail>

        In other words, if any one problem fails, the funcs['fail']
        entry will be False. This information can then be used
        directly in the pyOptSparse. 

        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to to get the solution for

        funcs : dict
            Dictionary into which the functions are saved.

        evalFuncs : iterable object containing strings
          If not None, use these functions to evaluate.

        sps : int, optional
            Spectral instance to use for funcions

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

        # We need to determine how many different masks we have, since
        # we will have to call getSolution for each *unique* function
        # mask. We can also do the error checking 
        groupMap = {}
        for f in evalFuncs:
            if f.lower() in self.sumbCostFunctions:            
                group = self.sumbCostFunctions[f][0]
                basicFunc = self.sumbCostFunctions[f][1]
                if group not in groupMap:
                    groupMap[group] = [[basicFunc, f]]
                else:
                    groupMap[group].append([basicFunc, f])
            else:
                if not ignoreMissing:
                    raise Error('Supplied function %s is not known to SUmb.'%f)

        # Now we loop over the unique groups calling the required
        # getSolution, there may be just one. No need for error
        # checking here.
        for group in groupMap:
            res = self.getSolution(sps, groupName=group)
            # g contains the "basic function" (index 0) and the actual
            # function name (index 1)
            for g in groupMap[group]:
                key = self.curAP.name + '_%s'% g[1]
                self.curAP.funcNames[g[1]] = key
                funcs[key] = res[g[0]]

    def checkSolutionFailure(self, aeroProblem, funcs):
        """
        Take in a an aeroProblem and check for failure. Then append the fail
        flag in funcs.
    
        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to to get the solution for

        funcs : dict
            Dictionary into which the functions are saved.
        """
        self.setAeroProblem(aeroProblem)
        # We also add the fail flag into the funcs dictionary. If fail
        # is already there, we just logically 'or' what was
        # there. Otherwise we add a new entry. 
        failFlag = self.curAP.solveFailed or self.curAP.fatalFail
        if 'fail' in funcs:
            funcs['fail'] = funcs['fail'] or failFlag
        else:
            funcs['fail'] = failFlag

    def evalFunctionsSens(self, aeroProblem, funcsSens, evalFuncs=None, sps=1):
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

        sps : int, optional
            Spectral instance to use for funcions

        Examples
        --------
        >>> funcSens = {}
        >>> CFDsolver.evalFunctionsSens(ap1, funcSens, ['cl', 'cd'])
        """

        # This is the one and only gateway to the getting derivatives
        # out of sumb. If you want a derivative, it should come from
        # here.

        self.setAeroProblem(aeroProblem)

        if evalFuncs is None:
            evalFuncs = self.curAP.evalFuncs
        else:
            evalFuncs = set(evalFuncs)

        # Do the functions one at a time:
        for f in evalFuncs:
            if f.lower() not in self.sumbCostFunctions:
                raise Error('Supplied %s function is not known to SUmb.'%f)

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
                

    def solveCL(self, aeroProblem, CLStar, alpha0=0,
                delta=0.5, tol=1e-3, autoReset=True):
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
            Initial guess for secant seach (deg)
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

        Returns
        -------
        None, but the correct alpha is stored in the aeroProblem
        """
        self.setAeroProblem(aeroProblem)

        anm2 = alpha0


        # Solve for the n-2 value:
        aeroProblem.alpha = anm2

        # Print alpha
        if self.comm.rank == 0:
            print ('Current alpha is: ', aeroProblem.alpha)

        self.__call__(aeroProblem, writeSolution=False)
        self.curAP.sumbData.callCounter -= 1
        sol = self.getSolution()
        fnm2 = sol['cl'] - CLStar
        if fnm2 < 0:
            anm1 = alpha0 + abs(delta)
        else:
            anm1 = alpha0 - abs(delta)

        minIterSave = self.getOption('minIterationNum')
        self.setOption('minIterationNum', 25)
        for iIter in range(20):
            # We need to reset the flow since changing the alpha leads
            # to problems with the NK solver
            if autoReset:
                self.resetFlow(aeroProblem)

            # Set current alpha
            aeroProblem.alpha = anm1

            # Print alpha
            if self.comm.rank == 0:
                print ('Current alpha is: ', aeroProblem.alpha)

            # Solve for n-1 value (anm1)
            self.__call__(aeroProblem, writeSolution=False)
            self.curAP.sumbData.callCounter -= 1
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

        # Restore the min iter option
        self.setOption('minIterationNum', minIterSave)

    def solveTrimCL(self, aeroProblem, trimFunc, trimDV, dvIndex,
                    CLStar, trimStar=0.0, alpha0=None, trim0=None, da=1e-3,
                    deta=1e-2, tol=1e-4, nIter=10):
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
        trim0 : float or NOne
            Starting trim value. If None, use what is in the DVGeo object
        da : float
            Initial alpha step for jacobian
        deta : float
            Initial  stet in the 'eta' or trim dv function
        tol : float
            Tolerance for trimCL solve solution
        nIter : int
            Maximum number of iterations. 
            """

        self.setAeroProblem(aeroProblem)

        # Set defaults if they are not given.                                                                                               
        if alpha0 is None:
            alpha0 = self.curAP.alpha
        if trim0 is None:
            trim0 = 0.0

        def Func(Xn, CLs):
            self.curAP.alpha = Xn[0]
            xGeo = self.DVGeo.getValues()
            xGeo[trimDV][dvIndex] = Xn[1]
            self.DVGeo.setDesignVars(xGeo)
            self.__call__(self.curAP)
            funcs = {}
            self.evalFunctions(self.curAP, funcs, evalFuncs=['cl',trimFunc])
            sol = self.getSolution()
            F = numpy.array([funcs['%s_cl'%self.curAP.name]-CLs, funcs['%s_%s'%(self.curAP.name, trimFunc)]])
            return F, sol

        # Generate initial point and jacobian

        Xn = numpy.array([alpha0, trim0])
        Fn, sol = Func(Xn, CLStar)

        # Next we generate the jacobian                                                                                                     
        J = numpy.zeros((2,2))

        # Perturb alpha                                                                                                                     
        Xn[0] += da
        Fpda, sol = Func(Xn, CLStar)
        J[:, 0] = (Fpda - Fn)/da
        Xn[0] -= da
        
        # Perturb eta:                                                                                                                      
        Xn[1] += deta
        Fpdeta, sol = Func(Xn, CLStar)
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
            Fnp1, sol = Func(Xnp1, CLStar)
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

    def solveSep(self, aeroProblem, sepStar, alpha0=0,
                delta=0.5, tol=1e-3, expansionRatio=1.2):
        """This is a safe-guarded secant search method to determine the alpha
        that yields a specified value of the sep sensor. Since this
        function is highly nonlinear we use a bisection search
        instead of a secant search. 

        Parameters
        ----------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to solve
        sepStar : float
            The desired target separation sensor value
        alpha0 : angle (deg)
            Initial guess 
        delta : angle (deg)
            Initial step. Only the magnitude is significant
        tol : float
            Desired tolerance for sepSensor

        Returns
        -------
        None, but the correct alpha is stored in the aeroProblem
        """ 
        self.setAeroProblem(aeroProblem)

        # There are two main parts of the algorithm: First we need to
        # find the alpha range bracketing the desired function
        # value. Then we use a safe-guarded secant search to zero into
        # that vlaue.

        minIterSave = self.getOption('minIterationNum')
        self.setOption('minIterationNum', 25)

        # Solve first problem
        aeroProblem.alpha = alpha0
        self.__call__(aeroProblem, writeSolution=False)
        sol = self.getSolution()
        f = sol['sepsensor'] - sepStar
        da = delta
        # Now try to find the interval
        for i in range(20):
            if f > 0.0:
                da = -numpy.abs(da)
            else:
                da = numpy.abs(da)
            # Increment the alpha and solve again
            aeroProblem.alpha += da
            self.__call__(aeroProblem, writeSolution=False)
            sol = self.getSolution()
            fnew = sol['sepsensor'] - sepStar

            if numpy.sign(f) != numpy.sign(fnew):
                # We crossed the zero:
                break
            else:
                f = fnew
                # Expand the delta:
                da *= expansionRatio

        # Now we know anm2 and anm1 bracket the zero. We now start the
        # secant search but make sure that we stay inside of these known bounds
        anm1 = aeroProblem.alpha
        anm2 = aeroProblem.alpha - da
        fnm1 = fnew
        fnm2 = f
        lowAlpha = min(anm1, anm2)
        highAlpha = max(anm1, anm2)
        for iIter in range(20):
            if iIter != 0:
                aeroProblem.alpha = anm1
                self.__call__(aeroProblem, writeSolution=False)
                sol = self.getSolution()
                fnm1 = sol['sepsensor'] - sepStar

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
        # Restore the min iter option
        self.setOption('minIterationNum', minIterSave)
            
    def writeSolution(self, outputDir=None, baseName=None, number=None):
        """This is a generic shell function that potentially writes
        the various output files. The intent is that the user or
        calling program can call this file and SUmb write all the
        files that the user has defined. It is recommneded that this
        function is used along with the associated logical flags in
        the options to determine the desired writing procedure

        Optional arguments

        outputDir: Use the supplied output directory

        baseName: Use this supplied string for the base filename. Typically
                  only used from an external solver.
        number: Use the user spplied number to index solutino. Again, only
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
                baseName = baseName + '_%3.3d'% self.curAP.sumbData.callCounter

        # Join to get the actual filename root
        base = os.path.join(outputDir, baseName)

        # Now call each of the 4 routines with the appropriate file name:
        if self.getOption('writevolumesolution'):
            self.writeVolumeSolutionFile(base + '_vol.cgns')
        if self.getOption('writesurfacesolution'):
            self.writeSurfaceSolutionFile(base + '_surf.cgns')

        self.writeLiftDistributionFile(base + '_lift.dat')
        self.writeSlicesFile(base + '_slices.dat')

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
        self.sumb.monitor.writegrid = True
        self.sumb.monitor.writevolume = False
        self.sumb.monitor.writesurface = False

        # Set filename in sumb
        self.sumb.inputio.solfile[:] = ''
        self.sumb.inputio.solfile[0:len(fileName)] = fileName

        self.sumb.inputio.newgridfile[:] = ''
        self.sumb.inputio.newgridfile[0:len(fileName)] = fileName

        # Actual fortran write call
        self.sumb.writesol()

    def writeVolumeSolutionFile(self, fileName, writeGrid=True):
        """Write the current state of the volume flow solution to a CGNS file.

        Parameters
        ----------
        fileName : str
            Name of the file. Should have .cgns extension.
        writeGrid : bool
            Flag specifying whether the grid should be included or if
            links should be used. Always writing the grid is
            recommended even in cases when it is not strictly necessary.
        """
        # Ensure extension is .cgns even if the user didn't specify
        fileName, ext = os.path.splitext(fileName)
        fileName += '.cgns'

        # Set Flags for writing
        self.sumb.monitor.writegrid = writeGrid
        self.sumb.monitor.writevolume = True
        self.sumb.monitor.writesurface = False

        # Set fileName in sumb
        self.sumb.inputio.solfile[:] = ''
        self.sumb.inputio.solfile[0:len(fileName)] = fileName

        self.sumb.inputio.newgridfile[:] = ''
        self.sumb.inputio.newgridfile[0:len(fileName)] = fileName

        # Actual fortran write call
        self.sumb.writesol()

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
        self.sumb.monitor.writegrid=False
        self.sumb.monitor.writevolume=False
        self.sumb.monitor.writesurface=True

        # Set fileName in sumb
        self.sumb.inputio.surfacesolfile[:] = ''
        self.sumb.inputio.surfacesolfile[0:len(fileName)] = fileName

        # Actual fortran write call
        self.sumb.writesol()

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
        self.sumb.writeliftdistributionfile(fileName)

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
        self.sumb.writeslicesfile(fileName)

    def writeForceFile(self, fileName, TS=0, groupName='all',
                       cfdForcePts=None):
        """This function collects all the forces and locations and
        writes them to a file with each line having: X Y Z Fx Fy Fz.
        This can then be used to set a set of structural loads in TACS
        for structural only optimization

        Like the getForces() routine, an external set of forces may be
        passed in on which to evaluate the forces. This is only
        typically used in an aerostructural case.

        """

        if self.mesh is None:
            mpiPrint('Error: A pyWarp mesh be specified to use writeForceFile',
                     comm=self.comm)
            return

        # Now we need to gather the data:
        if cfdForcePts is None:
            pts = self.comm.gather(self.getForcePoints(TS, groupName), root=0)
        else:
            pts = self.comm.gather(cfdForcePts)

        # Forces are still evaluated on the displaced surface so do NOT pass in pts.
        forces = self.comm.gather(self.getForces(groupName, TS=TS), root=0)
        conn   = self.comm.gather(self.mesh.getSurfaceConnectivity(groupName),
                                  root=0)

        # Write out Data only on root proc:
        if self.myid == 0:
            # First sum up the total number of nodes and elements:
            nPt = 0
            nCell = 0
            for iProc in xrange(len(pts)):
                nPt += len(pts[iProc])
                nCell += len(conn[iProc])//4

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
                for i in xrange(len(conn[iProc])//4):
                    f.write('%d %d %d %d\n'%(
                            conn[iProc][4*i+0]+nodeOffset,
                            conn[iProc][4*i+1]+nodeOffset,
                            conn[iProc][4*i+2]+nodeOffset,
                            conn[iProc][4*i+3]+nodeOffset))

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
        if obj in self.curAP.sumbData.adjoints:
            self.curAP.sumbData.adjoints[obj][:] = 0.0

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
        self.sumb.referencestate()
        self.sumb.setflowinfinitystate()

        strLvl = self.getOption('MGStartLevel')
        nLevels = self.sumb.inputiteration.nmglevels
        if strLvl < 0 or strLvl > nLevels:
            strLvl = nLevels

        self.sumb.inputiteration.mgstartlevel = strLvl
        self.sumb.iteration.groundlevel = strLvl
        self.sumb.iteration.currentlevel = strLvl
        self.sumb.monitor.niterold = 0
        self.sumb.monitor.nitercur = 0
        self.sumb.iteration.itertot = 0
        self.sumb.setuniformflow()
        self.sumb.nksolvervars.nksolvecount = 0
        self.sumb.killsignals.routinefailed =  False
        self.sumb.killsignals.fatalfail = False
        self.sumb.nksolvervars.freestreamresset = False

    def _setMask(self, groupName):
        """Set the wall masks for the supplied groupName"""
        
        if groupName is None:
            self.sumb.setfullmask()
            return

        try:
            famList = self.mesh.familyGroup[groupName]['families']
        except:
            raise Error("The supplied family group name has not "
                        "been added in the mesh object OR a mesh"
                        "object has not been supplied to sumb")

        # Set each of the families:
        self.sumb.setnmaskfams(len(famList))
        for i in range(len(famList)):
            self.sumb.setmaskfam(i+1, famList[i].upper())
        self.sumb.setmask()

    def getSolution(self, sps=1, groupName=None):
        """ Retrieve the solution variables from the solver. Note this
        is a collective function and must be called on all processors
        """
        # Set the mask for the group
        self._setMask(groupName)        
        
        # We should return the list of results that is the same as the
        # possibleObjectives list
        self.sumb.getsolution(sps)

        funcVals = self.sumb.costfunctions.functionvalue
        SUmbsolution = {
            'lift':funcVals[self.sumb.costfunctions.costfunclift-1],
            'drag':funcVals[self.sumb.costfunctions.costfuncdrag-1],
            'cl'  :funcVals[self.sumb.costfunctions.costfuncliftcoef-1],
            'cd'  :funcVals[self.sumb.costfunctions.costfuncdragcoef-1],
            'fx'  :funcVals[self.sumb.costfunctions.costfuncforcex-1],
            'fy'  :funcVals[self.sumb.costfunctions.costfuncforcey-1],
            'fz'  :funcVals[self.sumb.costfunctions.costfuncforcez-1],
            'cfx' :funcVals[self.sumb.costfunctions.costfuncforcexcoef-1],
            'cfy' :funcVals[self.sumb.costfunctions.costfuncforceycoef-1],
            'cfz' :funcVals[self.sumb.costfunctions.costfuncforcezcoef-1],
            'mx'  :funcVals[self.sumb.costfunctions.costfuncmomx-1],
            'my'  :funcVals[self.sumb.costfunctions.costfuncmomy-1],
            'mz'  :funcVals[self.sumb.costfunctions.costfuncmomz-1],
            'cmx' :funcVals[self.sumb.costfunctions.costfuncmomxcoef-1],
            'cmy' :funcVals[self.sumb.costfunctions.costfuncmomycoef-1],
            'cmz' :funcVals[self.sumb.costfunctions.costfuncmomzcoef-1],
            'cmzalphadot':funcVals[self.sumb.costfunctions.costfunccmzalphadot-1],
            'cmzalpha'   :funcVals[self.sumb.costfunctions.costfunccmzalpha-1],
            'cm0'        :funcVals[self.sumb.costfunctions.costfunccm0-1],
            'clalphadot' :funcVals[self.sumb.costfunctions.costfuncclalphadot-1],
            'clalpha'    :funcVals[self.sumb.costfunctions.costfuncclalpha-1],
            'cl0'        :funcVals[self.sumb.costfunctions.costfunccl0-1],
            'cfyalphadot':funcVals[self.sumb.costfunctions.costfunccfyalphadot-1],
            'cfyalpha'   :funcVals[self.sumb.costfunctions.costfunccfyalpha-1],
            'cfy0'       :funcVals[self.sumb.costfunctions.costfunccfy0-1],
            'cdalphadot' :funcVals[self.sumb.costfunctions.costfunccdalphadot-1],
            'cdalpha'    :funcVals[self.sumb.costfunctions.costfunccdalpha-1],
            'cd0'        :funcVals[self.sumb.costfunctions.costfunccd0-1],
            'cmzqdot'    :funcVals[self.sumb.costfunctions.costfunccmzqdot-1],
            'cmzq'       :funcVals[self.sumb.costfunctions.costfunccmzq-1],
            'clqdot'     :funcVals[self.sumb.costfunctions.costfuncclqdot-1],
            'clq'        :funcVals[self.sumb.costfunctions.costfuncclq-1],
            'cbend'      :funcVals[self.sumb.costfunctions.costfuncbendingcoef-1],
            'sepsensor'  :funcVals[self.sumb.costfunctions.costfuncsepsensor-1],
            'sepsensoravgx'  :funcVals[self.sumb.costfunctions.costfuncsepsensoravgx-1],
            'sepsensoravgy'  :funcVals[self.sumb.costfunctions.costfuncsepsensoravgy-1],
            'sepsensoravgz'  :funcVals[self.sumb.costfunctions.costfuncsepsensoravgz-1],
            'cavitation' :funcVals[self.sumb.costfunctions.costfunccavitation-1],
            }

        return SUmbsolution
        
    def printCurrentOptions(self):

        """
        Prints a nicely formatted dictionary of all the current SUmb
        options to the stdout on the root processor"""
        if self.comm.rank == 0:
            print('+---------------------------------------+')
            print('|          All SUmb Options:            |')
            print('+---------------------------------------+')
            # Need to assemble a temporary dictionary 
            tmpDict = {}
            for key in self.options:
                if key != 'defaults':
                    tmpDict[key] = self.getOption(key)
            pp(tmpDict)

    # =========================================================================
    #   The following routines are public functions however, they should
    #   not need to be used by a user using this class directly. They are
    #   primarly used internally OR by a solver that uses this class
    #   i.e. an Aerostructural solver
    # =========================================================================

    def _getCellGroupMask(self, groupName=None):
        """
        This function determines a mask (array of 1's and 0's) that
        determines if the cell in the globally reduced set of wall
        faces is a member of the groupName as defined by the warping.

        Parameters
        ----------
        groupName : str
            The name of the family group defined in the warping

        Returns
        -------
        mask : numpy array of integers
            The mask. Only returned on root proc. All other procs have
            numpy array of length 0.
            """

        localMask = self._getCellGroupMaskLocal(groupName)
        # Now we need to gather all of the local masks to the root
        # proc.
        tmp = self.comm.gather(localMask, root=0)
        globalMask = []
        if self.comm.rank == 0:
            for i in range(self.comm.size):
                globalMask.extend(tmp[i])
            globalMask = numpy.array(globalMask)
            
        return globalMask

    def _getCellGroupMaskLocal(self, groupName=None):
        """
        This function determines a mask (array of 1's and 0's) that
        determines if the cell in the local set of wall faces 
        is a member of the groupName as defined by the warping.

        Parameters
        ----------
        groupName : str
            The name of the family group defined in the warping

        Returns
        -------
        mask : numpy array of integers
            The mask. Returned on all procs. May be length 0.
            """

        # See if we can save ourselves some work if the mask is
        # already in the cache:
        if groupName is None:
            localMask = []
            for i in xrange(self.sumb.getnpatches()):
                patchsize = self.sumb.getpatchsize(i+1)
                nCells = patchsize[0]*patchsize[1]
                localMask.extend(numpy.ones(nCells, 'intc'))
        else:
            try:
                famList = self.mesh.familyGroup[groupName]['families']
            except:
                raise Error("The supplied family group name has not "
                            "been added in the mesh object.")

            # Now we just go through each patch and see if it in our familyList
            localMask = []
            for i in xrange(self.sumb.getnpatches()):
                tmp = numpy.zeros(256, 'c')
                self.sumb.getpatchname(i+1, tmp)
                patchname = ''.join([tmp[j] for j in range(256)]).lower().strip()
                patchsize = self.sumb.getpatchsize(i+1)
                nCells = (patchsize[0]-1)*(patchsize[1]-1)
                if patchname in famList:
                    localMask.extend(numpy.ones(nCells, 'intc'))
                else:
                    localMask.extend(numpy.zeros(nCells, 'intc'))
                
        return localMask

    def getSurfaceCoordinates(self, groupName='all'):
        """
        See MultiBlockMesh.py for more info
        """
        if self.mesh is not None:
            return self.mesh.getSurfaceCoordinates(groupName)
        else:
            raise Error('Must have a mesh object to use getSurfaceCoordinates.')  
    
    def getInitialSurfaceCoordinates(self, groupName='all'):
        """
        See MultiBlockMesh.py for more info
        """

        if self.mesh is not None:
            if self.DVGeo is not None:
                # if we have a geometry object, return the undeflected
                # shape generated directly from the design variables
                ptSetName = 'sumb_%s_coords'% self.curAP.name
                self.setSurfaceCoordinates(
                    self.DVGeo.update(ptSetName, config=self.curAP.name), self.groupName)
                self.updateGeometryInfo()
                return self.mesh.getSurfaceCoordinates(groupName)
            else:
                # otherwise, the initial mesh is the undeflected mesh, so 
                # return that
                return self.coords0
        else:
            raise Error('getInitialSurfaceCoordinateCoordinates requires a mesh object.')

    def setSurfaceCoordinates(self, coordinates, groupName='all'):
        """
        See MultiBlockMesh.py for more info
        """
        self._updateGeomInfo = True
        if self.mesh is not None:
            self.mesh.setSurfaceCoordinates(coordinates, groupName)

    def getSurfaceConnectivity(self, groupName='all'):
        """
        See MultiBlockMesh.py for more info
        """
        if self.mesh is not None:
            return self.mesh.getSurfaceConnectivity(groupName)

    def setAeroProblem(self, aeroProblem, releaseAdjointMemory=True):
        """Set the supplied aeroProblem to be used in SUmb"""
        ptSetName = 'sumb_%s_coords'% aeroProblem.name

        # If the aeroProblem is the same as the currently stored
        # reference; only check the DV changes
        if aeroProblem is self.curAP:
            if self.DVGeo is not None:
               
                # DVGeo appeared and we have not embedded points!
                if not ptSetName in self.DVGeo.points:
                    self.DVGeo.addPointSet(self.coords0, ptSetName)

                if not self.DVGeo.pointSetUpToDate(ptSetName):
                    coords = self.DVGeo.update(ptSetName, config=self.curAP.name)

                    if self.curAP.sumbData.disp is not None:
                        coords += self.curAP.sumbData.disp
                    self.setSurfaceCoordinates(coords, self.groupName)

            # Finally update other data
            self._setAeroProblemData()
            self.updateGeometryInfo()          

            return

        if self.comm.rank == 0:
            print('+'+'-'*70+'+')
            print('|  Switching to Aero Problem: %-41s|'% aeroProblem.name)
            print('+'+'-'*70+'+')

        # See if the aeroProblem has sumbData already, if not, create.
        try:
            aeroProblem.sumbData
        except AttributeError:
            aeroProblem.sumbData = sumbFlowCase()
            aeroProblem.ptSetName = ptSetName
            aeroProblem.surfMesh = None
            if self.mesh is not None:
                aeroProblem.surfMesh = self.getSurfaceCoordinates(self.groupName)
                
        if self.curAP is not None:
            # If we have already solved something and are now
            # switching, save what we need:
            self.curAP.stateInfo = self._getInfo()
            if self.mesh is not None:
                self.curAP.surfMesh = self.getSurfaceCoordinates(self.groupName)

        # If not done so already, embed the coordinates:
        if self.DVGeo is not None and ptSetName not in self.DVGeo.points:
            self.DVGeo.addPointSet(self.coords0, ptSetName)

        # We are now ready to associate self.curAP with the supplied AP
        self.curAP = aeroProblem
        self.curAP.adjointRHS = None

        # Now set the data from the incomming aeroProblem:
        stateInfo = aeroProblem.sumbData.stateInfo
        if stateInfo is not None:
            self._setInfo(stateInfo)

        # We have to update coordinates here as well:
        if self.DVGeo is not None:
            if not self.DVGeo.pointSetUpToDate(ptSetName):
                coords = self.DVGeo.update(ptSetName, config=self.curAP.name)
            else:
                coords = self.curAP.surfMesh
        else:
            coords = self.curAP.surfMesh
            
        if self.curAP.sumbData.disp is not None:
            coords += self.curAP.sumbData.disp

        self.setSurfaceCoordinates(coords, self.groupName)

        # Finally update other data
        self._setAeroProblemData()

        # Now we have to do a bunch of updates. This is fairly
        # expensive so switchign aeroProblems should not be done that
        #  frequently
        self.updateGeometryInfo()

        self.sumb.killsignals.routinefailed = False
        self.sumb.killsignals.fatalFail = False

        # Destroy the NK solver and the adjoint memory
        self.sumb.destroynksolver()
        if releaseAdjointMemory:
            self.releaseAdjointMemory()

      

    def _setAeroProblemData(self, firstCall=False):
        """
        After an aeroProblem has been associated with self.cuAP, set
        all the updated information in SUmb."""

        AP = self.curAP
        alpha = AP.alpha
        beta = AP.beta
        beta = AP.beta
        mach = AP.mach
        xRef = AP.xRef; yRef = AP.yRef; zRef = AP.zRef
        xRot = AP.xRot; yRot = AP.yRot; zRot = AP.zRot
        areaRef = AP.areaRef
        chordRef = AP.chordRef
        liftIndex = self.getOption('liftIndex')
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
                        " for SUmb.")
        if areaRef is None:
            raise Error("'areaRef' must be specified in aeroProblem"
                        " for SUmb.")
        if chordRef is None:
            raise Error("'chordRef' must be specified in aeroProblem"
                        " for SUmb.")

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

        if T is None or P is None or rho is None or V is None or mu is None:
            raise Error("Insufficient information is given in the "
                        "aeroProblem to determine physical state. "
                        "See AeroProblem documentation for how to "
                        "specify complete aerodynamic states.")

        # 1. Angle of attack:
        dToR = numpy.pi/180.0
        self.sumb.adjustinflowangle(alpha*dToR, beta*dToR, liftIndex)
        if self.getOption('printIterations') and self.comm.rank == 0:
            print('-> Alpha... %f '% numpy.real(alpha))

        # 2. Reference Points:
        self.sumb.inputphysics.pointref = [xRef, yRef, zRef]
        self.sumb.inputmotion.rotpoint = [xRot, yRot, zRot]

        # 3. Reference Areas
        self.sumb.inputphysics.surfaceref = areaRef
        self.sumb.inputphysics.lengthref = chordRef

        # 4. Set mach numbers
        if self.getOption('equationMode').lower()=='time spectral':
            self.sumb.inputphysics.mach = 0.0
            self.sumb.inputphysics.machcoef = mach
            self.sumb.inputphysics.machgrid = mach
            self.sumb.inputmotion.gridmotionspecified = True
        else:
            self.sumb.inputphysics.mach = mach
            self.sumb.inputphysics.machcoef = mach
            self.sumb.inputphysics.machgrid = 0.0

        # Set reference state information:
        if self.getOption('equationType') != 'euler':
            # RANS
            self.sumb.flowvarrefstate.pref = P
            self.sumb.flowvarrefstate.tref = T
            self.sumb.inputphysics.tempfreestream = T
            ReLength = 1.0
            self.sumb.inputphysics.reynolds = rho*V/mu
            self.sumb.inputphysics.reynoldslength = ReLength

            self.sumb.inputphysics.ssuthdim = SSuthDim
            self.sumb.inputphysics.musuthdim = muSuthDim
            self.sumb.inputphysics.tsuthdim = TSuthDim
            self.sumb.inputphysics.rgasdim = RGasDim
            
            # Update gamma only if it has changed from what currently is set
            if abs(self.sumb.inputphysics.gammaconstant - gammaConstant) > 1.0e-12:
                self.sumb.inputphysics.gammaconstant = gammaConstant
                self.sumb.updategamma() # NOTE! It is absolutely necessary to call this function, otherwise gamma is not properly updated.
            self.sumb.inputphysics.prandtl = Pr            
        else:
            # EULER
            self.sumb.flowvarrefstate.pref = P
            self.sumb.flowvarrefstate.tref = T
            self.sumb.inputphysics.tempfreestream = T
            self.sumb.inputphysics.reynolds = 1.0
            self.sumb.inputphysics.reynoldslength = 1.0

            self.sumb.inputphysics.rgasdim = RGasDim            
            # Update gamma only if it has changed from what currently is set
            if abs(self.sumb.inputphysics.gammaconstant - gammaConstant) > 1.0e-12:
                self.sumb.inputphysics.gammaconstant = gammaConstant
                self.sumb.updategamma() # NOTE! It is absolutely necessary to call this function, otherwise gamma is not properly updated.
            self.sumb.inputphysics.prandtl = Pr


        # 4. Periodic Parameters --- These are not checked/verified
        # and come directly from aeroProblem. Make sure you specify
        # them there properly!!
        if  self.getOption('alphaMode'):
            self.sumb.inputmotion.degreepolalpha = int(AP.degreePol)
            self.sumb.inputmotion.coefpolalpha = AP.coefPol
            self.sumb.inputmotion.omegafouralpha   = AP.omegaFourier
            self.sumb.inputmotion.degreefouralpha  = AP.degreeFourier
            self.sumb.inputmotion.coscoeffouralpha = AP.cosCoefFourier
            self.sumb.inputmotion.sincoeffouralpha = AP.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif  self.getOption('betaMode'):
            self.sumb.inputmotion.degreepolmach = int(AP.degreePol)
            self.sumb.inputmotion.coefpolmach = AP.coefPol
            self.sumb.inputmotion.omegafourbeta   = AP.omegaFourier
            self.sumb.inputmotion.degreefourbeta  = AP.degreeFourier
            self.sumb.inputmotion.coscoeffourbeta = AP.cosCoefFourier
            self.sumb.inputmotion.sincoeffourbeta = AP.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('machMode'):
            self.sumb.inputmotion.degreepolmach = int(AP.degreePol)
            self.sumb.inputmotion.coefpolmach = AP.coefPol
            self.sumb.inputmotion.omegafourmach   = AP.omegaFourier
            self.sumb.inputmotion.degreefourmach  = AP.degreeFourier
            self.sumb.inputmotion.coscoeffourmach = AP.cosCoefFourier
            self.sumb.inputmotion.sincoeffourmach = AP.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif  self.getOption('pMode'):
            ### add in lift axis dependence
            self.sumb.inputmotion.degreepolxrot = int(AP.degreePol)
            self.sumb.inputmotion.coefpolxrot = AP.coefPol
            self.sumb.inputmotion.omegafourxrot = AP.omegaFourier
            self.sumb.inputmotion.degreefourxrot  = AP.degreeFourier
            self.sumb.inputmotion.coscoeffourxrot = AP.cosCoefFourier
            self.sumb.inputmotion.sincoeffourxrot = AP.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('qMode'):
            self.sumb.inputmotion.degreepolzrot = int(AP.degreePol)
            self.sumb.inputmotion.coefpolzrot = AP.coefPol
            self.sumb.inputmotion.omegafourzrot = AP.omegaFourier
            self.sumb.inputmotion.degreefourzrot  = AP.degreeFourier
            self.sumb.inputmotion.coscoeffourzrot = AP.cosCoefFourier
            self.sumb.inputmotion.sincoeffourzrot = AP.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True
        elif self.getOption('rMode'):
            self.sumb.inputmotion.degreepolyrot = int(AP.degreePol)
            self.sumb.inputmotion.coefpolyrot = AP.coefPol
            self.sumb.inputmotion.omegafouryrot = AP.omegaFourier
            self.sumb.inputmotion.degreefouryrot  = AP.degreeFourier
            self.sumb.inputmotion.coscoeffouryrot = AP.cosCoefFourier
            self.sumb.inputmotion.sincoeffouryrot = AP.sinCoefFourier
            self.sumb.inputmotion.gridmotionspecified = True

        if not firstCall:
            self.sumb.referencestate()
            self.sumb.setflowinfinitystate()
            self.sumb.iteration.groundlevel = 1
            self.sumb.updateperiodicinfoalllevels()
            self.sumb.updategridvelocitiesalllevels()

    def getForces(self, groupName=None, TS=0):
        """ Return the forces on this processor.

        Parameters
        ----------
        groupName : str
            Group identifier to get only forces cooresponding to the
            desired group.  The group must have been added with
            addFamilyGroup in the mesh object.
        TS : int
            Spectral instance for which to get the forces
        pressure : bool
            Flag as to include the pressure forces
        viscous : bool
            Flag as to include the viscous forces

        Returns
        -------
        forces : array (N,3)
            Forces (or tractions depending on that forceAsTractions
            options) on this processor. Note that N may be 0, and an
            empty array of shape (0, 3) can be returned.
        """
        # Just to be sure, we'll set the mask
        self._setMask(groupName)

        [npts, ncell] = self.sumb.getforcesize()
        if npts > 0:
            forces = self.sumb.getforces(npts, TS+1).T
        else:
            forces = numpy.zeros((0,3),self.dtype)

        if groupName is not None:
            forces = self.mesh.sectionVectorByFamily(groupName, forces)

        return forces

    def getForcePoints(self, TS=0, groupName=None):
        [npts, ncell] = self.sumb.getforcesize()
        pts = numpy.zeros((npts, 3),self.dtype)
        if npts > 0:
            self.sumb.getforcepoints(pts.T, TS+1)

        if groupName is not None:
            pts = self.mesh.sectionVectorByFamily(groupName, pts)

        return pts

    def getForceConnectivity(self):
        [npts, ncells] = self.sumb.getforcesize()
        conn =  numpy.zeros((ncells, 4), dtype='intc')
        self.sumb.getforceconnectivity(numpy.ravel(conn))

        return conn

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

        return self.sumb.applypc(inVec, outVec)

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

        return self.sumb.applyadjointpc(inVec, outVec)

    def verifyAD(self):
        """
        Use Tapenade TGT debugger to verify AD
        """
        self.sumb.verifyad()

    def _addAeroDV(self, dv):
        """Add a single desgin variable that SUmb knows about.

        Parameters
        ----------
        dv : str
            dv name. Must be in the self.possibleAeroDVs list
            """

        if dv not in self.possibleAeroDVs:
            raise Error("%s was not one of the possible AeroDVs. "
                        "The complete list of DVs for SUmb is %s. "%(
                            dv, repr(set(self.possibleAeroDVs.keys()))))

        if dv not in self.aeroDVs:
            # A new DV add it:
            self.aeroDVs[dv] = len(self.aeroDVs)

    def _setupAdjoint(self, reform=False):
        """
        Setup the data structures required to solve the adjoint problem
        """

        # Destroy the NKsolver to free memory -- Call this even if the
        # solver is not used...a safeguard check is done in Fortran
        self.sumb.destroynksolver()
        self._setAeroDVs()

        if not self.adjointSetup or reform:
            # Create any PETSc variables if necessary
            self.sumb.createpetscvars()

            if self.getOption('useReverseModeAD'):
                raise Error('The old reverse mode AD routines have been '
                            'deprecated.')

            # Setup all required matrices in forward mode. (possibly none
            # of them)
            self.sumb.setupallresidualmatricesfwd()
            self.sumb.setuppetscksp()
            self.adjointSetup = True

    def releaseAdjointMemory(self):
        """
        release the PETSc Memory that have been allocated
        """
        if self.adjointSetup:
            self.sumb.destroypetscvars()
            self.adjointSetup = False

    def solveAdjoint(self, aeroProblem, objective, forcePoints=None,
                      structAdjoint=None, groupName=None):
        # May be switching aeroProblems here
        self.setAeroProblem(aeroProblem)

        # Possibly setup adjoint matrices/vector as required
        self._setupAdjoint()

        # # Check to see if the RHS Partials have been computed
        if objective not in self.curAP.sumbData.adjointRHS:
            RHS = self.computeJacobianVectorProductBwd(
                funcsBar={objective.lower():1.0}, wDeriv=True)
            self.curAP.sumbData.adjointRHS[objective] = RHS.copy()
        else:
            RHS = self.curAP.sumbData.adjointRHS[objective].copy()

        # Check to see if we need to agument the RHS with a structural
        # adjoint:
        if structAdjoint is not None and groupName is not None:
            phi = self.mesh.expandVectorByFamily(groupName, structAdjoint)
            agument = self.computeJacobianVectorProductBwd(
                fBar=phi, wDeriv=True)
            RHS -= agument

        # Check if objective is python 'allocated':
        if objective not in self.curAP.sumbData.adjoints:
            self.curAP.sumbData.adjoints[objective] = (
                numpy.zeros(self.getAdjointStateSize(), float))

        # Extract the psi:
        psi = self.curAP.sumbData.adjoints[objective]

        # Actually Solve the adjoint system...psi is updated with the
        # new solution.
        self.sumb.solveadjoint(RHS, psi, True)

        # Now set the flags and possibly reset adjoint
        if self.sumb.killsignals.adjointfailed:
            self.adjointFailed = True
            # Reset stored adjoint
            self.curAP.sumbData.adjoints[objective][:] = 0.0
        else:
            self.curAP.sumbData.adjoints[objective] = psi
            self.adjointFailed = False

    def _processAeroDerivatives(self, dIda):
        """This internal furncion is used to convert the raw array ouput from
        the matrix-free product bwd routine into the required
        dictionary format."""
        
	funcsSens = {}

        DVsRequired = list(self.curAP.DVNames.keys())
        for dv in DVsRequired:
            if dv.lower() in ['altitude', 'mach', 'P', 'T', 'reynolds']:
                # These design variables are *special*! The issue is
                # that SUmb takes as effective input P, T, mach and
                # reynolds and these are *not* the typical values we
                # use from the aeroproblem. So to ensure generic
                # treatment of the derivatives we have already added
                # P, T, mach and Reynolds to SUmb as "design
                # variables" --- this means we get the sensitivity of
                # the objective with respect to these 4
                # variables. Next we compute the derivatives of
                # P,T,mach,rho,V, and mu with respect to the actual
                # aeroproblem design variable 'dv'. We then chain rule
                # them together, along with the linearization fo
                # Re=rho*V/mu. 
                tmp = {}
	        self.curAP.evalFunctionsSens(tmp, ['P', 'T', 'mach', 'rho','V', 'mu'])

                # Extract the derivatives wrt the independent
                # parameters in SUmb       
                dIdP = dIda[self.aeroDVs['P']]
                dIdT = dIda[self.aeroDVs['T']]
                dIdMach = dIda[self.aeroDVs['mach']]
                # Chain rule the reynolds dependance back to what came from aeroproblem:
                dIdReynolds = dIda[self.aeroDVs['reynolds']]
                dIdV = self.curAP.rho/self.curAP.mu*dIdReynolds
                dIdrho = self.curAP.V/self.curAP.mu*dIdReynolds
                dIdmu = -self.curAP.rho*self.curAP.V/self.curAP.mu**2 * dIdReynolds

                # Chain-rule to get the final derivative:
                funcsSens[self.curAP.DVNames[dv]] = (
                    tmp[self.curAP['P']][self.curAP.DVNames[dv]]*dIdP +
                    tmp[self.curAP['T']][self.curAP.DVNames[dv]]*dIdT +
                    tmp[self.curAP['mach']][self.curAP.DVNames[dv]]*dIdMach + 
                    tmp[self.curAP['rho']][self.curAP.DVNames[dv]]*dIdrho + 
                    tmp[self.curAP['V']][self.curAP.DVNames[dv]]*dIdV + 
                    tmp[self.curAP['mu']][self.curAP.DVNames[dv]]*dIdmu)

            elif dv in self.possibleAeroDVs:
                if dv in self.possibleAeroDVs:
                    funcsSens[self.curAP.DVNames[dv]] = dIda[self.aeroDVs[dv]]
                    if dv == 'alpha':
                        funcsSens[self.curAP.DVNames[dv]] *= numpy.pi/180.0
     
	return funcsSens


    def _setAeroDVs(self):

        """ Do everything that is required to deal with aerodynamic
        design variables in SUmb"""

        DVsRequired = list(self.curAP.DVNames.keys())
        DVMap = {}
        for dv in DVsRequired:
            if dv.lower() in ['altitude', 'mach', 'P', 'T', 'reynolds']:
                # All these variables need to be compined 
                self._addAeroDV('P')
                self._addAeroDV('mach')
                self._addAeroDV('reynolds')
                self._addAeroDV('T')
            elif dv in self.possibleAeroDVs:
                self._addAeroDV(dv)
            else:
                raise Error("The design variable '%s' as specified in the"
                            " aeroProblem cannot be used with SUmb."% dv)

        # Reset all:
        for pdv in self.possibleAeroDVs:
            execStr = 'self.sumb.' + self.possibleAeroDVs[pdv] + '=-1'
            exec(execStr)

        # Set the required paramters for the aero-Only design vars:
        self.nDVAero = len(self.aeroDVs)
        self.sumb.adjointvars.ndesignextra = self.nDVAero
        self.sumb.adjointvars.dida = numpy.zeros(self.nDVAero)
        for adv in self.aeroDVs:
            execStr = 'self.sumb.' + self.possibleAeroDVs[adv] + \
                      '= %d'% self.aeroDVs[adv]
            # Leave this zero-based since we only need to use it in petsc
            exec(execStr)

    def getNAeroDV(self):
        """
        Return the number of aerodynamic variables
        """
        self._setAeroDVs()
        return self.nDVAero

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
        outVec = self.sumb.solveadjointforrhs(inVec, relTol)

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
        outVec = self.sumb.solvedirectforrhs(inVec, relTol)

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
        self.sumb.saveadjointmatrix(adjointMatrixName)
        self.sumb.saveadjointpc(pcMatrixName)
        self.sumb.saveadjointrhs(rhsName)

    def computeStabilityParameters(self):
        """
        run the stability derivative driver to compute the stability parameters
        from the time spectral solution
        """
        self.sumb.stabilityderivativedriver()

    def updateGeometryInfo(self):
        """Update the SUmb internal geometry info, if necessary."""

        if self._updateGeomInfo and self.mesh is not None:
            timeA = time.time()
            self.mesh.warpMesh()
            newGrid = self.mesh.getSolverGrid()
            self.sumb.killsignals.routinefailed = False
            self.sumb.killsignals.fatalFail = False
            self.updateTime = time.time()-timeA
            if newGrid is not None:
                self.sumb.setgrid(newGrid)
            self.sumb.updatecoordinatesalllevels()
            self.sumb.updatewalldistancealllevels()
            self.sumb.updateslidingalllevels()
            self.sumb.updatemetricsalllevels()
            self.sumb.updategridvelocitiesalllevels()
            self._updateGeomInfo = False
            self.sumb.killsignals.routinefailed = \
                self.comm.allreduce(
                bool(self.sumb.killsignals.routinefailed), op=MPI.LOR)
            self.sumb.killsignals.fatalfail = self.sumb.killsignals.routinefailed

    def getAdjointResNorms(self):
        '''
        Return the following adjoint residual norms:
        initRes Norm: Norm the adjoint RHS
        startRes Norm: Norm at the start of adjoint call (with possible non-zero restart)
        finalCFD Norm: Norm at the end of adjoint solve
        '''
        initRes  = self.sumb.adjointpetsc.adjresinit
        startRes = self.sumb.adjointpetsc.adjresstart
        finalRes = self.sumb.adjointpetsc.adjresfinal
        fail = self.sumb.killsignals.adjointfailed

        return initRes, startRes, finalRes, fail

    def getResNorms(self):
        """Return the initial, starting and final Res Norms. Typically
        used by an external solver."""
        return \
            numpy.real(self.sumb.nksolvervars.totalr0), \
            numpy.real(self.sumb.nksolvervars.totalrstart),\
            numpy.real(self.sumb.nksolvervars.totalrfinal)

    def setResNorms(self, initNorm=None, startNorm=None, finalNorm=None):
        """ Set one of these norms if not None. Typlically used by an
        external solver"""
        if initNorm is not None:
            self.sumb.nksolvervars.totalr0 = initNorm
        if startNorm is not None:
            self.sumb.nksolvervars.totalrstart = startNorm
        if finalNorm is not None:
            self.sumb.nksolvervars.finalNorm = finalNorm

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
        if the string is a valid function in SUmb. The function returns
        a boolean.

        Parameters
        ----------
        objective: string
            The objective to be tested.        
        """
        if objective.lower() in self.sumbCostFunctions.keys():
            return True
        else:
            return False

    # =========================================================================
    #   The following routines two routines are the workhorse of the
    #   forward and reverse mode routines. These can compute *ANY*
    #   product that is possible from the solver. 
    #   =========================================================================

    def computeJacobianVectorProductFwd(self, xDvDot=None, xSDot=None, xVDot=None, wDot=None,
                                    residualDeriv=False, funcDeriv=False, fDeriv=False):
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

        Returns
        -------
        dwdot, funcsdot, fDot : array, dict, array
            One or more of the these are return depending ont he *Deriv flags
        """

        if xDvDot is None and xSDot is None and xVDot is None and wDot is None:
            raise Error('computeJacobianVectorProductFwd: xDvDot, xSDot, xVDot and wDot cannot '
                        'all be None')
            
        self._setAeroDVs()

        # Default flags
        useState = False
        useSpatial = False

        # Process the Xs perturbation
        if xSDot is None:
            xsdot = numpy.zeros_like(self.coords0)
            useSpatial = True
        else:
            xsdot = xSDot

        # Process the Xv perturbation
        if xVDot is None:
            xvdot = numpy.zeros(self.getSpatialSize())
            useSpatial = True
        else:
            xvdot = xVDot

        # Process wDot perturbation
        if wDot is None:
            wdot = numpy.zeros(self.getStateSize())
        else:
            wdot = wDot
            useState = True

        # Process the extra variable perturbation....this comes from
        # xDvDot
        extradot = numpy.zeros(max(1, self.nDVAero))
        if xDvDot is not None:
            useSpatial = True
            for key in self.aeroDVs:
                execStr = 'mapping = self.sumb.%s'%self.possibleAeroDVs[key.lower()]
                exec(execStr)
                if key == 'alpha':
                    # convert angle of attack to radians
                    extradot[mapping] = xDvDot[key]*(numpy.pi/180.0)
                else:
                    extradot[mapping] = xDvDot[key]

        # For the geometric xDvDot perturbation we accumulate into the
        # already existing (and possibly nonzero) xsdot and xvdot
        if xDvDot is not None or xSDot is not None:
            if xDvDot is not None:
                xsdot += self.DVGeo.totalSensitivityProd(xDvDot, self.curAP.ptSetName)
            xvdot += self.mesh.warpDerivFwd(xsdot)
            useSpatial = True

        # Sizes for output arrays
        costSize = self.sumb.costfunctions.ncostfunction
        fSize, nCell = self.sumb.getforcesize()
        
        dwdot,tmp,fdot = self.sumb.computematrixfreeproductfwd(
            xvdot, extradot, wdot, useSpatial, useState, costSize, fSize)

        # Process the derivative of the functions
        funcsdot = {}
        for f in self.curAP.evalFuncs:
            mapping = self.sumbCostFunctions[f.lower()]
            funcsdot[f] = tmp[mapping - 1]
        
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
            Seed for the residuals (dwb in sumb)
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

        # ---------------------
        #  Check for funcsBar
        # ---------------------
        self.sumb.setfullmask()
        if funcsBar is None:
            funcsBar = numpy.zeros(self.sumb.costfunctions.ncostfunction)
        else:
            tmp = numpy.zeros(self.sumb.costfunctions.ncostfunction)

            # We have to make sure that the user has supplied
            # functions that have the same mask, otherwise, we can't
            # do it
            groupMasks = set()
            
            for f in funcsBar:
                if f.lower() in self.sumbCostFunctions:

                    mask = self.sumbCostFunctions[f][0]
                    basicFunc = self.sumbCostFunctions[f][1]
                    groupMasks.add(mask)

                    mapping = self.basicCostFunctions[basicFunc]
                    tmp[mapping-1] = funcsBar[f]

            if len(groupMasks) == 1:
                # We're ok..there was only one mask from the
                # functions...just pop it out of the set and then set
                # it.
                self._setMask(groupMasks.pop())
            elif len(groupMasks) > 1:
                raise Error("Attemping to compute a jacobian vector product "
                            "with multiple functions that have different "
                            "masks. This is not allowed.")
            
            # If there wasn't any actual funcsBar, then tmp is still
            # just zeros
            funcsBar = tmp

        # -------------------------
        #  Check for fBar (forces)
        # -------------------------
        if fBar is None:
            [nPts, nCell] = self.sumb.getforcesize()
            fBar = numpy.zeros((nPts, 3))
        else:
            fBar = fBar
            
        # Determine if we can do a cheaper state-variable only
        # computation or if we need to include the spatial terms:
        useSpatial = False
        useState = False
        if wDeriv:
            useState = True
        if xDvDeriv or xVDeriv or xSDeriv or xDvDerivAero:
            useSpatial = True

        # Do actual Fortran call. 
        xvbar, extrabar, wbar = self.sumb.computematrixfreeproductbwd(
            resBar, funcsBar, fBar.T, useSpatial, useState, self.getSpatialSize(), 
            max(1, self.nDVAero))
     
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
                xsbar = self.mesh.getdXs(self.groupName)

                if xSDeriv:
                    returns.append(xsbar)

                # Process all the way back to the DVs, provided
                # DVGeo exists
                if xDvDeriv: 
                    xdvbar = {}
                    if self.DVGeo is not None and self.DVGeo.getNDV() > 0:
                        xdvbar.update(self.DVGeo.totalSensitivity(
                            xsbar, self.curAP.ptSetName, self.comm, config=self.curAP.name))
                    else:
                        if self.comm.rank == 0:
                            SUMBWarning("No DVGeo object is present or it has no "
                                        "design variables specified. No geometric "
                                        "derivatives computed.")

                    # Include aero derivatives here:
                    xdvbar.update(self._processAeroDerivatives(extrabar))
                    returns.append(xdvbar)
            else:
                raise Error("Could not complete requested xDvDeriv or xSDeriv "
                            "derivatives since no mesh is present")
                
        # Include the aerodynamic variables if requested to do so.
        if xDvDerivAero:
            xdvaerobar = {}
            xdvaerobar.update(self._processAeroDerivatives(extrabar))
            returns.append(xdvaerobar)

        # Single return (most frequent) is 'clean', otherwise a tuple.
        return tuple(returns) if len(returns) > 1 else returns[0]

    def sectionVectorByFamily(self, *args, **kwargs):
        return self.mesh.sectionVectorByFamily(*args, **kwargs)

    def expandVectorByFamily(self, *args, **kwargs):
        return self.mesh.expandVectorByFamily(*args, **kwargs)

    def finalizeAdjoint(self):
        """
        destroy the PESTc memory
        """
        self.releaseAdjointMemory()
        self.sumb.releasememadjoint()

    def getStateSize(self):
        """Return the number of degrees of freedom (states) that are
        on this processor. This is (number of states)*(number of
        cells)*(number of spectral instacnes)"""

        nstate = self.sumb.flowvarrefstate.nw
        ncells = self.sumb.adjointvars.ncellslocal[0]
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return nstate*ncells*ntime

    def getAdjointStateSize(self):
        """Return the number of ADJOINT degrees of freedom (states)
        that are on this processor. The reason this is different from
        getStateSize() is that if frozenTurbulence is used for RANS,
        the nonlinear system has 5+neq turb states per cell, while the
        adjoint still has 5."""
        if self.getOption('frozenTurbulence'):
            nstate = self.sumb.flowvarrefstate.nwf
        else:
            nstate = self.sumb.flowvarrefstate.nw

        ncells = self.sumb.adjointvars.ncellslocal[0]
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return nstate*ncells*ntime

    def getSpatialSize(self):
        """Return the number of degrees of spatial degrees of freedom
        on this processor. This is (number of nodes)*(number of
        spectral instances)*3"""

        nnodes = self.sumb.adjointvars.nnodeslocal[0]
        ntime  = self.sumb.inputtimespectral.ntimeintervalsspectral

        return 3*nnodes*ntime

    def getStates(self):
        """Return the states on this processor. Used in aerostructural
        analysis"""

        return self.sumb.getstates(self.getStateSize())

    def setStates(self, states):
        """Set the states on this processor. Used in aerostructural analysis
        and for switching aeroproblems
        """
        self.sumb.setstates(states)

    def _getInfo(self):
        """Get the haloed state vector, pressure (and viscocities). Used to
        save "state" between aeroProblems

        """
        return self.sumb.getinfo(self.sumb.getinfosize())

    def _setInfo(self, info):
        """Set the haloed state vector, pressure (and viscocities). Used to
        restore "state" between aeroProblems
        """
        self.sumb.setinfo(info)

    def setAdjoint(self, adjoint, objective=None):
        """Sets the adjoint vector externally. Used in coupled solver"""
        if objective is not None:
            self.curAP.sumbData.adjoints[objective] = adjoint.copy()

    def getAdjoint(self, objective):
        """ Return the adjoint values for objective if they
        exist. Otherwise just return zeros"""

        if objective in self.curAP.sumbData.adjoints:
            return self.curAP.sumbData.adjoints[objective]
        else:
            return numpy.zeros(self.getAdjointStateSize(), self.dtype)

    def getResidual(self, aeroProblem, res=None, releaseAdjointMemory=True):
        """Return the residual on this processor. Used in aerostructural
        analysis"""
        self.setAeroProblem(aeroProblem, releaseAdjointMemory)
        if res is None:
            res = numpy.zeros(self.getStateSize())
        res = self.sumb.getres(res)

        return res

    def computeArea(self, aeroProblem, funcs, axis, groupName=None, TS=0):
        """
        Compute the projected area of the surface mesh

        Input Arguments:
           axis, numpy array, size(3): The projection vector
               along which to determine the shadow area
           groupName, str: The group from which to obtain the coordinates.
               This name must have been obtained from addFamilyGroup() or
               be the default 'all' which contains all surface coordiantes

        Output Arguments:
            Area: The resulting area
            """
        self.setAeroProblem(aeroProblem)
        cfdForcePts = self.getForcePoints(TS)
        if len(cfdForcePts) > 0:
            areas = self.sumb.getareas(cfdForcePts.T, TS+1, axis).T
        else:
            areas = numpy.zeros((0,3), self.dtype)
        # end if

        if groupName is not None:
            areas = self.mesh.sectionVectorByFamily(groupName, areas)
        # end if

        # Now we do an mpiallreduce with sum:
        area = self.comm.allreduce(numpy.sum(areas), op=MPI.SUM)

        key = self.curAP.name + '_area'
        self.curAP.funcNames['area'] = key
        funcs[key] = area

    def computeAreaSensitivity(self, aeroProblem, funcsSens, axis, groupName=None, TS=0):
        """
        Compute the projected area of the surface mesh

        Input Arguments:
           axis, numpy array, size(3): The projection vector
               along which to determine the shadow area
           groupName, str: The group from which to obtain the coordinates.
               This name must have been obtained from addFamilyGroup() or
               be the default 'all' which contains all surface coordiantes

        Output Arguments:
            Area: The resulting area
            """
        self.setAeroProblem(aeroProblem)
        cfdForcePts = self.getForcePoints(TS)

        if len(cfdForcePts) > 0:
            da = self.sumb.getareasensitivity(cfdForcePts.T, TS+1, axis).T
        else:
            da = numpy.zeros((0,3), self.dtype)

        if groupName is not None:
            da = self.mesh.sectionVectorByFamily(groupName, da)
            da = self.mesh.expandVectorByFamily(groupName, da)

        funcsSens[self.curAP.name + '_area'] = self.DVGeo.totalSensitivity(
            da, ptSetName=self.curAP.ptSetName, comm=self.comm, config=self.curAP.name)

    def setOption(self, name, value):
        """
        Set Solver Option Value
        """
        name = name.lower()

        # Check to see if we have a deprecated option. Print a useful
        # warning that this is deprecated.
        if name in self.deprecatedOptions.keys():
            if self.comm.rank == 0:
                SUMBWarning("Option '%-29s\' is a deprecated SUmb Option |"% name)
            return

        # Try the option in the option dictionary to make sure we are setting a valid option
        defOptions = self.options['defaults']
        try:
            defOptions[name]
        except:
            if self.comm.rank == 0:
                SUMBWarning("Option '%-30s' is not a valid SUmb Option |"%name)
            return

        # Now we know the option exists, lets check if the type is ok:
        if isinstance(value, self.options[name][0]):
            # Just set:
            self.options[name] = [type(value),value]
        else:
            raise Error("Datatype for Option %-35s was not valid \n "
                        "Expected data type is %-47s \n "
                        "Received data type is %-47s"% (
                            name, self.options[name][0], type(value)))

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
                    self.sumb.monitorvariables(varStr)
                if name == 'surfacevariables':
                    self.sumb.surfacevariables(varStr)
                if name == 'volumevariables':
                    self.sumb.volumevariables(varStr)
                if name == 'isovariables':
                    self.sumb.isovariables(varStr)

            if name == 'isosurface':
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
                    # end for
                # end for
                val = numpy.array(val)

                self.sumb.initializeisosurfacevariables(val)
                for i in xrange(len(val)):
                    self.sumb.setisosurfacevariable(var[i], i+1)

            # end if
            if name == 'metricconversion':
                self.sumb.flowvarrefstate.lref = value
                self.sumb.flowvarrefstate.lrefspecified = True
                self.metricConversion = value
            # end if
            

            if name == "turbresscale":
                # If value is None no value has been specified by the user. None is the default value.
                # Do nothing as it will be updated with _updateTurbResScale from __init__                
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
                                raise Error("Datatype for Option %-35s was not valid. Expected list of <type 'float'>. Received data type is %-47s"% (name, type(value[0])))
                        else:
                            raise Error("Option %-35s of %-35s contains %-35s elements. Min and max number of elements are 1 and 4 respectively"% (name, type(value), len(value)))
                    else:
                        raise Error("Datatype for Option %-35s not valid. Expected data type is <type 'float'> or <type 'list'>. Received data type is %-47s"% (name, type(value)))

                    execStr = 'self.sumb.'+self.optionMap[name]['location'] + '=' + str(tmp_turbresscalar)
                    exec(execStr)   

            return
        # end if

        # All other options do genericaly by setting value in module:
        # Check if there is an additional mapping to what actually
        # has to be set in the solver

        temp = copy.copy(self.optionMap[name]) # This is the dictionary
        temp.pop('location')
        try:
            temp.pop('len')
        except:
            pass

        # If temp has anything left in it, we MUST be able to match to
        # one of them.

        if len(temp) == 0:
            pass
        else:
            #Convert the value to lower case:
            value = self.optionMap[name][value.lower()]

        # If value is a string, put quotes around it and make it
        # the correct length, otherwise convert to string
        if isinstance(value, str):
            spacesToAdd = self.optionMap[name]['len'] - len(value)
            value = '\'' + value + ' '*spacesToAdd + '\''
        else:
            value = str(value)
        # end if

        # Exec str is what is actually executed:
        execStr = 'self.sumb.'+self.optionMap[name]['location'] + '=' + value

        exec(execStr)

    def getOption(self, name):
        # Redefine the getOption def from the base class so we can
        # make sure the name is lowercase

        if name.lower() in self.options['defaults']:
            return self.options[name.lower()][1]
        else:
            raise Error('%s is not a valid option name'% name)

    def _getDefOptions(self):
        """
        There are many options for SUmb. These technically belong in
        the __init__ function but it gets far too long so we split
        them out.
        """
        defOpts = {
            # Log File instead of stdout
            'logfile':[str, ''],

            # Common Paramters
            'gridfile':[str, 'default.cgns'],
            'restartfile':[str, 'default_restart.cgns'],
            'solrestart':[bool, False],

            # Output Parameters
            'storerindlayer':[bool, True],
            'outputdirectory':[str, './'],
            'writesymmetry':[bool, True],
            'writefarfield':[bool, False],
            'writesurfacesolution':[bool,True],
            'writevolumesolution':[bool,True],
            'nsavevolume':[int,1],
            'nsavesurface':[int,1],
            'solutionprecision':[str,'single'],
            'gridprecision':[str,'double'],
            'isosurface':[dict, {}],
            'isovariables':[list, []],
            'viscoussurfacevelocities':[bool, True],
            'slicefiletractions':[bool, False],

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
            'usewallfunctions':[bool, False],
            'useapproxwalldistance':[bool, True],
            'walltreatment':[str, 'linear pressure extrapolation'],
            'dissipationscalingexponent':[float, 0.67],
            'vis4':[float, 0.0156],
            'vis2':[float, 0.25],
            'vis2coarse':[float, 0.5],
            'restrictionrelaxation':[float, .80],
            'liftindex':[int, 2],
            'lowspeedpreconditioner':[bool, False],

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

            # Unsteady Paramters
            'timeintegrationscheme':[str, 'bdf'],
            'timeaccuracy':[int, 2],
            'ntimestepscoarse':[int, 48],
            'ntimestepsfine':[int, 400],
            'deltat':[float, .010],
            'useale':[bool, True],

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
            'coeffconvcheck':[bool, False],
            'miniterationnum':[int, 0],

            # Newton-Krylov Paramters
            'usenksolver':[bool, False],
            'nklinearsolver':[str, 'gmres'],
            'nkswitchtol':[float, 2.5e-4],
            'nksubspacesize':[int, 60],
            'nklinearsolvetol':[float, 0.3],
            'nkuseew':[bool, True],
            'nkpc':[str, 'additive schwartz'],
            'nkadpc':[bool, False],
            'nkviscpc':[bool, False],
            'nkasmoverlap':[int, 1],
            'nkpcilufill':[int, 2],
            'nklocalpcordering':[str, 'rcm'],
            'nkjacobianlag':[int, 20],
            'rkreset':[bool, False],
            'nrkreset':[int, 5],
            'applypcsubspacesize':[int, 10],
            'nkinnerpreconits':[int, 1],
            'nkouterpreconits':[int, 1],
            'nkls':[str, 'cubic'],
            
            # Load Balance/partitioning parameters
            'blocksplitting':[bool, True],
            'loadimbalance':[float, 0.1],
            'loadbalanceiter':[int, 10],
            'partitiononly':[bool, False],

            # Misc Paramters
            'metricconversion':[float, 1.0],
            'autosolveretry':[bool, False],
            'autoadjointretry':[bool, False],
            'storehistory':[bool, False],
            'numbersolutions':[bool, True],
            'printiterations':[bool, True],
            'printtiming':[bool, True],
            'setmonitor':[bool, True],
            'printwarnings':[bool, True],
            'monitorvariables':[list, ['cpu','resrho','cl', 'cd']],
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
            'usereversemodead':[bool, False],
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

    def _getOptionMap(self):
        """ The SUmb option map"""

        optionMap = {
            # Common Paramters
            'gridfile':{'location':'inputio.gridfile',
                        'len':self.sumb.constants.maxstringlen},
            'restartfile':{'location':'inputio.restartfile',
                           'len':self.sumb.constants.maxstringlen},
            'solrestart':{'location':'inputio.restart'},
            'storerindlayer':{'location':'inputio.storerindlayer'},
            'writesymmetry':{'location':'inputio.writesymmetry'},
            'writefarfield':{'location':'inputio.writefarfield'},
            'slicefiletractions':{'location':'inputio.slicefiletractions'},
            'nsavevolume':{'location':'inputiteration.nsavevolume'},
            'nsavesurface':{'location':'inputiteration.nsavesurface'},
            'viscoussurfacevelocities':{'location':'inputio.viscoussurfacevelocities'},
            'solutionprecision':{'single':self.sumb.inputio.precisionsingle,
                                 'double':self.sumb.inputio.precisiondouble,
                                 'location':'inputio.precisionsol'},
            'gridprecision':{'single':self.sumb.inputio.precisionsingle,
                             'double':self.sumb.inputio.precisiondouble,
                             'location':'inputio.precisiongrid'},

            # Physics Paramters
            'discretization':{'central plus scalar dissipation':
                              self.sumb.inputdiscretization.dissscalar,
                              'central plus matrix dissipation':
                              self.sumb.inputdiscretization.dissmatrix,
                              'central plus cusp dissipation':
                              self.sumb.inputdiscretization.disscusp,
                              'upwind':
                              self.sumb.inputdiscretization.upwind,
                              'location':'inputdiscretization.spacediscr'},
            'coarsediscretization':{'central plus scalar dissipation':
                                    self.sumb.inputdiscretization.dissscalar,
                                    'central plus matrix dissipation':
                                    self.sumb.inputdiscretization.dissmatrix,
                                    'central plus cusp dissipation':
                                    self.sumb.inputdiscretization.disscusp,
                                    'upwind':
                                    self.sumb.inputdiscretization.upwind,
                                    'location':'inputdiscretization.spacediscrcoarse'},
            'limiter':{'vanalbeda':self.sumb.inputdiscretization.vanalbeda,
                       'minmod':self.sumb.inputdiscretization.minmod,
                       'nolimiter':self.sumb.inputdiscretization.nolimiter,
                       'location':'inputdiscretization.limiter'},
            'smoother':{'runge kutta':self.sumb.inputiteration.rungekutta,
                        'lu sgs':self.sumb.inputiteration.nllusgs,
                        'lu sgs line':self.sumb.inputiteration.nllusgsline,
                        'dadi':self.sumb.inputiteration.dadi,
                        'location':'inputiteration.smoother'},

            'equationtype':{'euler':self.sumb.inputphysics.eulerequations,
                            'laminar ns':self.sumb.inputphysics.nsequations,
                            'rans':self.sumb.inputphysics.ransequations,
                            'location':'inputphysics.equations'},
            'equationmode':{'steady':self.sumb.inputphysics.steady,
                            'unsteady':self.sumb.inputphysics.unsteady,
                            'time spectral':self.sumb.inputphysics.timespectral,
                            'location':'inputphysics.equationmode'},
            'flowtype':{'internal':self.sumb.inputphysics.internalflow,
                        'external':self.sumb.inputphysics.externalflow,
                        'location':'inputphysics.flowtype'},
            'turbulencemodel':{'baldwin lomax':self.sumb.inputphysics.baldwinlomax,
                               'sa':self.sumb.inputphysics.spalartallmaras,
                               'sae':self.sumb.inputphysics.spalartallmarasedwards,
                               'k omega wilcox':self.sumb.inputphysics.komegawilcox,
                               'k omega modified':self.sumb.inputphysics.komegamodified,
                               'ktau':self.sumb.inputphysics.ktau,
                               'menter sst':self.sumb.inputphysics.mentersst,
                               'v2f':self.sumb.inputphysics.v2f,
                               'location':'inputphysics.turbmodel'},
            'turbulenceorder':{'first order':1,
                               'second order':2,
                               'location':'inputdiscretization.orderturb'},
            'turbresscale':{'location':'inputiteration.turbresscale'},
            'usewallfunctions':{'location':'inputphysics.wallfunctions'},
            'useapproxwalldistance':{'location':'inputdiscretization.useapproxwalldistance'},
                    'reynoldsnumber':{'location':'inputphysics.reynolds'},
            'reynoldslength':{'location':'inputphysics.reynoldslength'},
            'walltreatment':{'linear pressure extrapolation':
                             self.sumb.inputdiscretization.linextrapolpressure,
                             'constant pressure extrapolation':
                             self.sumb.inputdiscretization.constantpressure,
                             'quadratic pressure extrapolation':
                             self.sumb.inputdiscretization.quadextrapolpressure,
                             'normal momentum':
                             self.sumb.inputdiscretization.normalmomentum,
                             'location':'inputdiscretization.wallbctreatment'},
            'dissipationscalingexponent':{'location':'inputdiscretization.adis'},
            'vis4':{'location':'inputdiscretization.vis4'},
            'vis2':{'location':'inputdiscretization.vis2'},
            'vis2coarse':{'location':'inputdiscretization.vis2coarse'},
            'restrictionrelaxation':{'location':'inputiteration.fcoll'},
            'forcesastractions':{'location':'inputphysics.forcesastractions'},
            'lowspeedpreconditioner':{'location':'inputdiscretization.lowspeedpreconditioner'},

            # Common Paramters
            'ncycles':{'location':'inputiteration.ncycles'},
            'ncyclescoarse':{'location':'inputiteration.ncyclescoarse'},
            'nsubiterturb':{'location':'inputiteration.nsubiterturb'},
            'nsubiter':{'location':'inputiteration.nsubiterations'},
            'cfl':{'location':'inputiteration.cfl'},
            'cflcoarse':{'location':'inputiteration.cflcoarse'},
            'mgcycle':{'location':'localmg.mgdescription',
                       'len':self.sumb.constants.maxstringlen},
            'mgstartlevel':{'location':'inputiteration.mgstartlevel'},
            'resaveraging':{'noresaveraging':self.sumb.inputiteration.noresaveraging,
                            'alwaysresaveraging':self.sumb.inputiteration.alwaysresaveraging,
                            'alternateresaveraging':self.sumb.inputiteration.alternateresaveraging,
                            'location':'inputiteration.resaveraging'},
            'smoothparameter':{'location':'inputiteration.smoop'},
            'cfllimit':{'location':'inputiteration.cfllimit'},

            # Unsteady Params
            'timeintegrationscheme':{'bdf':self.sumb.inputunsteady.bdf,
                                     'explicitrk':self.sumb.inputunsteady.explicitrk,
                                     'implicitrk':self.sumb.inputunsteady.implicitrk,
                                     'md':self.sumb.inputunsteady.md,
                                     'location':'inputunsteady.timeintegrationscheme'},
            'timeaccuracy':{'location':'inputunsteady.timeaccuracy'},
            'ntimestepscoarse':{'location':'inputunsteady.ntimestepscoarse'},
            'ntimestepsfine':{'location':'inputunsteady.ntimestepsfine'},
            'deltat':{'location':'inputunsteady.deltat'},
            'useale':{'location':'inputunsteady.useale'},
            # Time Spectral Paramters
            'timeintervals':{'location':'inputtimespectral.ntimeintervalsspectral'},
            'alphamode':{'location':'inputtsstabderiv.tsalphamode'},
            'betamode':{'location':'inputtsstabderiv.tsbetamode'},
            'machmode':{'location':'inputtsstabderiv.tsmachmode'},
            'pmode':{'location':'inputtsstabderiv.tspmode'},
            'qmode':{'location':'inputtsstabderiv.tsqmode'},
            'rmode':{'location':'inputtsstabderiv.tsrmode'},
            'altitudemode':{'location':'inputtsstabderiv.tsaltitudemode'},
            'windaxis':{'location':'inputtsstabderiv.usewindaxis'},
            'alphafollowing':{'location':'inputtsstabderiv.tsalphafollowing'}, 
            'tsstability':{'location':'inputtsstabderiv.tsstability'},

            # Convergence Paramters
            'l2convergence':{'location':'inputiteration.l2conv'},
            'l2convergencerel':{'location':'inputiteration.l2convrel'},
            'l2convergencecoarse':{'location':'inputiteration.l2convcoarse'},
            'maxl2deviationfactor':{'location':'inputiteration.maxl2deviationfactor'},
            'coeffconvcheck':{'location':'monitor.coeffconvcheck'},
            'miniterationnum':{'location':'inputiteration.miniternum'},

            # Newton-Krylov Paramters
            'usenksolver':{'location':'nksolvervars.usenksolver'},
            'nklinearsolver':{'gmres':'gmres',
                              'tfqmr':'tfqmr',
                              'location':'nksolvervars.ksp_solver_type',
                              'len':self.sumb.constants.maxstringlen},
            'nkuseew':{'location':'nksolvervars.nkuseew'},
            'nkswitchtol':{'location':'nksolvervars.nk_switch_tol'},
            'nksubspacesize':{'location':'nksolvervars.ksp_subspace'},
            'nklinearsolvetol':{'location':'nksolvervars.ksp_rtol'},
            'nkpc':{'additive schwartz':'asm',
                    'multigrid':'mg',
                    'location':'nksolvervars.global_pc_type',
                    'len':self.sumb.constants.maxstringlen},
            'nkasmoverlap':{'location':'nksolvervars.asm_overlap'},
            'nkpcilufill':{'location':'nksolvervars.local_pc_ilu_level'},
            'nklocalpcordering':{'natural':'natural',
                                 'rcm':'rcm',
                                 'nested dissection':'nd',
                                 'one way dissection':'1wd',
                                 'quotient minimum degree':'qmd',
                                 'location':
                                 'nksolvervars.local_pc_ordering',
                                 'len':self.sumb.constants.maxstringlen},
            'nkmaxlinearkspits':{'location':'nksolvervars.ksp_max_it'},
            'nkjacobianlag':{'location':'nksolvervars.jacobian_lag'},
            'rkreset':{'location':'nksolvervars.rkreset'},
            'nrkreset':{'location':'nksolvervars.nrkreset'},
            'nkadpc':{'location':'nksolvervars.nkadpc'},
            'nkviscpc':{'location':'nksolvervars.nkviscpc'},
            'applypcsubspacesize':{'location':'nksolvervars.applypcsubspacesize'},
            'nkinnerpreconits':{'location':'nksolvervars.innerpreconits'},
            'nkouterpreconits':{'location':'nksolvervars.outerpreconits'},
            'nkls':{'none':self.sumb.nksolvervars.nolinesearch,
                    'cubic':self.sumb.nksolvervars.cubiclinesearch,
                    'non monotone':self.sumb.nksolvervars.nonmonotonelinesearch,
                    'location':'nksolvervars.nkls'
                },

            # Load Balance Paramters
            'blocksplitting':{'location':'inputparallel.splitblocks'},
            'loadimbalance':{'location':'inputparallel.loadimbalance'},
            'loadbalanceiter':{'location':'inputparallel.loadbalanceiter'},

            # Misc Paramters
            'printiterations':{'location':'inputiteration.printiterations'},
            'printwarnings':{'location':'inputiteration.printwarnings'},
            'printtiming':{'location':'inputadjoint.printtiming'},
            'setmonitor':{'location':'inputadjoint.setmonitor'},

            # Adjoint Params
            'adjointl2convergence':{'location':'inputadjoint.adjreltol'},
            'adjointl2convergencerel':{'location':'inputadjoint.adjreltolrel'},
            'adjointl2convergenceabs':{'location':'inputadjoint.adjabstol'},
            'adjointdivtol':{'location':'inputadjoint.adjdivtol'},
            'approxpc':{'location':'inputadjoint.approxpc'},
            'adpc':{'location':'inputadjoint.adpc'},
            'viscpc':{'location':'inputadjoint.viscpc'},
            'frozenturbulence':{'location':'inputadjoint.frozenturbulence'},
            'usediagtspc':{'location':'inputadjoint.usediagtspc'},
            'restartadjoint':{'location':'inputadjoint.restartadjoint'},
            'adjointsolver':{'gmres':'gmres',
                             'tfqmr':'tfqmr',
                             'richardson':'richardson',
                             'bcgs':'bcgs',
                             'ibcgs':'ibcgs',
                             'location':'inputadjoint.adjointsolvertype',
                             'len':self.sumb.constants.maxstringlen},
            'adjointmaxiter':{'location':'inputadjoint.adjmaxiter'},
            'adjointsubspacesize':{'location':'inputadjoint.adjrestart'},
            'adjointmonitorstep':{'location':'inputadjoint.adjmonstep'},
            'dissipationlumpingparameter':{'location':'inputdiscretization.sigma'},
            'preconditionerside':{'left':'left',
                                  'right':'right',
                                  'location':'inputadjoint.adjointpcside',
                                  'len':self.sumb.constants.maxstringlen},
            'matrixordering':{'natural':'natural',
                              'rcm':'rcm',
                              'nested dissection':'nd',
                              'one way dissection':'1wd',
                              'quotient minimum degree':'qmd',
                              'location':
                              'inputadjoint.matrixordering',
                              'len':self.sumb.constants.maxstringlen},
            'globalpreconditioner':{'additive schwartz':'asm',
                                    'multigrid':'mg',
                                    'location':'inputadjoint.precondtype',
                                    'len':self.sumb.constants.maxstringlen},
            'localpreconditioner':{'ilu':'ilu',
                                   'location':'inputadjoint.localpctype',
                                   'len':self.sumb.constants.maxstringlen},
            'ilufill':{'location':'inputadjoint.filllevel'},
            'applyadjointpcsubspacesize':{
                'location':'inputadjoint.applyadjointpcsubspacesize'},
            'asmoverlap':{'location':'inputadjoint.overlap'},
            'innerpreconits':{'location':'inputadjoint.innerpreconits'},
            'outerpreconits':{'location':'inputadjoint.outerpreconits'},
            'firstrun':{'location':'inputadjoint.firstrun'},
            'verifystate':{'location':'inputadjoint.verifystate'},
            'verifyspatial':{'location':'inputadjoint.verifyspatial'},
            'verifyextra':{'location':'inputadjoint.verifyextra'},
            'usematrixfreedrdw':{'location':'inputadjoint.usematrixfreedrdw'},

            # Parameters for functions
            'sepsensoroffset':{'location':'costfunctions.sepsensoroffset'},
            'sepsensorsharpness':{'location':'costfunctions.sepsensorsharpness'},
            }

        return optionMap

    def _getSpecialOptionLists(self):
        """
        Lists of special options
        """

        # These "ignore_options" are NOT actually, ignored, rather,
        # they DO NOT GET SET IN THE FORTRAN CODE. Rather, they are
        # used strictly in Python

        ignoreOptions = [
            'defaults',
            'storehistory',
            'numbersolutions',
            'writesurfacesolution',
            'writevolumesolution',
            'autosolveretry',
            'autoadjointretry',
            'usereversemodead',
            'partitiononly',
            'liftindex',
            'logfile'
             ]

        # Deprecated options. These should not be used, but old
        # scripts can continue to run
        deprecatedOptions = {'finitedifferencepc':'Use the ADPC option.',
                             'writesolution':'Use writeSurfaceSolution and writeVolumeSolution options instead.',
                             'reynoldsnumber':'Put required information in aeroProblem class',
                             'reynoldslength':'Put required information in aeroProblem class'
                             }

        specialOptions = ['surfacevariables',
                          'volumevariables',
                          'monitorvariables',
                          'metricconversion',
                          'outputdirectory',
                          'isovariables',
                          'isosurface',
                          'turbresscale',
                          ]

        return ignoreOptions, deprecatedOptions, specialOptions

    def _getObjectivesAndDVs(self):

        possibleAeroDVs = {
            'alpha':'adjointvars.ndesignaoa',
            'beta':'adjointvars.ndesignssa',
            'mach':'adjointvars.ndesignmach',
            'machgrid':'adjointvars.ndesignmachgrid',
            'P':'adjointvars.ndesignpressure',
            'reynolds':'adjointvars.ndesignreynolds',
            'T':'adjointvars.ndesigntemperature',
            'rotX':'adjointvars.ndesignrotx',
            'rotY':'adjointvars.ndesignroty',
            'rotZ':'adjointvars.ndesignrotz',
            'rotcenx':'adjointvars.ndesignrotcenx',
            'rotceny':'adjointvars.ndesignrotceny',
            'rotcenz':'adjointvars.ndesignrotcenz',
            'xRef':'adjointvars.ndesignpointrefx',
            'yYef':'adjointvars.ndesignpointrefy',
            'zRef':'adjointvars.ndesignpointrefz',
            'chordRef':'adjointvars.ndesignlengthref',
            'areaRef':'adjointvars.ndesignsurfaceref',
            'disserror':'adjointvars.ndesigndisserror'
            }

        # This is SUmb's internal mapping for cost functions
        sumbCostFunctions = {
            'lift':self.sumb.costfunctions.costfunclift,
            'drag':self.sumb.costfunctions.costfuncdrag,
            'cl'  :self.sumb.costfunctions.costfuncliftcoef,
            'cd'  :self.sumb.costfunctions.costfuncdragcoef,
            'fx'  :self.sumb.costfunctions.costfuncforcex,
            'fy'  :self.sumb.costfunctions.costfuncforcey,
            'fz'  :self.sumb.costfunctions.costfuncforcez,
            'cfx' :self.sumb.costfunctions.costfuncforcexcoef,
            'cfy' :self.sumb.costfunctions.costfuncforceycoef,
            'cfz' :self.sumb.costfunctions.costfuncforcezcoef,
            'mx'  :self.sumb.costfunctions.costfuncmomx,
            'my'  :self.sumb.costfunctions.costfuncmomy,
            'mz'  :self.sumb.costfunctions.costfuncmomz,
            'cmx':self.sumb.costfunctions.costfuncmomxcoef,
            'cmy':self.sumb.costfunctions.costfuncmomycoef,
            'cmz':self.sumb.costfunctions.costfuncmomzcoef,
            'cm0':self.sumb.costfunctions.costfunccm0,
            'cmzalpha':self.sumb.costfunctions.costfunccmzalpha,
            'cmzalphadot':self.sumb.costfunctions.costfunccmzalphadot,
            'cl0':self.sumb.costfunctions.costfunccl0,
            'clalpha':self.sumb.costfunctions.costfuncclalpha,
            'clalphadot':self.sumb.costfunctions.costfuncclalphadot,
            'cfy0':self.sumb.costfunctions.costfunccfy0,
            'cfyalpha':self.sumb.costfunctions.costfunccfyalpha,
            'cfyalphdDot':self.sumb.costfunctions.costfunccfyalphadot,
            'cd0':self.sumb.costfunctions.costfunccd0,
            'cdalpha':self.sumb.costfunctions.costfunccdalpha,
            'cdalphadot':self.sumb.costfunctions.costfunccdalphadot,
            'cmzq':self.sumb.costfunctions.costfunccmzq,
            'cmzqdot':self.sumb.costfunctions.costfunccmzqdot,
            'clq':self.sumb.costfunctions.costfuncclq,
            'clqdot':self.sumb.costfunctions.costfuncclqdot,
            'cbend':self.sumb.costfunctions.costfuncbendingcoef,
            'sepsensor':self.sumb.costfunctions.costfuncsepsensor,
            'sepsensoravgx':self.sumb.costfunctions.costfuncsepsensoravgx,
            'sepsensoravgy':self.sumb.costfunctions.costfuncsepsensoravgy,
            'sepsensoravgz':self.sumb.costfunctions.costfuncsepsensoravgz,
            'cavitation':self.sumb.costfunctions.costfunccavitation,
            }

        return possibleAeroDVs, sumbCostFunctions
        
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


    def __del__(self):
        if self.logFile:
            self.sumb.closeLog()

class sumbFlowCase(object):
    """
    Class containing the data that SUmb requires to be saved to an
    aeroProblem to permit the analysis of multiple flow cases
    """
    def __init__(self):
        self.stateInfo = None
        self.adjoints = {}
        self.adjointRHS = {}
        self.coords = None
        self.callCounter = -1
        self.disp = None


