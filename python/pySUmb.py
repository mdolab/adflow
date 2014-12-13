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
        i = 18
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
        self.possibleObjectives, self.possibleAeroDVs, self.sumbCostFunctions = \
                                 self._getObjectivesAndDVs()

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
        # geometric manipulation object
        self.mesh = None
        self.DVGeo = None

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

        if self.getOption('partitionOnly'):
            return

        # Finally complete loading
        self.sumb.dummyreadparamfile()
        self.sumb.partitionandreadgrid()
        self.sumb.preprocessing()
        self.sumb.initflow()
        self.sumb.preprocessingadjoint()

    def setMesh(self, mesh):
        """
        Set the mesh object to SUmb to do geometric deformations

        Parameters
        ----------
        mesh : multiBlockMesh object
            The pyWarp mesh object for doing the warping
        """

        self.mesh = mesh

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
        self.coords0 = self.getSurfaceCoordinates()

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
            The family (as defined in pyWarp) to use for the slices. Not
            currently supported.
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

        # If this problem has't been solved yet, reset flow to this
        # flight condition
        if self.curAP.sumbData.stateInfo is None:
            self.resetFlow(aeroProblem, releaseAdjointMemory)

        # Possibly release adjoint memory if not already done so.
        if releaseAdjointMemory:
            self.releaseAdjointMemory()

        # Save aeroProblem, and other information into the current flow case
        self.curAP.sumbData.adjointRHS = None
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
        if MDCallBack is None:
            self.sumb.solver()
        else:
            self.sumb.solverunsteadymd(MDCallBack)

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

        # Just call the regular getSolution() command and extract the
        # ones we need:
        res = self.getSolution(sps)

        for f in evalFuncs:
            if f in self.possibleObjectives:
                key = self.curAP.name + '_%s'% f
                self.curAP.funcNames[f] = key
                funcs[key] = res[f]
            else:
                if not ignoreMissing:
                    raise Error('Supplied function is not known to SUmb.')

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

        # 'Extra' derivatives --- ie the aerodynamic/reference
        # derivatives. We must provide the derivatives that have been
        # asked for in the aeroProblem. Before we start, lets setup
        # the aeroDVs we need
        self.setAeroProblem(aeroProblem)

        if evalFuncs is None:
            evalFuncs = self.curAP.evalFuncs
        else:
            evalFuncs = set(evalFuncs)

        # Do the functions one at a time:
        for f in evalFuncs:
            if f not in self.possibleObjectives:
                raise Error('Supplied function is not known to SUmb.')
            if self.comm.rank == 0:
                print('Solving adjoint: %s'%f)

            key = self.curAP.name + '_%s'% f
            ptSetName = self.curAP.ptSetName

            # Set dict structure for this derivative
            funcsSens[key] = {}

            # Solve adjoint equation (if necessary)
            self.solveAdjoint(aeroProblem, f)

            # Geometric derivatives
            if self.DVGeo is not None and self.DVGeo.getNDV() > 0:
                dIdpt = self.totalSurfaceDerivative(f)
                dIdx = self.DVGeo.totalSensitivity(
                    dIdpt, ptSetName=ptSetName, comm=self.comm, config=self.curAP.name)
                funcsSens[key].update(dIdx)

            # Compute total aero derivatives
            funcsSens[key].update(self.totalAeroDerivative(f))

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

    def solveTrimCL(self, aeroProblem, CLStar, trimStar, alpha0, tail0, 
                    trimValue, trimDV, dvIndex=0, tol=1e-3, nIter=10):
        """Solve the trim-Cl problem: Find the angle of attack and
        trim actuator value that satifies a given lift and moment
        value. It is solved using a Broyden method. A DVGeometry and
        mesh object must be set for this routine to be used. 

        Parameter
        ---------
        aeroProblem : pyAero_problem class
            The aerodynamic problem to solve
        CLStar : float
            The desired target CL
        trimStar : float
            The desired target moment value. Usually this is zero
        alpha0 : angle (deg)
            Initial guess for angle of attach 
        tail0 : angle (deg)
            Initial guess for tail-anlge design varaible 
        trimValue : str
            The SUmb function to use for moment. This will normally be either
            'cmy' or 'cmz' depending on the orientation. 
        trimDV : str
             DVGeo function name that is used for trimming
        dvIndex : int 
             Index of trimDV to use in DVGeo. By default this is 0 which means use
             use the first DV in the DVGeo design variable.
        tol : float
             Tolerance for broyden solution
        nIter : int
             Maximum number of iterations
             
        Returns
        -------
        None. However, the alpha value is updated in the aeroProblem and tail angle
        variable is updated in DVGeo.
           """

        self.setAeroProblem(aeroProblem)
        
        def broyden3(F, xin, tol, iter, alpha=0.4):
            """Broyden's second method. See Scipy Documetation --> This is
            copied verbatim. Modified to add tolerance"""
            zy = []
            def myF(F, xm):
                return numpy.matrix(F(tuple(xm.flat))).T

            def updateG(z, y):
                "G:=G+z*y.T"
                zy.append((z, y))

            def norm(v):
                """Returns an L2 norm of the vector."""
                return numpy.sqrt(numpy.sum((numpy.array(v)**2).flat))
        
            def Gmul(f):
                "G=-alpha*1+z*y.T+z*y.T ..."
                s = -alpha*f
                for z, y in zy:
                    s = s+z*(y.T*f)
                return s

            xm = numpy.matrix(xin).T
            Fxm = myF(F, xm)
            #    Gm=-alpha*numpy.matrix(numpy.identity(len(xin)))
            for n in range(iter):
                #deltaxm=-Gm*Fxm
                deltaxm = Gmul(-Fxm)
                xm = xm+deltaxm
                Fxm1 = myF(F, xm)
                deltaFxm = Fxm1 - Fxm
                Fxm = Fxm1
                #Gm=Gm+(deltaxm-Gm*deltaFxm)*deltaFxm.T/norm(deltaFxm)**2
                updateG(deltaxm-Gmul(deltaFxm), deltaFxm/norm(deltaFxm)**2)
                if numpy.linalg.norm(Fxm) < tol:
                    break

            return xm.flat

        def function(X):
            if self.comm.rank == 0:
                print ('X:',X)
            # Set Alpha
            self.curAP.alpha = X[0]

            # Set Trim Anlge
            xGeo = self.DVGeo.getValues()
            xGeo[trimDV][dvIndex] = X[1]
            self.DVGeo.setDesignVars(xGeo)

            # Solve Problem:
            self.__call__(aeroProblem, writeSolution=False)
    
            # Extract Solution
            sol = self.getSolution()

            # Build 'F'
            F = [sol['cl']-CLStar, sol[trimValue]-trimStar]

            return F
        
        # Run Broyden search
        broyden3(function, [alpha0, tail0], tol, nIter)

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
        self.sumb.inputiteration.groundlevel = strLvl
        self.sumb.inputiteration.currentlevel = strLvl
        self.sumb.monitor.niterold = 0
        self.sumb.monitor.nitercur = 0
        self.sumb.iteration.itertot = 0
        self.sumb.setuniformflow()
        self.sumb.nksolvervars.nksolvecount = 0
        self.sumb.killsignals.routinefailed =  False
        self.sumb.killsignals.fatalfail = False
        self.sumb.nksolvervars.freestreamresset = False
    def getSolution(self, sps=1, groupName=None):
        """ Retrieve the solution variables from the solver. Note this
        is a collective function and must be called on all processors
        """

        # Get the mask for the group
        mask = self._getCellGroupMaskLocal(groupName)
        
        # We should return the list of results that is the same as the
        # possibleObjectives list
        self.sumb.getsolutionmask(sps, mask)

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
                            "been added in pyWarp.")

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
            return numpy.empty((0, 3), self.dtype)

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
                if not ptSetName in self.DVGeo.points:
                    # DVGeo appeared and we have not embedded points!
                    self.DVGeo.addPointSet(self.coords0, ptSetName)
                if not self.DVGeo.pointSetUpToDate(ptSetName):
                    self.setSurfaceCoordinates(
                        self.DVGeo.update(ptSetName, config=self.curAP.name), 'all')
                    self.updateGeometryInfo()
            # Finally update other data
            self._setAeroProblemData()

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
            aeroProblem.surfMesh = self.getSurfaceCoordinates('all')

        if self.curAP is not None:
            # If we have already solved something and are now
            # switching, save what we need:
            self.curAP.stateInfo = self._getInfo()
            self.curAP.surfMesh = self.getSurfaceCoordinates('all')

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
                self.setSurfaceCoordinates(self.DVGeo.update(ptSetName, config=self.curAP.name), 'all')
            else:
                self.setSurfaceCoordinates(self.curAP.surfMesh, 'all')
        else:
            self.setSurfaceCoordinates(self.curAP.surfMesh, 'all')

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

        # Finally update other data
        self._setAeroProblemData()

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
        else:
            T = AP.T
            P = AP.P
            rho = AP.rho
            V = AP.V
            mu = AP.mu

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

        #update the flow vars
        if not firstCall:
            self.sumb.updateflow()

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
            self.sumb.flowvarrefstate.pref = P
            self.sumb.flowvarrefstate.tref = T
            self.sumb.inputphysics.tempfreestream = T
            ReLength = 1.0
            self.sumb.inputphysics.reynolds = rho*V/mu
            self.sumb.inputphysics.reynoldslength = ReLength
        else:
            self.sumb.flowvarrefstate.pref = P
            self.sumb.flowvarrefstate.tref = T
            self.sumb.inputphysics.tempfreestream = T
            self.sumb.inputphysics.reynolds = 1.0
            self.sumb.inputphysics.reynoldslength = 1.0

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
        # end if

        if not firstCall:
            self.sumb.referencestate()
            self.sumb.setflowinfinitystate()
            self.sumb.updateperiodicinfoalllevels()
            self.sumb.updategridvelocitiesalllevels()

    def getForces(self, groupName=None, TS=0, pressure=True, viscous=True):
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

        [npts, ncell] = self.sumb.getforcesize()
        if npts > 0:
            forcesp, forcesv = self.sumb.getforces(npts, TS+1)
            forcesp = forcesp.T
            forcesv = forcesv.T

            forces = numpy.zeros_like(forcesp)
            if pressure:
                forces += forcesp
            if viscous:
                forces += forcesv
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
        # For now, just create all the petsc variables
        if not self.adjointSetup or reform:
            self.releaseAdjointMemory()
            self.sumb.createpetscvars()

            if self.getOption('useReverseModeAD'):
                self.sumb.setupallresidualmatrices()
            else:
                self.sumb.setupallresidualmatricesfwd()

            # Create coupling matrix struct whether we need it or not
            [npts, ncells] = self.sumb.getforcesize()
            nTS = self.sumb.inputtimespectral.ntimeintervalsspectral
            forcePoints = numpy.zeros((nTS, npts, 3), self.dtype)
            for i in xrange(nTS):
                forcePoints[i] = self.getForcePoints(TS=i)

            self.sumb.setupcouplingmatrixstruct(forcePoints.T)

            # Setup the KSP object
            self.sumb.setuppetscksp()

            # Set the flag
            self.adjointSetup = True

    def printMatrixInfo(self, dRdwT=False, dRdwPre=False, dRdx=False,
                        dRda=False, dSdw=False, dSdx=False,
                        printLocal=False, printSum=False, printMax=False):

        # Call sumb matrixinfo function
        self.sumb.matrixinfo(dRdwT, dRdwPre, dRdx, dRda, dSdw, dSdx,
                             printLocal, printSum, printMax)

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

        obj, aeroObj = self._getObjective(objective)

        # Setup adjoint matrices/vector as required
        self._setupAdjoint()

        # Check to see if the RHS Partials have been computed
        if not self.curAP.adjointRHS == obj:
            self.computeObjPartials(objective, forcePoints)

        # Check to see if we need to agument the RHS with a structural
        # adjoint:
        if structAdjoint is not None and groupName is not None:
            if self.getOption('usereversemodead'):
                raise Error("Reverse mode AD no longer supported with "
                "aerostructural analysis. Use Forward mode AD for the adjoint")

            phi = self.mesh.expandVectorByFamily(groupName, structAdjoint)
            self.sumb.agumentrhs(numpy.ravel(phi))

        # Check if objective is allocated:
        if obj not in self.curAP.sumbData.adjoints.keys():
            self.curAP.sumbData.adjoints[obj] = numpy.zeros(self.getAdjointStateSize(), float)
        self.sumb.setadjoint(self.curAP.sumbData.adjoints[obj])

        # Actually Solve the adjoint system
        self.sumb.solveadjointtransposepetsc()

        # Possibly try another solve
        if self.sumb.killsignals.adjointfailed and self.getOption('restartAdjoint'):
            # Only retry if the following conditions are met:

            # 1. restartAdjoint is true -> that is we were starting
            # from a non-zero starting point

            # 2. The stored adjoint must have been already set at
            # least once; that is we've already tried one solve

            self.curAP.sumbData.adjoints[obj][:] = 0.0
            if self.getOption('autoAdjointRetry'):
                self.sumb.solveadjointtransposepetsc()

        # Now set the flags and possibly reset adjoint
        if self.sumb.killsignals.adjointfailed == False:
            self.curAP.sumbData.adjoints[obj] = \
                self.sumb.getadjoint(self.getAdjointStateSize())
            self.adjointFailed = False
        else:
            self.adjointFailed = True

            # Reset stored adjoint
            self.curAP.sumbData.adjoints[obj][:] = 0.0

    def totalSurfaceDerivative(self, objective):
        # The adjoint vector is now calculated so perform the
        # following operation to produce dI/dX_surf:
        # (p represents partial, d total)
        # dI/dX_s = pI/pX_s - (dXv/dXs)^T * ( dRdX_v^T * psi)
        #
        # The derivative wrt the surface captures the effect of ALL
        # GLOBAL Multidisciplinary variables -- any DV that changes
        # the surface.

        # NOTE: do dRdxvPsi MUST be done first since this
        # allocates spatial memory if required.
        dIdxs_2 = self.getdRdXvTPsi(objective, 'all')

        # Direct partial derivative contibution
        dIdxs_1 = self.getdIdx(objective, groupName='all')

        # Total derivative of the obective with surface coordinates
        dIdXs = dIdxs_1 - dIdxs_2

        return dIdXs

    def totalAeroDerivative(self, obj, extraSens=None):
        """
        This function returns the total derivative of the obj with
        respect to the aerodynamic variables defined the currently set
        aeroproblem. It essentially processes the output of
        _totalAeroDerivative to the dictionary return format

        Parameters
        ----------
        objectives : list, set
            The list of objectives to get derivatives for

        Returns
        -------
        funcsSens : dict
            The dictionary of the derivatives of obj wrt the aerodynamic
            variables
            """

        # Get the list of derivatives and then process:
        res = self._totalAeroDerivative(obj)
        funcsSens = {}

        DVsRequired = list(self.curAP.DVNames.keys())
        for dv in DVsRequired:
            if dv in self.possibleAeroDVs:
                funcsSens[self.curAP.DVNames[dv]] = res[self.aeroDVs[dv]]
            elif dv in ['altitude']:
                # Extract the derivatives wrt the independent
                # parameters in SUmb
                dIdP = res[self.aeroDVs['P']]
                dIdT = res[self.aeroDVs['T']]

                # Get the derivative of THESE values with respect
                # to the altitude:
                tmp = {}
                self.curAP.evalFunctionsSens(tmp, ['P', 'T'])
                # Chain-rule to get the final derivative:
                funcsSens[self.curAP.DVNames[dv]] = (
                    tmp[self.curAP['P']][self.curAP.DVNames['altitude']]*dIdP +
                    tmp[self.curAP['T']][self.curAP.DVNames['altitude']]*dIdT)
            # end if (dv type)
        # end for (dv loop)

        return funcsSens

    def _setAeroDVs(self):
        """ Do everything that is required to deal with aerodynamic
        design variables in SUmb"""

        DVsRequired = list(self.curAP.DVNames.keys())
        DVMap = {}
        for dv in DVsRequired:
            if dv in self.possibleAeroDVs:
                self._addAeroDV(dv)
            elif dv.lower() in ['altitude']:
                # We internally need to add the following:
                self._addAeroDV('T')
                self._addAeroDV('P')
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

        # Run the small amount of fortran code for adjoint initialization
        self.sumb.preprocessingadjoint()

    def _totalAeroDerivative(self, objective):
        """
        This function computes the total aerodynamic derivatives for
        the given objective. It returns a flattened array of the
        aerodynamic variables which may more may not be actually the
        DV's that the user wants.

        The adjoint vector is now calculated. This function computes
        dI/dX_aero = pI/pX_aero - dR/dX_aero^T * psi. The "aero"
        variables are intrinsic ONLY to the aero discipline. Nothing
        in the structural process should depend on these functions
        directly.
        """

        obj, aeroObj = self._getObjective(objective)

        if obj in self.curAP.sumbData.adjoints:
            psi = self.curAP.sumbData.adjoints[obj]
        else:
            raise Error("%s adjoint for current aeroProblem "
                        "%s is not computed."% obj)

        # Direct partial derivative contibution
        dIda_1 = self.getdIda(objective)

        # dIda contribution for drda^T * psi
        dIda_2 = self.getdRdaPsi(psi)

        # Total derivative of the obective wrt aero-only DVs
        dIda = dIda_1 - dIda_2

        # Alpha derivative needs to be scaled to degrees to be
        # consistent:
        if self.sumb.adjointvars.ndesignaoa != -1:
            dIda[self.sumb.adjointvars.ndesignaoa] *= numpy.pi/180.0

        return dIda

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

    def saveAdjointMatrix(self, fileName):
        """ Save the adjoint matrix to a binary petsc file for
        possible future external testing

        Parameters
        ----------
        fileName : str
            Filename to use.
        """
        if self.adjointSetup:
            self.sumb.saveadjointmatrix(fileName)
        else:
            raise Error('Cannot save matrix since adjoint not setup.')

    def saveAdjointPC(self, fileName):
        """ Save the adjoint preconditioning matrix to a binary petsc
        file for possible future external testing

        Parameters
        ----------
        fileName : str
            Filename to use.
        """

        if self.adjointSetup and self.getOption('approxpc'):
            self.sumb.saveadjointpc(fileName)
        else:
            raise Error('Cannot save PC matrix since adjoint not setup.')

    def saveAdjointRHS(self, fileName):
        """ Save the adjoint RHS to a binary petsc
        file for possible future external testing

        Parameters
        ----------
        fileName : str
            Filename to use.
        """
        if self.adjointSetup:
            self.sumb.saveadjointrhs(fileName)
        else:
            raise Error('Cannot save RHS since adjoint not setup.')

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
        initCFD Norm: Norm the adjoint starts with (zero adjoint)
        startCFD Norm: Norm at the start of adjoint call
        finalCFD Norm: Norm at the end of adjoint call
        '''
        startRes = self.sumb.adjointpetsc.adjreshist[0]
        finalIt  = self.sumb.adjointpetsc.adjconvits
        finalRes = self.sumb.adjointpetsc.adjreshist[finalIt-1]
        fail = self.sumb.killsignals.adjointfailed

        return startRes, finalRes, fail

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

    def getdRdXvTPsi(self, objective, groupName=None):
        """
        Compute the product of (dR/dXv)^T * psi for the objective
        given in 'objective'. If the mesh is present this will also
        pass it through the warping code and complete the following
        calculation:

        dXs = (dXv/dXs)^T * (dR/dXv)^T * psi

        Parameters
        ----------
        groupName : str
            Family name to use to section out just part of dXs
        objective : str
            Object to use. Must already be computed
            """

        # Get objective
        obj, aeroObj = self._getObjective(objective)

        if obj in self.curAP.sumbData.adjoints.keys():
            psi = self.curAP.sumbData.adjoints[obj]
        else:
            raise Error('%s adjoint for current aeroProblem is not computed.'%
                    obj)

        # Now call getdrdxvtpsi WITH the psi vector:
        dxvSolver = self.sumb.getdrdxvtpsi(self.getSpatialSize(), psi)

        # If we are doing a prescribed motion TS motion, we need to
        # convert this back to a single instance
        if self._prescribedTSMotion():
            ndof_1_instance = self.sumb.adjointvars.nnodeslocal[0]*3
            dxvSolver = self.sumb.spectralprecscribedmotion(
                dxvSolver, ndof_1_instance)

        if groupName is not None:
            self.mesh.warpDeriv(dxvSolver)
            dxs = self.mesh.getdXs(groupName)
            return dxs
        else:
            return dxvSolver

    def _prescribedTSMotion(self):
        """Determine if we have prescribed motion timespectral analysis"""

        if (self.getOption('alphamode') or self.getOption('betamode') or 
            self.getOption('machmode') or self.getOption('pmode') or 
            self.getOption('qmode') or self.getOption('rmode') or 
            self.getOption('altitudemode')):
            return True
        else:
            return False

    def getdRdXvTVec(self, inVec, groupName):
        """
        Compute the product of (dXv/dXs)^T * (dR/dXv)^T * inVec. It is
        assumed the mesh is present and groupName is defined.

        Parameters
        ----------
        inVec : array
            Array of size getStateAdjointSize()
        groupName : str
            Family name to use to section out just part of dXs
            """

        dxvSolver = self.sumb.getdrdxvtpsi(self.getSpatialSize(), inVec)
        self.mesh.warpDeriv(dxvSolver)
        dxs = self.mesh.getdXs(groupName)

        return dxs

    def getdRdaPsi(self,  psi):

        if self.nDVAero > 0:
            dIda = self.sumb.getdrdapsi(self.nDVAero, psi)
        else:
            dIda = numpy.zeros((0))

        return dIda

    def getdRdwTVec(self, inVec, outVec):
        """ Compute the result: outVec = dRdw^T * inVec

        Parameters
        ----------
        inVec : arrary
            Arrary of size getAdjointStateSize()
        outVec : array
            Ouput result. Same size as inVec.
        """
        outVec = self.sumb.getdrdwtvec(inVec, outVec)

        return outVec

    def getdFdxAero(self, iDV, groupName=None):
        """Potential aerodynamic variable dependence on forces. This
        is zero for all aerodynamic variables in SUmb"""
        return None

    def getdFdxVec(self, groupName, vec):
        # Calculate dFdx * vec and return the result
        vec = self.mesh.expandVectorByFamily(groupName, vec)
        if len(vec) > 0:
            vec = self.sumb.getdfdxvec(numpy.ravel(vec))
        vec = self.mesh.sectionVectorByFamily(groupName, vec)

        return vec

    def getdFdxTVec(self, groupName, vec):
        # Calculate dFdx^T * vec and return the result
        vec = self.mesh.expandVectorByFamily(groupName, vec)
        if len(vec) > 0:
            vec = self.sumb.getdfdxtvec(numpy.ravel(vec))
        vec = self.mesh.sectionVectorByFamily(groupName, vec)

        return vec

    def computeObjPartials(self, objective, forcePoints=None):

        obj, aeroObj = self._getObjective(objective)
        if aeroObj:
            objNum = self.sumbCostFunctions[obj]
            if self.getOption('useReverseModeAD'):
                # Note: Computeobjective partials MUST be called with the full
                # force pt list.
                if forcePoints is None:
                    [npts, ncells] = self.sumb.getforcesize()
                    nTS = self.sumb.inputtimespectral.ntimeintervalsspectral
                    forcePoints = numpy.zeros((nTS, npts, 3), self.dtype)
                    for i in xrange(nTS):
                        forcePoints[i] = self.getForcePoints(TS=i)

                self.sumb.computeobjpartials(
                    objNum, forcePoints.T, True, True)
            else:
                self.sumb.computeobjectivepartialsfwd(objNum)

            # Store the current RHS
            self.curAP.sumbData.adjointRHS = obj
        else:
            self.sumb.zeroobjpartials(True, True)

    def getdIdx(self, objective, forcePoints=None, TS=0, groupName=None):
        self._setupAdjoint()

        # Compute the partials
        self.computeObjPartials(objective, forcePoints)
        dXv = numpy.zeros(self.getSpatialSize())
        self.sumb.getdidx(dXv)

        # If we are doing a prescribed motion TS motion, we need to
        # convert this back to a single instance
        if self._prescribedTSMotion():
            ndof_1_instance = self.sumb.adjointvars.nnodeslocal[0]*3
            dXv = self.sumb.spectralprecscribedmotion(dXv, ndof_1_instance)
        # end if

        if groupName is not None:
            # We have a decision to make here: If we have euler
            # analysis, we can do a "surfOnly" meshDerivative since
            # there is no information on the interior anyway. However,
            # if we have a viscous analysis, then we DO have to do a
            # proper mesh warp, its fairly costly, but worth it.

            if self.getOption('equationType') == 'euler':
               self.mesh.warpDeriv(dXv, surfOnly=True)
            else:
                self.mesh.warpDeriv(dXv, surfOnly=False)
            # end if

            dxs = self.mesh.getdXs(groupName)
            return dxs
        else:
            return dXv

    def getdIda(self, objective, forcePoints=None):

        obj, aeroObj = self._getObjective(objective)

        if self.nDVAero > 0:
            self.computeObjPartials(objective, forcePoints)
            if aeroObj:
                dIdaLocal = self.sumb.adjointvars.dida
            else:
                dIdaLocal = numpy.zeros_like(self.sumb.adjointvars.dida)

            # We must MPI all reuduce
            dIda = self.comm.allreduce(dIdaLocal, op=MPI.SUM)
        else:
            dIda = numpy.zeros((0))

        return dIda

    def getdIdw(self, dIdw, objective, forcePoints=None):
        obj, aeroObj = self._getObjective(objective)
        if aeroObj:
            self.computeObjPartials(objective, forcePoints)
            dIdw = self.sumb.getdidw(dIdw)

        return dIdw

    def computeMatrixFreeProductFwd(self, xDVdot=None, wDot=None, residualDeriv=True, funcDeriv=True):
        
        if xDVdot is None and wDot is None:
            raise Error('computeMatrixFreeProductFwd: xDVdot and wDot cannot both be None')
        
        self._setAeroDVs()

        if xDVdot is not None:
            # Do the sptatial design variables -> Go through the geometry + mesh warp
            xsdot = self.DVGeo.totalSensitivityProd(xDVdot, self.curAP.ptSetName)
            xvdot = self.mesh.warpDerivFwd(xsdot)
            # Do the aerodynamic design variables
            if self.nDVAero > 0:
                extradot = numpy.zeros(self.nDVAero)
                for key in self.aeroDVs:
                    execStr = 'mapping = self.sumb.%s'%self.possibleAeroDVs[key.lower()]
                    exec(execStr)
                    if key == 'alpha':
                        # convert angle of attack to radians
                        extradot[mapping] = xDVdot[key]*(numpy.pi/180.0)
                    else:
                        extradot[mapping] = xDVdot[key]
            else:
                extradot = numpy.zeros(1)
            useSpatial = True
        else:
            xvdot = numpy.zeros(self.getSpatialSize())
            if self.nDVAero > 0:
                extradot = numpy.zeros(self.nDVAero)
            else:
                extradot = numpy.zeros(1)
            useSpatial = False

        if wDot is not None:
            useState = True
        else:
            wDot = numpy.zeros(self.getStateSize())
            useState = False

        costSize = self.sumb.costfunctions.ncostfunction
        dwdot, tmp = self.sumb.computematrixfreeproductfwd(
            xvdot, extradot, wDot, useSpatial, useState, costSize)

        funcsdot = {}
        for f in self.curAP.evalFuncs:
            mapping = self.sumbCostFunctions[self.possibleObjectives[f]]
            funcsdot[f] = tmp[mapping - 1]
            
        if residualDeriv and funcDeriv:
            return dwdot, funcsdot
        elif residualDeriv:
            return dwdot
        else:
            return funcsdot

    def computeMatrixFreeProductBwd(self, resBar=None, funcsBar=None, wDeriv=True, xDvDeriv=True):

        if resBar is None and funcsBar is None:
            raise Error("computeMatrixFreeProductBwd: resBar and funcsBar"
                        " cannot both be None")
        
        if resBar is None:
            resBar = numpy.zeros(self.getStateSize())

        if funcsBar is None:
            funcsBar = numpy.zeros(self.sumb.costfunctions.ncostfunction)
        else:
            tmp = numpy.zeros(self.sumb.costfunctions.ncostfunction)
            # Extract out the seeds
            for f in funcsBar:
                mapping = self.sumbCostFunctions[self.possibleObjectives[f.lower()]]
                tmp[mapping-1] = funcsBar[f]
            funcsBar = tmp

        useSpatial = False
        useState = False
        if wDeriv:
            useState = True
        if xDvDeriv:
            useSpatial = True
            
        if self.nDVAero > 0:
            aeroDVsize = self.nDVAero
        else:
            aeroDVsize = 1
        
        xvbar, extrabar, wbar = self.sumb.computematrixfreeproductbwd(
            resBar, funcsBar, useSpatial, useState, self.getSpatialSize(), aeroDVsize)
            
        if xDvDeriv:
            xdvbar = {}
            self.mesh.warpDeriv(xvbar)
            xsbar = self.mesh.getdXs('all')
            xdvbar.update(self.DVGeo.totalSensitivity(xsbar, 
                self.curAP.ptSetName, self.comm, config=self.curAP.name))
               
            # We also need to add in the aero derivatives here
            for key in self.aeroDVs:
                execStr = 'mapping = self.sumb.%s'%self.possibleAeroDVs[key]
                exec(execStr)
                if key == 'alpha':
                    # convert angle of attack to degrees
                    xdvbar[key] = extrabar[mapping]*(numpy.pi/180.0)
                else:
                    xdvbar[key] = extrabar[mapping]
            
        if wDeriv and xDvDeriv:
            return wbar, xdvbar
        elif wDeriv:
            return wbar
        else:
            return xdvbar

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
        self.sumb.setadjoint(adjoint)
        if objective is not None:
            obj, aeroObj = self._getObjective(objective)
            self.curAP.sumbData.adjoints[obj] = adjoint.copy()

    def getAdjoint(self, objective):
        """ Return the adjoint values for objective if they
        exist. Otherwise just return zeros"""
        obj, aeroObj = self._getObjective(objective)

        if obj in self.curAP.sumbData.adjoints:
            return self.curAP.sumbData.adjoints[obj]
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

    def computedSdwTVec(self, inVec, outVec, groupName):
        """This function computes: outVec = outVec + dFdw^T*inVec"""
        phi = self.mesh.expandVectorByFamily(groupName, inVec)
        outVec = self.sumb.getdfdwtvec(numpy.ravel(phi), outVec)

        return outVec

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

    def _getObjective(self, objective):
        """Check to see if objective is one of the possible
        objective. If it is, return the obj value for SUmb and
        True. Otherwise simply return the objective string and
        False"""
        if objective in self.possibleObjectives.keys():
            obj = self.possibleObjectives[objective]
            aeroObj = True
        else:
            obj = objective
            aeroObj = False

        return obj, aeroObj

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

        # Try to the option in the option dictionary
        defOptions = self.options['defaults']
        try:
            defOptions[name]
        except:
            if self.comm.rank == 0:
                SUMBWarning("Option '%-30s' is not a valid SUmb Option |"%name)
            return

        # Now we know the option exists, lets check if the type is ok:
        if type(value) == self.options[name][0]:
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
        # mane sure the name is lowercase

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
            'solutionprecision':[str,'single'],
            'gridprecision':[str,'double'],
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
            'turbresscale':[float, 10000.0],

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
            'usematrixfreedrdw':[bool, False],
            'usematrixfreedrdx':[bool, False],

            # ADjoint debugger
            'firstrun':[bool, True],
            'verifystate':[bool, True],
            'verifyspatial':[bool, True],
            'verifyextra':[bool, True],
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
            'turbresscale':{'location':'inputiteration.turbresscale'},
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
                                     'inplicitrk':self.sumb.inputunsteady.implicitrk,
                                     'md':self.sumb.inputunsteady.md,
                                     'location':'inputunsteady.timeintegrationscheme'},
            'timeaccuracy':{'location':'inputunsteady.timeaccuracy'},
            'ntimestepscoarse':{'location':'inputunsteady.ntimestepscoarse'},
            'ntimestepsfine':{'location':'inputunsteady.ntimestepsfine'},
            'deltat':{'location':'inputunsteady.deltat'},

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
            'usematrixfreedrdx':{'location':'inputadjoint.usematrixfreedrdx'},
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
            'liftindex'
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
                          ]

        return ignoreOptions, deprecatedOptions, specialOptions

    def _getObjectivesAndDVs(self):

        possibleObjectives = {
            'lift':'Lift',
            'drag':'Drag',
            'cl':'Cl', 'cd':'Cd',
            'fx':'Fx', 'fy':'Fy','fz':'Fz',
            'cfx':'cFx','cfy':'cFy', 'cfz':'cFz',
            'mx':'Mx', 'my':'My', 'mz':'Mz',
            'cmx':'cMx', 'cmy':'cMy', 'cmz':'cMz',
            'cm0':'cM0',
            'cmzalpha':'cMzAlpha',
            'cmzalphadot':'cMzAlphaDot',
            'cl0':'cl0',
            'clalpha':'clAlpha',
            'clalphadot':'clAlphaDot',
            'cd0':'cd0',
            'cdalpha':'cdAlpha',
            'cdalphadot':'cdAlphaDot',
            'cmzq':'cmzq',
            'cmzqdot':'cmzqdot',
            'clq':'clq',
            'clqdot':'clqDot',
            'cbend':'cBend',
            'sepsensor':'sepsensor',
            'cavitation':'cavitation',
            }

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
            'pointrefx':'adjointvars.ndesignpointrefx',
            'pointrefy':'adjointvars.ndesignpointrefy',
            'pointrefz':'adjointvars.ndesignpointrefz',
            'chordRef':'adjointvars.ndesignlengthref',
            'areaRef':'adjointvars.ndesignsurfaceref',
            'disserror':'adjointvars.ndesigndisserror'
            }

        # This is SUmb's internal mapping for cost functions
        sumbCostFunctions = {
            'Lift':self.sumb.costfunctions.costfunclift,
            'Drag':self.sumb.costfunctions.costfuncdrag,
            'Cl'  :self.sumb.costfunctions.costfuncliftcoef,
            'Cd'  :self.sumb.costfunctions.costfuncdragcoef,
            'Fx'  :self.sumb.costfunctions.costfuncforcex,
            'Fy'  :self.sumb.costfunctions.costfuncforcey,
            'Fz'  :self.sumb.costfunctions.costfuncforcez,
            'cFx' :self.sumb.costfunctions.costfuncforcexcoef,
            'cFy' :self.sumb.costfunctions.costfuncforceycoef,
            'cFz' :self.sumb.costfunctions.costfuncforcezcoef,
            'Mx'  :self.sumb.costfunctions.costfuncmomx,
            'My'  :self.sumb.costfunctions.costfuncmomy,
            'Mz'  :self.sumb.costfunctions.costfuncmomz,
            'cMx':self.sumb.costfunctions.costfuncmomxcoef,
            'cMy':self.sumb.costfunctions.costfuncmomycoef,
            'cMz':self.sumb.costfunctions.costfuncmomzcoef,
            'cM0':self.sumb.costfunctions.costfunccm0,
            'cMzAlpha':self.sumb.costfunctions.costfunccmzalpha,
            'cMzAlphaDot':self.sumb.costfunctions.costfunccmzalphadot,
            'cl0':self.sumb.costfunctions.costfunccl0,
            'clAlpha':self.sumb.costfunctions.costfuncclalpha,
            'clAlphaDot':self.sumb.costfunctions.costfuncclalphadot,
            'cfy0':self.sumb.costfunctions.costfunccfy0,
            'cfyAlpha':self.sumb.costfunctions.costfunccfyalpha,
            'cfyAlphaDot':self.sumb.costfunctions.costfunccfyalphadot,
            'cd0':self.sumb.costfunctions.costfunccd0,
            'cdAlpha':self.sumb.costfunctions.costfunccdalpha,
            'cdAlphaDot':self.sumb.costfunctions.costfunccdalphadot,
            'cmzq':self.sumb.costfunctions.costfunccmzq,
            'cmzqDot':self.sumb.costfunctions.costfunccmzqdot,
            'clq':self.sumb.costfunctions.costfuncclq,
            'clqDot':self.sumb.costfunctions.costfuncclqdot,
            'cBend':self.sumb.costfunctions.costfuncbendingcoef,
            'sepsensor':self.sumb.costfunctions.costfuncsepsensor,
            'cavitation':self.sumb.costfunctions.costfunccavitation,
            }

        return possibleObjectives, possibleAeroDVs, sumbCostFunctions

class sumbFlowCase(object):
    """
    Class containing the data that SUmb requires to be saved to an
    aeroProblem to permit the analysis of multiple flow cases
    """
    def __init__(self):
        self.stateInfo = None
        self.adjoints = {}
        self.coords = None
        self.callCounter = -1
        self.adjointRHS = None

