#!/usr/bin/python
"""
Overset Check: A minimal interface for checking overset connectivity

Copyright (c) 2016 by Dr. G. K. W. Kenway
All rights reserved. Not to be used for commercial purposes.

Developers:
-----------
- Dr. Gaetan K.W. Kenway (GKK)
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
from .pyADflow import ADFLOW

class Error(Exception):
    """
    Format the error message in a box to make it clear this
    was a expliclty raised exception.
    """
    def __init__(self, message):
        msg = '\n+'+'-'*78+'+'+'\n' + '| oversetCheck Error: '
        i = 21
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

# =============================================================================
# ADFLOW Class
# =============================================================================
class OversetCheck(ADFLOW):
    """
    Create the OversetCheck object.

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
        """
    def __init__(self, comm=None, options=None, debug=False):

        # Load the compiled module using MExt, allowing multiple
        # imports
        curDir = os.path.dirname(os.path.realpath(__file__))
        self.adflow = MExt.MExt('libadflow', [curDir], debug=debug)._module

        # Information for base class:
        name = 'ADFLOW'
        category = 'Three Dimensional CFD'
        informs = {}

        # If 'options' is not None, go through and make sure all keys
        # are lower case:
        if options is not None:
            for key in options.keys():
                options[key.lower()] = options.pop(key)
        else:
            raise Error("The 'options' keyword argument must be passed "
                        "adflow. The options dictionary must contain (at least) "
                        "the gridFile entry for the grid")

        # Load all the option/objective/DV information:
        defOpts = self._getDefOptions()
        self.optionMap = self._getOptionMap()
        self.ignoreOptions, self.deprecatedOptions, self.specialOptions = \
                           self._getSpecialOptionLists()

        self.possibleAeroDVs, self.basicCostFunctions = (
            self._getObjectivesAndDVs())

        # This is the real solver so dtype is 'd'
        self.dtype = 'd'

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

        # Initialize the inherited aerosolver
        AeroSolver.__init__(self, name, category, defOpts, informs,
                            options=options)

        # Initialize petec in case the user has not already
        self.adflow.initializepetsc()

        # Set the stand-alone adflow flag to false...this changes how
        # terminate calls are handled.
        self.adflow.iteration.standalonemode = False

        # Set the frompython flag to true... this also changes how
        # terminate calls are handled
        self.adflow.killsignals.frompython = True

        # Set default values
        self.adflow.setdefaultvalues()
        self.adflow.inputio.autoparameterupdate = False
        
        # Make sure all the params are ok
        for option in self.options:
            if option != 'defaults':
                self.setOption(option.lower(), self.options[option][1])

        dummyAP = AeroProblem(name='dummy', mach=0.5, altitude=10000.0,
                              areaRef=1.0, chordRef=1.0, alpha=0.0, degreePol=0,
                              coefPol=[0, 0], degreeFourier=1,
                              omegaFourier=6.28, sinCoefFourier=[0, 0],
                              cosCoefFourier=[0, 0])

        self.curAP = dummyAP
        self._setAeroProblemData(firstCall=True)

        # Finally complete loading
        self.adflow.dummyreadparamfile()
        self.adflow.partitionandreadgrid(False)
        self.adflow.preprocessingcustomoverset()
        # self.adflow.initflow()
        # self.adflow.preprocessingpart2()
        # self.adflow.initflowpart2()
        # self.adflow.preprocessingadjoint()
