#! /usr/bin/env python

# File: dummy_mpi.py
# By: Patrick LeGresley
# Last Modified: 3/23/2004

# This is a 'dummy' MPI module that works in a manner similar to the dummy mpi
# calls in TFLO.  This is in no way complete but provides enough functionality
# for the current Python wrapping of TFLO.

__version__ = "0.1"

class communicator:

    def comm_create(self,arg):
        return communicator()

    def __getitem__(self,arg):
        return 1

    def __init__(self):
        self.rank = 0
        self.size = 1

    def __int__(self):
        return 1

    def bcast(self,value):
        return value
 
    def allreduce(self,value,operation):
        return value

    def barrier(self):
        pass

rank = 0
size = 1
COMM_WORLD = communicator()
SUM = None

