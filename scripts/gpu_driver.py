"""
@File    :   cpu_driver.py
@Date    :   07/22/2026
@Author  :   Safa Bakhshi
@Description : run GPU ADflow on test grid and save the restart file
"""

# ==============================================================================
# Standard Python modules
# ==============================================================================
import os
import argparse
import time


# ==============================================================================
# External Python modules
# ==============================================================================
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker



# ==============================================================================
# Extension modules
# ==============================================================================


gridFile = "INPUT/test_mesh_nb1_nc32.cgns"



alpha = 0
mach = 0.5
bRef = 12.53
areaRef = 37.1982321487
chordRef = 2.96
MGCycle = "sg"
MGSTART = -1
altitude = 10000
name = "fc"



aeroOptions = {
    # Common Parameters
    "gridFile": gridFile,
    # "restartFile": gridFile,
    # "gridFile": "test_mesh_nb1_nc2_restart.cgns",
    # "restartFile": "INPUT/test_mesh_nb1_nc32_restart.cgns",
    "outputDirectory": "./INPUT/",
    "solutionPrecision": "double",
    "equationType": "Euler",
    # "equationType": "laminar NS",
    # "equationType": "RANS",
    "CFL": 0.1,
    "CFLCoarse": 0.1,
    "MGCycle": MGCycle,
    "MGStartLevel": -1,
    "nCyclesCoarse": 250,
    "nCycles": 10,
    "nsubiterturb": 5,
    "printIntro": False,
    "printAllOptions": False,
    # "vis2": 0.0,
    # "vis4": 0.0,
    "useblockettes": False,
    # "liftIndex": 3
    # "writeSolution": False,
    # "usetestwithbcs":updatebcs,
    # Debug viscous flux
    "useQCR":False,
    # "eulerWallTreatment":"linear pressure extrapolation",
    # Write Restart File
    "writeSurfaceSolution" : False,
    "writeVolumeSolution" : True,
    "useANKSolver":False,
    "useNKSolver":False,
    "smoother":"Runge-Kutta",
}

ap = AeroProblem(
    name=name,
    alpha=alpha,
    beta=0.0,
    mach=mach,
    altitude=altitude,
    areaRef=areaRef,
    chordRef=chordRef,
    evalFuncs=["cl", "cd"],
)

# Create solver
CFDSolver = ADFLOW(options=aeroOptions, debug=True, comm=MPI.COMM_WORLD)
resCPU1 = CFDSolver.getResidual(ap)


# Prepare CUDA API
ndimw = len(resCPU1)
nw = CFDSolver.adflow.flowvarrefstate.nw
nIters = 4

# One-time GPU setup + host->device copy, paid ONCE outside the timed loop below.
CFDSolver.adflow.cudaapi.setupcudaapi()

# Timed region: the core GPU residual only (calculateCudaResidual). No host<->device
# copies happen inside this loop, so profiling it measures kernel time, not transfers.
CFDSolver.adflow.cudaresidual.calculatecudaresidual(True, 1, nw)
print(f"Iteration {0} - L2 Total Res {CFDSolver.adflow.cudaresidual.getdeviceressum():.5e}")
for i in range(nIters):
    #
    CFDSolver.adflow.cudaresidual.calculatecudaresidual(True, 1, nw)
    CFDSolver.adflow.cudasmoothers.rungekuttasmoother()
    monres = CFDSolver.adflow.cudaresidual.getdeviceressum()
    print(f"Iteration {i+1} - L2 Total Res {monres:.5e}")


# Copy back to host
CFDSolver.adflow.cudablock.copycudablocktohost()


# Debug
CFDSolver.adflow.cudadebug.debugcomparecudadata()

# CFDSolver.adflow.cudablock.copycudastatestohost()
# CFDSolver.adflow.cudablock.copycudaresidualtohost()


# Destroy CUDA Memory
# CFDSolver.adflow.cudaapi.destroycudaapi()


# resCPU2 = CFDSolver.getResidual(ap)


# CFDSolver.writeSolution(outputDir="./INPUT/",baseName="test_mesh_nb1_nc32_GPU_sol")




