"""
@File    :   cpu_driver.py
@Date    :   07/22/2026
@Author  :   Safa Bakhshi
@Description : run CPU ADflow on test grid and save the restart file
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
    "CFL": 1.5,
    "CFLCoarse": 1.25,
    "MGCycle": MGCycle,
    "MGStartLevel": -1,
    "nCyclesCoarse": 250,
    "nCycles": 100000,
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

CFDSolver(ap)

# CFDSolver.writeSolution(outputDir="./INPUT/",baseName="test_mesh_nb1_nc32_restart")




