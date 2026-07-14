import argparse
from mpi4py import MPI
from baseclasses import AeroProblem
from adflow import ADFLOW
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker

updatebcs = True
nproc = 1
ncell = 32
# gridFile = f"./INPUT/test_mesh_{ncell}_{nproc}.cgns"
gridFile = "INPUT/test_mesh_nb1_nc32_perturbed_restart.cgns"
# gridFile = "test_mesh_nb1_nc2_restart.cgns"
n = 32
outputDirectory = "./"

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
    "restartFile": gridFile,
    # "gridFile": "test_mesh_nb1_nc2_restart.cgns",
    # "restartFile": "INPUT/test_mesh_nb1_nc32_restart.cgns",
    "outputDirectory": outputDirectory,
    "solutionPrecision": "double",
    # "equationType": "Euler",
    "equationType": "laminar NS",
    # "equationType": "RANS",
    "CFL": 1.5,
    "CFLCoarse": 1.25,
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
}


ap = AeroProblem(
    name=name,
    alpha=alpha,
    beta=45,
    mach=mach,
    altitude=altitude,
    areaRef=areaRef,
    chordRef=chordRef,
    evalFuncs=["cl", "cd"],
)

# Create solver
CFDSolver = ADFLOW(options=aeroOptions, comm=MPI.COMM_WORLD)
CFDSolver.setAeroProblem(ap)
# CFDSolver(ap)
# CFDSolver.adflow.initializeflow.setuniformflow()
resCPU1 = CFDSolver.getResidual(ap)
# resCPU1 *= 0.0
# resCPU1 = CFDSolver.adflow.nksolver.getres(resCPU1)
# states = CFDSolver.getStates()
# ndimw = states.shape[0]
# np.random.seed(10)
# states  += np.random.rand(ndimw)
# CFDSolver.setStates(states)
# resCPU2 = CFDSolver.getResidual(ap)

ndimw = len(resCPU1)
nw = CFDSolver.adflow.flowvarrefstate.nw
nIters = 10

# One-time GPU setup + host->device copy, paid ONCE outside the timed loop below.
CFDSolver.adflow.cudablock.allocatecudablock()
CFDSolver.adflow.cudablock.copycudaconstantstodevice()
CFDSolver.adflow.cudablock.copycudablocktodevice()

# Timed region: the core GPU residual only (calculateCudaResidual). No host<->device
# copies happen inside this loop, so profiling it measures kernel time, not transfers.
for _ in range(nIters):
    CFDSolver.adflow.cudaresidual.calculatecudaresidual(True, 1, nw)

# Read the residual back to the host once and check correctness.
resGPU = CFDSolver.adflow.cudaresidual.getcudaresidual(ndimw)
print(len(resGPU))

# print("===========================================")
# i = np.arange(n)
# j = np.arange(n)
# I,J = np.meshgrid(i,j)
# for i in range(1,n):
#     x = n*n
#     err = np.abs(resGPU[x*(i-1)*5:x*i*5]-resCPU1[x*(i-1)*5:x*i*5])
#     err = err[::5] + err[1::5] + err[2::5] + err[3::5] + err[4::5]
#     err = err.reshape(n,n) + 1e-16
#     # print(i,np.linalg.norm(err,ord=np.inf),np.linalg.norm(np.abs(resGPU[x*(i-1)*5:x*i*5]-resCPU1[x*(i-1)*5:x*i*5]),ord=np.inf))
#     if np.linalg.norm(err,ord=np.inf) > 1e-12:
#         # if i == 1:
#             # print(i,err[0,:])
#             # print(i,err[1,:])
#             # print(i,err[2,:])
#             # print(i,err[3,:])
#         fig, ax = plt.subplots()
#         cs = ax.contourf(I,J,err,locator=ticker.LogLocator())
#         fig.colorbar(cs)
#         plt.savefig("k%d.pdf"%i)
#         plt.clf()
#     # print(np.linalg.norm(resGPU[x*(i-1)*5:x*i*5]-resCPU1[x*(i-1)*5:x*i*5],ord=np.inf))
#     # print(np.linalg.norm(resGPU[x*(i-1)*5:x*i*5]-resCPU1[x*(i-1)*5:x*i*5],ord=1))

# print("===========================================")

# print(resGPU[-20:]-resCPU1[-20:])
# for i in range(1,32):

#     print(i,np.linalg.norm(resCPU1[32*32*(i-1)*5:32*32*i*5]-resGPU[32*32*(i-1)*5:32*32*i*5]))
#     print(32*32*i*5)
# print(resGPU[:8*5])

# print(resCPU1[:30])
print("cpu norm:", np.linalg.norm(resCPU1, ord=1))
print("gpu norm:", np.linalg.norm(resGPU, ord=1))
print("error between cpu gpu infinity norm:", np.linalg.norm(resCPU1 - resGPU, ord=np.inf))
resGPU = np.reshape(resGPU, (32 * 32 * 32, 5))
resGPU = np.sum(resGPU, axis=1)
resGPU = np.reshape(resGPU, (32, 32, 32))
resCPU1 = np.reshape(resCPU1, (32 * 32 * 32, 5))
resCPU1 = np.sum(resCPU1, axis=1)
resCPU1 = np.reshape(resCPU1, (32, 32, 32))

i = np.arange(32)
j = np.arange(32)
I, J = np.meshgrid(i, j)
for k in range(32):
    plt.clf()
    plt.contourf(I, J, np.log10(np.abs(resGPU[:, :, k] - resCPU1[:, :, k])))
    plt.colorbar()
    plt.savefig("resdiff_%d.pdf" % k)
