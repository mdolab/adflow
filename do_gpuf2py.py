from adflow import ADFLOW
import numpy as np
from mpi4py import MPI
solver = ADFLOW(options={
    "gridFile":"input_files/cube_4x4x4.cgns"
    })
# np.random.seed(MPI.COMM_WORLD.rank)
n = 1024
a = np.random.rand(n,n)
b = np.random.rand(n,n)
c = solver.adflow.matmult.my_matmult_api(a,b,n)
print(f"error on rank %d: %.16e"%(0,np.amax(c-a@b)))