.. _sumb_performance:

Performance
===========

This section is intended to give a rough measure of the performance
that users may expect from SUmb. It is intended to give users some
indication if their simulation is performing as well as it should 
be. 

In general, simulations should be run on the *minimum* number of
processors possible. This will generally mean that memory usage will
be at, or near the maximum per core on the specific machine you are
using. Smaller processor counts minimize parallel losses, especially
from the breakdown of the effectiveness of the ASM/ILU preconditioner.

For the following analysis options:

* RANS
* SA turbulence model
* ILUFil = 2
* ASMOverlap = 2
* adjointSubspaceSize = 150
* viscPC = False
* useMatrixFreedRdw = True

SUmb will be able to run ~40\,000 cells/Gb of main memory. This
figure includes some additional overhead for the mesh movement
algorithm as well. Increasing the amount of ILU fill, the ASM overlap,
or subspace size will increase memory usage. The viscPC option will
increase the the total memory by a factor of approximately
2-3. Storing dRdw will increase memory usage by about 50%. Increasing
any of these options will increase the memory usage and reduce the
number of cells that will fit in 1Gb. 

Run-time performance can vary widely depending on the particular case
one is solving. A useful, but somewhat imperfect metric that can be
used to measure solution performance is: CPPH = # of cells converged/(proc *
hour). We will arbitrary say a solution is converged when the L2 norm
of the residual (totalR in the monitoring output) has dropped by 8
orders of magnitude. This is equivalent of setting L2convergence=1e-8
with the NKSolver. 

For simple cases with good meshes (think isolated wing with a pyHyp)
mesh and a modest number of cells (<1M), CPPH can exceed 1 million. As
a concrete example:

* 450k mesh
* 4 processors (Desktop machine)
* 400 sec solution time
* CPPH = 1 012 500

For much larger and more difficult case, such as complete
configuration: wing, body, nacelle, pylon, h-stab and v-stab, the CPPH
may be much lower ~200\,000. Another example:

* 5 200 000 cell mesh
* 64 processors 
* 1200 sec solution time
* CPPH = 244 000

These values are typical of using the NKSolver in transonic flow from
freestream condition with a reasonably efficient Runge-Kutta/DADI
start-up procedure before the NKSolver starts. For optimization,
subsequent solutions should take less time, often taking less than 1/2
of the time towards the end of an optimization.

In general, Euler cases will require somewhat less memory than the
RANS cases. However, with the matrix-free adjoint solver the fraction
reduction will not be much more than 5/6 of the RANS memory. Run time
for Euler cases will be typically be 2-5 times smaller than an
equivalently sized RANS case (number of cells). 
