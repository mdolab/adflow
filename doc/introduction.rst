.. _adflow_introduction:

Introduction
============

``ADflow`` is a multi-block structured flow solver initially developed in 
the Stanford University under the sponsorship of the Department of
Energy Advanced Strategic Computing (ASC) Initiative . It solves the 
compressible Euler, laminar Navier-Stokes and Reynolds-Averaged Navier-Stokes 
equations. Although its primary objective in this program was to compute the
flows in the rotating components of jet engines, ADflow has been developed as 
a completely general solver and it is therefore applicable to a variety of 
other types of problems, including external aerodynamic flows.

``ADflow`` is a parallel code, suited for running on massively parallel platforms. 
The parallelization is hidden as much as possible from the end user, i.e. 
there is only one grid file, one volume solution file and one surface solution file. The
only thing the end user needs to do is to specify the number of processors 
he/she wants to use via the mpirun (or equivalent) command.

A summary of the various features that can be found in ``ADflow`` is given below:

* Compressible, URANS flow solver with various turbulence modeling options (Spalart-Allmaras, k-w, SST, v2-f)

* Multiblock structured approach with arbitrary connectivity. One-to-one mesh point matching with subfacing
  (C-0 multiblock) and point mismatched abutting meshes at block interfaces (C-1 multiblock) are allowed.
  CGNS I/O (mesh and solution) as well as native, MPI-IO parallel I/O option with back and forth conversion
  utilities.

* Massively parallel (both CPU and memory scalable) implementation using MPI.

* ALE Deforming grid implementation using ``pyWarp``

* Interface to conservative and consistent load and displacement transfer for aeroelastic computations.

* Multigrid, Runge-Kutta solver for the mean flow and DD-ADI solution methodology for the turbulence
  equations.

* Central difference discretization (second order in space) with various options for artificial dissipation.
  
* Adaptive wall functions for poor quality meshes.

* Unsteady time integration using second- or third-order (in time) backwards difference formulae (BDF) or a
  time-spectral approach for time-periodic flows.

* Fully parallel, scalable pre-processor responsible load balancing.
