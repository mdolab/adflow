!
!      ******************************************************************
!      *                                                                *
!      * File:          iteration.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-13-2003                                      *
!      * Last modified: 09-20-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       module iteration
!
!      ******************************************************************
!      *                                                                *
!      * This module contains the iteration parameters mainly used in   *
!      * solver.                                                        *
!      *                                                                *
!      ******************************************************************
!
       use constants, only: intType, realType
       implicit none

       ! groundLevel:  Current ground level of the computation. Needed
       !               to determine what kind of action must be
       !               undertaken. E.G. On the coarse grids no solution
       !               will be written.
       ! currentLevel: MG level at which the compution currently resides.
       ! rkStage:      Current runge kutta stage. Needed to determine
       !               whether or not the artificial dissipation terms
       !               must be computed.

       integer(kind=intType) :: groundLevel, currentLevel
       integer(kind=intType) :: rkStage, Subit

       ! nStepsCycling: Number of steps in the current cycling strategy
       ! cycling:       The corresponding array defining the multigrid
       !                cycling strategy.

       integer(kind=intType) :: nStepsCycling
       integer(kind=intType), dimension(:), allocatable :: cycling

       ! iterTot: Total number of iterations on the current grid;
       !          a restart is not included in this count.

       integer(kind=intType) :: iterTot

       ! rFil : coefficient to control the fraction of the dissipation
       !        residual of the previous runge-kutta stage.

       real(kind=realType) :: rFil,rfilb

       ! t0Solver: Reference time for the solver.

       real(kind=realType) :: t0Solver

       ! converged:             Whether or not the solution has been
       !                        converged.
       ! exchangePressureEarly: Whether or not the pressure must be
       !                        exchanged early, i.e. before the
       !                        boundary conditions are applied.
       !                        This must be done for a correct treatment
       !                        of normal momentum boundary condition,
       !                        but it requires an extra call to the
       !                        halo routines.

       logical :: converged
       logical :: exchangePressureEarly

       ! standAloneMode:   Whether or not an executable in stand alone
       !                   mode is built.
       ! changing_Grid:    Whether or not the grid changes in time.
       !                   In stand alone mode this only happens when
       !                   moving parts are present. In a
       !                   multi-disciplinary environment more options
       !                   are possible, i.e. deforming meshes.
       ! deforming_Grid:   Whether or not the grid deforms; this can
       !                   only happen for a multi-disciplinary,
       !                   usually aero-elastic problem.
       ! changingOverset:  Whether or not the overset connectivity needs
       !                   to be updated at each time step, due to 
       !                   moving or deforming grids.
 
       logical :: standAloneMode, changing_Grid, deforming_Grid
       logical :: changingOverset

       ! nOldSolAvail:     Number of available old solutions for
       !                   the time integration.
       ! nOldLevels:       Number of old levels needed in the time
       !                   integration scheme.
       ! coefTime(0:nOld): The coefficients in the time integrator
       !                   for unsteady applications.

       integer(kind=intType) :: nOldSolAvail, nOldLevels
       real(kind=realType), dimension(:), allocatable :: coefTime

       ! iterType: The type of iteration performed. Will be one of RK,
       ! DADI, ANK or NK ( or None on the 0th evaluation)
       character(len=6) :: iterType

       ! approxTotalIts : A rough approximation of the total number of
       ! function evaluations. An RK or DADI multi grid iteration
       ! counts as 1.  ANK steps count as 1 + number of KSP
       ! iterations. NK steps count the total number of function
       ! evalautions either for mat-vecs or during a line search. It
       ! is this value that is checked again nCycles for doing too
       ! much work.
       integer(kind=intType) :: approxTotalIts
       
       ! Variable for monitoring the current CFL depending on the type
       ! of iteration
       real(kind=realTYpe) :: CFLMonitor

       ! *******************************
       ! Added by HDN
       ! *******************************
       ! nALEMeshes:                Number of ALE levels for intermediate mesh
       !                            between two steps
       ! nALEsteps:                 Number of ALE steps at one time step
       ! coefTimeALE(nALEsteps):    The weighting coefficients to average the fluxes
       ! coefMeshALE(nALEmeshes,2): The coefficients to interpolate the mesh
       integer(kind=intType)                            :: nALEMeshes, nALEsteps
       real(kind=realType), dimension(:),   allocatable :: coefTimeALE
       real(kind=realType), dimension(:,:), allocatable :: coefMeshALE

       ! timeSpectralGridsNotWritten: Whether or not grid files have
       !                              already been written in time
       !                              spectral mode. In this way
       !                              it is avoided that files are
       !                              written multiple times.
       ! oldSolWritten(nOldLevels-1): Logicals to indicate whether
       !                              or not old solution levels
       !                              have been written in
       !                              unsteady mode.

       logical :: timeSpectralGridsNotWritten

       logical, dimension(:), allocatable :: oldSolWritten

       external signalwritecallback

       end module iteration
