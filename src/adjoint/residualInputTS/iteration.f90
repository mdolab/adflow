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
       use precision
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
       integer(kind=intType) :: rkStage

       ! nStepsCycling: Number of steps in the current cycling strategy
       ! cycling:       The corresponding array defining the multigrid
       !                cycling strategy.

       integer(kind=intType) :: nStepsCycling
       integer(kind=intType), dimension(:), allocatable :: cycling

       ! nMGVar: Number of variables to which the multigrid must be
       !         applied. For the Euler and laminar Navier-Stokes
       !         equations this is the number of flow variables; for
       !         RANS this is either the total number of independent
       !         variables (coupled solver) or the number of flow
       !         variables (segregated solver).
       ! nt1MG:  Starting index for the turbulent variables in MG.
       ! nt2MG:  Ending index for the turbulent variables in MG.
       !         For a segregated solver these values are such
       !         that nothing is done on the turbulent equations.

       integer(kind=intType) :: nMGVar, nt1MG, nt2MG

       ! restrictEddyVis: Whether or not the eddy viscosity must
       !                  be restricted to the coarser levels.
       ! turbSegregated:  Whether or not the turbulent equations
       !                  are solved segregatedly from the mean
       !                  flow equations.
       ! turbCoupled:     Whether or not the turbulent equations are
       !                  solved in a coupled manner with the mean
       !                  flow equations. The reason why both
       !                  turbCoupled and turbSegregated are used is
       !                  that everything must work for Euler and
       !                  laminar NS as well.

       logical :: restrictEddyVis, turbSegregated, turbCoupled

       ! iterTot: Total number of iterations on the current grid;
       !          a restart is not included in this count.

       integer(kind=intType) :: iterTot

       ! rFil : coefficient to control the fraction of the dissipation
       !        residual of the previous runge-kutta stage.

       real(kind=realType) :: rFil

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
       ! PV3Initialized:   Whether or not PV3 has been initialized,
       !                   for use in multidisciplinary problems where
       !                   solver is called multiple times
 
       logical :: standAloneMode, changing_Grid, deforming_Grid
       logical :: changingOverset, PV3Initialized = .false.

       ! nOldSolAvail:     Number of available old solutions for
       !                   the time integration.
       ! nOldLevels:       Number of old levels needed in the time
       !                   integration scheme.
       ! coefTime(0:nOld): The coefficients in the time integrator
       !                   for unsteady applications.

       integer(kind=intType) :: nOldSolAvail, nOldLevels
       real(kind=realType), dimension(:), allocatable :: coefTime

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

       end module iteration
