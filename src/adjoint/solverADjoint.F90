!
!     ******************************************************************
!     *                                                                *
!     * File:          solverADjoint.f90                               *
!     * Authors:       Juan J. Alonso, Joaquim R. R. A. Martins,       *
!     *                C.A.(Sandy) Mader                               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-17-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine solverADjoint
!
!     ******************************************************************
!     *                                                                *
!     * This subroutine is the main driver for the discrete adjoint    *
!     *  solver for SUmb.  The approach to the setup and               *
!     *  solution of the adjoint problem is based on the use of        *
!     *  automatic differentiation techniques (Tapenade from INRIA)    *
!     *  applied to a version of the residual computation for a single *
!     *  cell in the SUmb mesh.  This residual computation is meant to *
!     *  be in all ways identical to that performed by the actual flow *
!     *  solver (no approximations), including all flux routines,      *
!     *  limiters (if used) and all boundary conditions.               *
!     *                                                                *
!     * This version of the residual computation is differentiated (in *
!     *  reverse mode) using Tapenade and is used to compute the exact *
!     *  matrix for the discrete adjoint problem.  No approximations   *
!     *  are made.  In principle, this approach is able to handle      *
!     *  flow solver discretizations of arbitrary complexity in an     *
!     *  automatic fashion and can therefore dramatically reduce the   *
!     *  time of development of adjoint solvers while ensuring that    *
!     *  the sensitivity information obtained is exactly consistent    *
!     *  with (infinitely accurate) finite differencing of the solver. *
!     *                                                                *
!     * For efficiency reasons, this version of the adjoint solver     *
!     *  stores the adjoint matrix.  This is a significant overhead    *
!     *  (approximately an order of magnitude) in memory with respect  *
!     *  to the original solver.  We have decided to put up with this  *
!     *  for the time being as the approach has the potential to be    *
!     *  used, exactly, for governing equations as complex as entire   *
!     *  multiple equation turbulence models and even complex MHD      *
!     *  formulations.                                                 *
!     *                                                                *
!     * The resulting linear system is solved using Krylov subspace    *
!     *  methods embedded in the PETSc toolkit from Argonne National   *
!     *  laboratories.                                                 *
!     *                                                                *
!     * This is an implementation on the cell-centered version of      *
!     *  the SUmb solver.  It is not intended to be comprehensive.     *
!     *  The intent is to support the entire SUmb vertex solver once   *
!     *  it is completed and operational.                              *
!     *                                                                *
!     * This driver assumes that the flow solution has already been    *
!     *  computed and is stored in the block module. It also assumes   *
!     *  that PETSc has already been initialized.                      *
!     *                                                                *
!     * The adjoint solver is applied only to the finest grid level    *
!     *  and first time instance.                                      *
!     *                                                                *
!     ******************************************************************
!
      use ADjointVars
      use ADjointPETSc    !petsc_comm_world,petsc_rank
      use communication   ! myID, nProc
      use cgnsGrid        ! cgnsNFamilies
      use iteration       ! groundLevel
      use inputTimeSpectral !nIntervalsTimeSpectral
      use inputTSStabDeriv !TSstability
      implicit none
!
!     Local variables.
!
      integer(kind=intType) :: level

      integer(kind=intType) :: npts,nTS,obj_num
      real(kind=realType), allocatable,dimension(:,:,:) :: points

      
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
     if(myID == 0) then
	print *,'In SolverADjoint...'
     endif
#ifdef USE_NO_PETSC

      ! Sources compiled without PETSc support. Since PETSc is required
      ! to solve the ADjoint system of equations, the adjoint solver
      ! cannot be run.

      call terminate("solverADjoint", &
                     "Routine should not be called if no PETSc support &
                     &is selected.")

#else
      ! Set the relevant grid level.
      level = 1
      groundLevel = level

!
!     ******************************************************************
!     *                                                                *
!     * Preprocessing tasks.                                           *
!     *                                                                *
!     ******************************************************************
!
      call InitializePETSc

      !setup the extra design variabled for dRda
      call setupExtraDVsFortran
      
      ! Assertion testing, memory allocation and global number indexing.
      !print *,'preprocessing...'
      call preprocessingADjoint(level)

      !destroy the NK solver if it exists so that we don't have memory conflicts.
      !print *,'destroying...'
      call destroyNKSolver()

      !print *,'creating'
      !create all of the various matricies
      call createPETScVars()
      
      !print *,'setting up residual derivatives..'
      call setupAllResidualMatrices


      ! Reordered for ASM preconditioner
      ! Create the Krylov subspace linear solver context,
      ! the preconditioner context, and set their various options.

      call setupPETScKsp(level)
      
      !Get the surface points we will use to compute the objective
      call getForceSize(npts,nTS)
      
      allocate(points(3,npts,nTS))
      call getForcePoints(points,npts,nTS)

      obj_num = 1
      !compute the objective partial derivative
      call computeobjpartials(obj_num,points,npts,nTS)

      !Solve the system
      call solveAdjointTransposePETSc()

      call destroyPETScVars()

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      deallocate(points)
      if (allocated(dIda)) then
         deallocate(dida)
      end if
      if( PETScRank==0 ) &
        print "(a)", "# ... Adjoint solution complete;"

#endif

      end subroutine solverADjoint
