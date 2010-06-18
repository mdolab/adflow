!
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScVars.F90                             *
!     * Author:        Andre C. Marta., C.A.(Sandy) Mader              *
!     * Starting date: 12-15-2005                                      *
!     * Last modified: 02-01-2008                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine createPETScVars
!
!     ******************************************************************
!     *                                                                *
!     * Create all the necessary PETSc objects to solve the discrete   *
!     * ADjoint equations, i.e., matrices, vectors and ksp.            *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      implicit none
!
!     Local variables.
!
      real(kind=realType), dimension(2) :: time
      real(kind=realType)               :: timeAdjLocal, timeAdj
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Creating PETSc objects..."

      ! Get the initial time.

      call cpu_time(time(1))

      ! Create the global adjoint matrix dR/dW and
      ! the auxiliar gradient matrices dR/da (extra design vars),
      ! dRdSigma (electrical conductivity) and 
      ! dRdx, dRdy, dRdz (spatial derivatives)

      call createPETScMat

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        print "(a)", "# ... matrices created;"

      ! Create the global right-hand side vector dJ/dW,
      ! the adjoint vector solution psi,
      ! the partial gradients dJ/da, dJ/dSigma and dJ/dx,dJ/dy,dJ/dz,
      ! and the total gradients dI/da, dI/dSigma and dI/dx,dI/dy,dI/dz.

      call createPETScVec

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        print "(a)", "# ... vectors created;"

      ! Create the global krylov object,
   
      call createPETScKsp

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        print "(a)", "# ... ksp created;"


      ! Get new time and compute the elapsed time.

      call cpu_time(time(2))
      timeAdjLocal = time(2)-time(1)

      ! Determine maximum time using MPI reduce
      ! with operation mpi_max.

      call mpi_reduce(timeAdjLocal, timeAdj, 1, sumb_real, &
                      mpi_max, 0, SUMB_PETSC_COMM_WORLD, PETScIerr)

      if( PETScRank==0 ) &
        write(*,20) "Creating PETSc objects time (s) = ", timeAdj

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output formats.

   10 format(a)
   20 format(a,1x,f8.2)

#endif

      end subroutine createPETScVars
