!
!     ******************************************************************
!     *                                                                *
!     * File:          warpingInitializePETSc.F90                      *
!     * Author:        C.A.(Sandy) Mader                               *
!     * Starting date: 05-21-2009                                      *
!     * Last modified: 05-21-2009                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine warpingInitializePETSc
!
!     ******************************************************************
!     *                                                                *
!     * Initialize some PETSc values for use in the Warping routines.  *
!     *                                                                *
!     ******************************************************************
!
      use warpingPETSc
      use communication
      implicit none
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC


      ! Determine the communicator size and the processor rank.

      call MPI_Comm_size(PETSC_COMM_WORLD, PETScSize, PETScIerr)
      call MPI_Comm_rank(PETSC_COMM_WORLD, PETScRank, PETScIerr)

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "Initializing PETSc..."

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)

#endif

      end subroutine warpingInitializePETSc
