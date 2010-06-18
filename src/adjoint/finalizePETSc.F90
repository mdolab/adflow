!
!     ******************************************************************
!     *                                                                *
!     * File:          finalizePETSc.F90                               *
!     * Author:        Andre C. Marta                                  *
!     * Starting date: 08-24-2005                                      *
!     * Last modified: 03-01-2007                                      *
!     *                                                                *
!     ******************************************************************
!
      subroutine finalizePETSc
!
!     ******************************************************************
!     *                                                                *
!     * Finalize PETSc by calling the appropriate routine              *
!     * PetscFinalize provided in the PETSc library. This              *
!     * automatically calls MPI_Finalize().                            *
!     *                                                                *
!     ******************************************************************
!
      use ADjointPETSc
      implicit none
!
!     ******************************************************************
!     *                                                                *
!     * Begin execution.                                               *
!     *                                                                *
!     ******************************************************************
!
#ifndef USE_NO_PETSC

      ! This is not necessary since MPI_Finalize is called in SUmb.c

      ! PetscFinalize - Checks for options to be called at the
      ! conclusion of the prog. MPI_Finalize() is called only if the
      ! user had not called MPI_Init() before calling PetscInitialize().
      !
      ! Synopsis
      !
      ! #include "petsc.h"   
      ! call PetscFinalize(PetscErrorCode ierr)
      !
      ! Collective on SUMB_PETSC_COMM_WORLD
      !
      ! see .../petsc/docs/manualpages/Sys/PetscFinalize.html
      ! or PETSc users manual, pp.17

!      call PetscFinalize(PETScIerr)

!      if( PETScIerr/=0 ) &
!        call terminate("finalizePETSc", "Could not finalize PETSc")

      ! Release the memory for the convergence residual history.

      if( allocated(adjResHist) ) deallocate(adjResHist)

      ! Send some feedback to screen.

      if( PETScRank==0 ) &
        write(*,10) "finalizing PETSc..."

      ! Flush the output buffer and synchronize the processors.

      call f77flush()
      call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)

      ! Output format.

   10 format(a)

#endif

      end subroutine finalizePETSc
