!
!     ******************************************************************
!     *                                                                *
!     * File:          initializePETSc.F90                             *
!     * Author:        C.A.(Sandy) Mader, Andre C. Marta               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 01-14-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine initializePETSc
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Initialize PETSc by calling the appropriate routine provided   *
  !     * provided in the PETSc library - PetscInitialize. This          *
  !     * automatically calls MPI_Init() if MPI has not been previously  *
  !     * initialized. It also initializes some auxiliar variables       *
  !     * declared in module ADjointPETSc.                               *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use inputADjoint
  use communication
  implicit none

#ifndef USE_NO_PETSC

  call initPETScWrap()

  SUMB_PETSC_COMM_WORLD = SUmb_comm_world
  PETSC_COMM_SELF = SUmb_comm_self

  ! Allocate memory for the convergence residual history.
  if (.not. allocated(adjResHist))then
     allocate(adjResHist(adjMaxIter))
  endif

  ! Flush the output buffer and synchronize the processors.
  call f77flush()
  call mpi_barrier(SUMB_PETSC_COMM_WORLD, PETScIerr)
#endif
end subroutine initializePETSc
