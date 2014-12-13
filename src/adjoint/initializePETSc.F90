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

  ! Call the C-version of the petsc initialize routine

  use ADjointPETSc
  use inputADjoint
  use communication
  implicit none

  PETSC_COMM_WORLD= SUMB_COMM_WORLD
  call initPETScWrap()
  !call PetscInitialize(PETSC_NULL_CHARACTER, PETScIErr)
  !call ECHK(PETScIerr, __FILE__, __LINE__)
end subroutine initializePETSc
