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

  use ADjointPETSc, only : petsc_comm_world
  use communication, only : sumb_comm_world
  implicit none

  PETSC_COMM_WORLD= SUMB_COMM_WORLD
  call initPETScWrap()

end subroutine initializePETSc
