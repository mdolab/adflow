subroutine initializePETSc

  ! Call the C-version of the petsc initialize routine

  use ADjointPETSc, only : petsc_comm_world
  use communication, only : sumb_comm_world
  implicit none

  PETSC_COMM_WORLD= SUMB_COMM_WORLD
  call initPETScWrap()

end subroutine initializePETSc
