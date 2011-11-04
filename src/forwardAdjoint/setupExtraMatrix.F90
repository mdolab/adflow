subroutine setupExtraMatrix
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Assembles the dRda Matrix                                      *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use inputADjoint
  use communication
  implicit none
  !
  !     Local variables.

  logical :: useAD

  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
#ifndef USE_NO_PETSC

  useAD = .True.
  print *,'calling setupexraResidualMatrix'
  call setupExtraResidualMatrix(drda,useAD)

#endif
end subroutine setupExtraMatrix
