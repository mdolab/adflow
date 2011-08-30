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

  call createStatePETScVars
  call createSpatialPETScVars
  call createExtraPETScVars
  call createCouplingPETScVars
  call createPETScKsp

#endif

end subroutine createPETScVars
