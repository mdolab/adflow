!
!     ******************************************************************
!     *                                                                *
!     * File:          agumentRHS.f90                                  *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 10-08-2010                                      *
!     * Last modified: 10-08-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine agumentRHS(ndof,phi)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the current adjoint vector by dRdXv to get a vector   *
  !     * of length Xv. This is collective communication part            *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc 
  use ADjointVars 
  use communication
  use inputADjoint

  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in) :: phi(ndof)

  integer(kind=intType) :: ierr,ndims
  real(kind=realType) :: val
  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DETERMINE,phi,phic,ierr)
  !call EChk(ierr,__FILE__,__LINE__)
  ! Dump the result into adjointRHS
  call MatMultTranspose(dFdw,phic,adjointRHS,ierr)
  !call EChk(ierr,__FILE__,__LINE__)

  call vecDestroy(phic,ierr)
  !call EChk(ierr,__FILE__,__LINE__)
end subroutine agumentRHS


