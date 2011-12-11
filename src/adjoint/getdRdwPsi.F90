!
!     ******************************************************************
!     *                                                                *
!     * File:          getdRdXvPsi.F90                                 *
!     * Author:        Gaetan Kenway                                   *
!     * Starting date: 10-08-2010                                      *
!     * Last modified: 10-08-2010                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine getdRdwTPsi(ndof,dRdwPsi)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the current adjoint vector by dRdw^T                  *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(out)  :: dRdwPsi(ndof)

  integer(kind=intType) :: ierr

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,dRdwPsi,wVec,ierr)
  call MatMult(dRdwT,psi,dRdwPsi,ierr)
  call VecDestroy(wVec,ierr)


end subroutine getdRdwTPsi

