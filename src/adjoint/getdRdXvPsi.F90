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
subroutine getdRdXvPsi(ndof,dXv)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the current adjoint vector by dRdXv to get a vector   *
  !     * of length Xv. This is collective communication part            *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc, only: dRdX,psi,gridVec
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(out)  :: dXv(ndof)

  integer(kind=intType) :: ierr

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,dXv,gridVec,ierr)
  print *,'myid,0:',myid,ierr
  call MatMultTranspose(dRdX,psi,GridVec,ierr)
  print *,'myid,1:',myid,ierr
  call VecDestroy(GridVec,ierr)
  print *,'myid,2:',myid,ierr

end subroutine getdRdXvPsi


