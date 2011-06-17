
subroutine getdFdxVec(ndof,vec_in,vec_out)
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the current adjoint vector by dRdXv to get a vector   *
  !     * of length Xv. This is collective communication part            *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc, only: dFdx,fVec1,fVec2
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in)  :: vec_in(ndof)
  real(kind=realType), intent(out)  :: vec_out(ndof)
  real(kind=realType) :: val,val2
  integer(kind=intType) :: ierr,rowstart,rowend,i
  vec_out(:) = 0.0

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_in,fVec1,ierr)
  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_out,fVec2,ierr)

  call MatMult(dFdx,fVec1,fVec2,ierr)
  call VecNorm(fVec2,NORM_2,val,ierr)

  call VecDestroy(fVec1,ierr)
  call VecDestroy(fVec2,ierr)


end subroutine getdFdxVec

