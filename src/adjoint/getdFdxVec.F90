
subroutine getdFdxVec(ndof,vec_in,vec_out)
#ifndef USE_NO_PETSC

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply vec_in by dFdx to produce vec_out                     *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc, only: dFdx,fVec1,fVec2, PETSC_VIEWER_STDOUT_SELF
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in)  :: vec_in(ndof)
  real(kind=realType), intent(out)  :: vec_out(ndof)
  integer(kind=intType) :: ierr
  vec_out(:) = 0.0

  call VecPlaceArray(fVec1,vec_in,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(fVec2,vec_out,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatMult(dFdx,fVec1,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine getdFdxVec

subroutine getdFdxTVec(ndof,vec_in,vec_out)
#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply vec_in by dFdx to produce vec_out                     *
  !     *                                                                *
  !     ******************************************************************
  !
  use communication
  use ADjointPETSc, only: dFdx,fVec1,fVec2
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in)  :: vec_in(ndof)
  real(kind=realType), intent(out)  :: vec_out(ndof)
  integer(kind=intType) :: ierr
  vec_out(:) = 0.0
  
  call VecPlaceArray(fVec1,vec_in,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(fVec2,vec_out,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call MatMultTranspose(dFdx,fVec1,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine getdFdxTVec
