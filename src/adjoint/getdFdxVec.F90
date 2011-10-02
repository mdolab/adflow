
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
  use ADjointPETSc, only: dFdx,fVec1,fVec2
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in)  :: vec_in(ndof)
  real(kind=realType), intent(out)  :: vec_out(ndof)
  integer(kind=intType) :: ierr
  vec_out(:) = 0.0

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_in,fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_out,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)
     
  call MatMult(dFdx,fVec1,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine getdFdxVec

subroutine getdFdxVec_NULL
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
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType) ndof
  real(kind=realType) :: vec_in(0)
  real(kind=realType) :: vec_out(0)
  integer(kind=intType) :: ierr

  ndof = 0
  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_in,fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_out,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)
     
  call MatMult(dFdx,fVec1,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine getdFdxVec_NULL

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
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in)  :: vec_in(ndof)
  real(kind=realType), intent(out)  :: vec_out(ndof)
  integer(kind=intType) :: ierr
  vec_out(:) = 0.0

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_in,fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_out,fVec2,ierr)
  call EChk(PETScIerr,__FILE__,__LINE__)
     
  call MatMultTranspose(dFdx,fVec1,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine getdFdxTVec


subroutine getdFdxTVec_NULL
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
  use warpingPETSC 
  use inputADjoint
  implicit none

  integer(kind=intType) ndof
  real(kind=realType) :: vec_in(0)
  real(kind=realType) :: vec_out(0)
  integer(kind=intType) :: ierr

  ndof = 0
  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_in,fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(SUMB_COMM_WORLD,ndof,PETSC_DECIDE,vec_out,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)
     
  call MatMultTranspose(dFdx,fVec1,fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec1,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(fVec2,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine getdFdxTVec_NULL
