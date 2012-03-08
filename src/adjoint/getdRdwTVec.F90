subroutine getdRdwTVec(in_vec, out_vec, ndof)
#ifndef USE_NO_PETSC

  use communication
  use ADjointPETSc
  implicit none

  ! Input/Output
  integer(kind=intType), intent(in) :: ndof
  real(kind=realType), intent(in) :: in_vec(ndof)
  real(kind=realType), intent(inout) :: out_vec

  ! Working Variables
  integer(kind=intType) :: ierr

  ! We will use empty generic vectors w_like1, w_like2

  call VecPlaceArray(w_like1, in_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(w_like2, out_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call MatMult(dRdwT, w_like1, w_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(w_like1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(w_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

#endif
end subroutine getdRdwTVec
