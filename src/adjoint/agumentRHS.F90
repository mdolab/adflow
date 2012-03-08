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
#ifndef USE_NO_PETSC 
 !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the structural adjoint vector phi, by dFdw^T to       *
  !     * produce the right hand side agumentation                       *
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
  integer(kind=intType) :: ierr

  call VecPlaceArray(fVec1, phi, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Dump the result into adjointRHS
  call MatMultTranspose(dFdw,fVec1,adjointRHS,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecResetArray(fVec1, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine agumentRHS


subroutine getdFdwTVec(in_vec, in_dof, out_vec, out_dof)
#ifndef USE_NO_PETSC 
 !
  !     ******************************************************************
  !     *                                                                *
  !     * Multiply the input_vector in_vec, by dFdw^T and return the     *
  !     * resulting vector                                               *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use communication

  implicit none

  ! Input/Ouput
  integer(kind=intType), intent(in) :: in_dof, out_dof
  real(kind=realType), intent(in) :: in_vec(in_dof)
  real(kind=realType), intent(inout) :: out_vec(out_dof)

  ! Working
  integer(kind=intType) :: ierr
  
  ! Put petsc wrapper around arrays
  call VecPlaceArray(fVec1, in_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call VecPlaceArray(w_like1, out_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Dump the result into adjointRHS since the we want to ADD this
  ! result to w_like1 below
  call MatMultTranspose(dFdw,fVec1,adjointRHS,ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! do: w_like1 = w_like1 + adjointRHS
  call VecAxpy(w_like1, one, adjointRHS, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call vecResetArray(fVec1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(w_like1,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif

end subroutine getdFdwTVec




