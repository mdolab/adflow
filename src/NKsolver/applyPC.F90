subroutine applyPC(in_vec, out_vec, ndof)
#ifndef USE_NO_PETSC

  ! Apply the NK PC to the in_vec. This subroutine is ONLY used as a
  ! preconditioner for a global Aero-Structural Newton-Krylov Method

  use NKSolverVars
  use communication
  implicit none

  ! Input/Output
  integer(kind=intType) :: ndof
  real(kind=realType), dimension(ndof),intent(in)    :: in_vec
  real(kind=realTYpe), dimension(ndof),intent(inout) :: out_vec

  ! Working Variables
  integer(kind=intType) :: ierr

  ! Setup the NKsolver if not already done so
  if (.not. NKSolverSetup) then
     call setupNKSolver
  end if
  
  ! We possibly need to re-form the jacobian
  if (mod(NKsolveCount,jacobian_lag) == 0) then
     call FormJacobian()
  end if

  ! Place the two arrays in the vector
  call VecPlaceArray(w_like1, in_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(w_like2, out_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  ! Set the base vec
  call setwVec(wVec)
  
  call MatMFFDSetBase(dRdW, wVec, PETSC_NULL_OBJECT, ierr)
  call EChk(ierr,__FILE__,__LINE__)

   ! This needs to be a bit better...
  call KSPSetTolerances(newtonKrylovKSP, 1e-8,1e-16,10.0,applyPCSubSpaceSize,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Actually do the Linear Krylov Solve
  call KSPSolve(newtonKrylovKSP, w_like1, w_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Reset the array pointers:
  call VecResetArray(w_like1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(w_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  NKSolveCount = NKSolveCount + 1
#endif
end subroutine applyPC

subroutine applyAdjointPC(in_vec, out_vec, ndof)
#ifndef USE_NO_PETSC
  ! Apply the Adjoint PC to the in_vec. This subroutine is ONLY used as a
  ! preconditioner for a global Aero-Structural Newton-Krylov Method

  use communication
  use ADjointPETSc
  use inputAdjoint
  implicit none

  ! Input/Output
  integer(kind=intType) :: ndof
  real(kind=realType), dimension(ndof), intent(in)    :: in_vec
  real(kind=realTYpe), dimension(ndof), intent(inout) :: out_vec
  
  ! Working Variables
  integer(kind=intType) :: ierr

  ! Hijack adjoint and adjointRes with in_vec and out_vec
  call VecPlaceArray(w_like1, in_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecPlaceArray(w_like2, out_vec, ierr)
  call EChk(ierr,__FILE__,__LINE__)
 
  ! This needs to be a bit better...
  call KSPSetTolerances(adjointKSP,1e-8,1e-16,10.0,applyAdjointPCSubSpaceSize,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Actually do the Linear Krylov Solve
  call KSPSolve(adjointKSP, w_like1, w_like2,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Reset the array pointers:
  call VecResetArray(w_like1, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecResetArray(w_like2, ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine applyAdjointPC
