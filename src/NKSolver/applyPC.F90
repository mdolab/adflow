subroutine applyPC(in_vec, out_vec, ndof)
#ifndef USE_NO_PETSC

  ! Apply the NK PC to the in_vec. This subroutine is ONLY used as a
  ! preconditioner for a global Aero-Structural Newton-Krylov Method

  use NKSolverVars
  use communication
  implicit none
#include "include/petscversion.h"
  ! Input/Output
  integer(kind=intType) :: ndof
  real(kind=realType), dimension(ndof), intent(in)    :: in_vec
  real(kind=realTYpe), dimension(ndof), intent(inout) :: out_vec

  ! Working Variables
  integer(kind=intType) :: ierr

  ! Setup the NKsolver if not already done so
  if (.not. NK_solverSetup) then
     call setupNKSolver
  end if
  
  ! We possibly need to re-form the jacobian
  if (mod(NK_iter, NK_jacobianLag) == 0) then 
     call FormJacobianNK()
  end if

  ! Place the two arrays into two vectos. We reuse 'work' and 'g'. 
  call VecPlaceArray(work, in_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(g, out_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
  
  ! Set the base vec
  call setwVec(wVec)
  
  call MatMFFDSetBase(dRdW, wVec, PETSC_NULL_OBJECT, ierr)
  call EChk(ierr, __FILE__, __LINE__)

   ! This needs to be a bit better...
  call KSPSetTolerances(NK_KSP, 1e-8, 1e-16, 10.0, &
       applyPCSubSpaceSize, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Actually do the Linear Krylov Solve
  call KSPSolve(NK_KSP, work, g, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Reset the array pointers:
  call VecResetArray(work, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(g, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  NK_iter = NK_iter + 1
#endif
end subroutine applyPC

subroutine applyAdjointPC(in_vec, out_vec, ndof)
#ifndef USE_NO_PETSC
  ! Apply the Adjoint PC to the in_vec. This subroutine is ONLY used as a
  ! preconditioner for a global Aero-Structural Krylov Method

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
  call VecPlaceArray(psi_like1, in_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecPlaceArray(psi_like2, out_vec, ierr)
  call EChk(ierr, __FILE__, __LINE__)
 
  ! Set KSP_NORM Type to none. Implictly turns off convergence
  ! check. Since we just want to run a fixed number of iterations this
  ! is fine. The should be set regardless of the KSPType.

  call KSPSetNormType(adjointKSP, KSP_NORM_NONE, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! This needs to be a bit better...
  call KSPSetTolerances(adjointKSP, PETSC_DEFAULT_REAL, &
       PETSC_DEFAULT_REAL, PETSC_DEFAULT_REAL, &
       applyAdjointPCSubSpaceSize, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Actually do the Linear Krylov Solve
  call KSPSolve(adjointKSP, psi_like1, psi_like2, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  ! Reset the array pointers:
  call VecResetArray(psi_like1, ierr)
  call EChk(ierr, __FILE__, __LINE__)

  call VecResetArray(psi_like2, ierr)
  call EChk(ierr, __FILE__, __LINE__)
#endif
end subroutine applyAdjointPC
