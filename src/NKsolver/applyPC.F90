subroutine applyPC(in_vec, out_vec, N)
#ifndef USE_NO_PETSC
  ! Apply the NK PC to the in_vec. This subroutine is ONLY used as a
  ! preconditioner for a global Aero-Structural Newton-Krylov Method
  use flowVarRefState
  use NKSolverVars, only: dRdw, NKSolverSetup, global_ksp, wVec, &
       NKSolveCount, jacobian_lag
  use communication 
  use inputIteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input/Output
  integer(kind=intType) :: N
  real(kind=realType), dimension(N), intent(in) :: in_vec(N)
  real(kind=realTYpe), dimension(N), intent(out):: out_vec(N)
  
  ! PETSc 
  Vec VecA,VecB

  ! Working Variables
  integer(kind=intType) :: ierr

  ! Put a petsc wrapper around the input and output vectors
  call VecCreateMPIWithArray(sumb_comm_world, N, PETSC_DETERMINE, in_vec, &
       VecA, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecSetBlockSize(vecA, nw, ierr);
  call EChk(ierr,__FILE__,__LINE__)

  call VecCreateMPIWithArray(sumb_comm_world, N, PETSC_DETERMINE, out_vec, &
       VecB, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecSetBlockSize(vecB, nw, ierr);
  call EChk(ierr,__FILE__,__LINE__)

  if (not(NKSolverSetup)) then
     call setupNKSolver
  end if
  
  if (mod(NKsolveCount,jacobian_lag) == 0) then
     call FormJacobian()
  end if

  ! Set the base vec
  call setwVec(wVec)
  
  call MatMFFDSetBase(dRdW, wVec, PETSC_NULL_OBJECT, ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  call KSPSetTolerances(global_ksp,.1,.00000000001,10.0,10,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Actually do the Linear Krylov Solve
  call KSPSolve(global_ksp, vecA, vecB,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Destroy the three petsc vectors
  call VecDestroy(VecA, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecDestroy(VecB, ierr)
  call EChk(ierr,__FILE__,__LINE__)

  NKSolveCount = NKSolveCount + 1
#endif
end subroutine applyPC
