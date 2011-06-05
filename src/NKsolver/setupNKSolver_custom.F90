subroutine setupNKsolver_custom

  ! Setup the PETSc objects for the Newton-Krylov
  ! solver. destroyNKsolver can be used to destroy the objects created
  ! in this function

  use communication
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: ksp,dRdw,dRdwPre,ctx,wVec,rVec,deltaW,&
       NKsolvedOnce,nksolversetup,ksp_subspace,ksp_solver_type
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw,totalCells
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  external FormFunction2

  if (not(NKSolverSetup)) then

     !  Create the linear solver context
     call KSPCreate(SUMB_PETSC_COMM_WORLD,ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     !  Create residual and state vectors
     nDimW = nw * nCellsLocal * nTimeIntervalsSpectral
     call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimw,PETSC_DETERMINE,wVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call VecSetBlockSize(wVec,nw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Use the wVec Template to create deltaW and rVec
     call VecDuplicate(wVec, rVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDuplicate(wVec, deltaW, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Setup Pre-Conditioning Matrix
     totalCells = nCellsLocal*nTimeIntervalsSpectral
     allocate( nnzDiagonal(totalCells),nnzOffDiag(totalCells))

     call drdwPCPreAllocation(nnzDiagonal,nnzOffDiag,totalCells)
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWPre, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag)

     ! Setup Matrix-Free dRdw matrix and its function
     call MatCreateMFFD(sumb_comm_world,nDimW,nDimW,&
          PETSC_DETERMINE,PETSC_DETERMINE,dRdw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatMFFDSetFunction(dRdw,FormFunction2,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the mat_row_oriented option to false so that dense
     ! subblocks can be passed in in fortran column-oriented format
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set operators for the solver
     call KSPSetOperators(ksp,dRdw,dRdWPre, DIFFERENT_NONZERO_PATTERN,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set possible options
     call KSPSetFromOptions(ksp, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set Solver Type
     call KSPSetType(ksp,ksp_solver_type,ierr);
     call EChk(ierr,__FILE__,__LINE__)

     ! Set Subspace Size
     call KSPGMRESSetRestart(ksp, ksp_subspace,ierr); 
     call EChk(ierr,__FILE__,__LINE__)

     ! Set PC Side as RIGHT only
     call KSPSetPreconditionerSide(ksp,PC_RIGHT,ierr);
     call EChk(ierr,__FILE__,__LINE__)

     NKSolverSetup = .True.
     NKSolvedOnce = .False.
  end if
end subroutine setupNKsolver_custom

