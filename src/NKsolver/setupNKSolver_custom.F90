subroutine setupNKsolver_custom

  ! Setup the PETSc objects for the Newton-Krylov
  ! solver. destroyNKsolver can be used to destroy the objects created
  ! in this function

  use communication
  use inputTimeSpectral
  use flowVarRefState
  use iteration
  use inputPhysics
  use stencils
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: ksp,dRdw,dRdwPre,dRdwPseudo,ctx,&
       wVec,rVec,diagV,deltaW,&
       NKsolvedOnce,nksolversetup,ksp_subspace,ksp_solver_type

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw,totalCells
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType) :: n_stencil
  integer(kind=intType), dimension(:,:), allocatable :: stencil

  external FormFunction2,mymult

  if (not(NKSolverSetup)) then
     
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

     call VecDuplicate(wVec, diagV, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Setup Pre-Conditioning Matrix
     totalCells = nCellsLocal*nTimeIntervalsSpectral
     allocate( nnzDiagonal(totalCells),nnzOffDiag(totalCells))

     call initialize_stencils
     if (not(viscous)) then
        n_stencil = N_euler_drdw
        allocate(stencil(n_stencil,3))
        stencil = euler_drdw_stencil
     else
        n_stencil = N_visc_pc
        allocate(stencil(n_stencil,3))
        stencil = visc_pc_stencil
     end if

     call statePreAllocation(nnzDiagonal,nnzOffDiag,nDimW/nw,stencil,n_stencil)
  
     call MatCreateMPIBAIJ(SUMB_PETSC_COMM_WORLD, nw,             &
          nDimW, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dRdWPre, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag,stencil)

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
     
     ! Create the matrix shell for the pseudo-transient part:
     call MatCreateShell(sumb_comm_world,nDimW,nDimW,&
          PETSC_DETERMINE,PETSC_DETERMINE,ctx,dRdwPseudo,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Set the mat-mult operation...that's all we need
     call MatShellSetOperation(dRdwPseudo,MATOP_MULT,mymult,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     !  Create the linear solver context
     call KSPCreate(SUMB_PETSC_COMM_WORLD,ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set operators for the solver
     !call KSPSetOperators(ksp,dRdw,dRdWPre, DIFFERENT_NONZERO_PATTERN,ierr)
     !call EChk(ierr,__FILE__,__LINE__)

     call KSPSetOperators(ksp,dRdwPseudo,dRdWPre, &
          DIFFERENT_NONZERO_PATTERN,ierr)
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
     call KSPSetPCSide(ksp,PC_RIGHT,ierr);
     call EChk(ierr,__FILE__,__LINE__)
   
     NKSolverSetup = .True.
     NKSolvedOnce = .False.
     if(equations == RANSEquations)  then
        turbCoupled = .True.
     end if
  end if
end subroutine setupNKsolver_custom


subroutine MyMult(matrix,X,F,ierr)

  !   Input Parameters:
  !.  X - input vector
  !
  !   Output Parameter:
  !.  F - function vector
  !
  use precision
  use NKSolverVars, only: dRdw,diagV
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
  
  ! Input/Output Vars
  Mat matrix
  Vec X,F
  integer(kind=intType) :: ierr

  ! Do a matmult followed by an addition

  call MatMult(dRdw,X,F,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! VecAXPY : Computes y = alpha x + y. 
  ! VecAXPY(Vec y,PetscScalar alpha,Vec x)

  !call VecAXPY(F,1.0,diagV,ierr)
  !call EChk(ierr,__FILE__,__LINE__)

end subroutine MyMult

