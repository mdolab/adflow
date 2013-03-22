subroutine setupNKsolver
#ifndef USE_NO_PETSC
  ! Setup the PETSc objects for the Newton-Krylov
  ! solver. destroyNKsolver can be used to destroy the objects created
  ! in this function

  use blockPointers
  use communication
  use inputTimeSpectral
  use flowVarRefState
  use iteration
  use inputPhysics
  use stencils
  use InputAdjoint, only: viscPC
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: dRdw, dRdwPre,dRdwPseudo, ctx, wVec, rVec, deltaW,&
       NKsolvecount, nksolversetup, work, g, w_like1, w_like2, scaleVec, &
       newtonKrylovKSP

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: ierr,nDimw,totalCells
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType) :: n_stencil
  integer(kind=intType), dimension(:,:), pointer :: stencil
  integer(kind=intType) :: i,j,k,nn,i2,j2,k2,d2,l,sps, level
  external FormFunction_mf
  external mykspmonitor
  if (.not. NKSolverSetup) then
     nDimW = nw * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

     call VecCreate(SUMB_COMM_WORLD, wVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
     call EChk(ierr, __FILE__, __LINE__)
 
     call VecSetBlockSize(wVec, nw, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecSetType(wVec, VECMPI, ierr) 
     call EChk(ierr,__FILE__,__LINE__)

     !  Create duplicates for residual and delta
     call VecDuplicate(wVec, rVec, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDuplicate(wVec, deltaW, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Create the two additional work vectors for the line search:
     call VecDuplicate(wVec,g,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDuplicate(wVec,work,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDuplicate(wVec,scaleVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Create two empty w-like vectors
     if (PETSC_VERSION_MINOR == 2) then
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nDimW,PETSC_DETERMINE,&
             PETSC_NULL_SCALAR,w_like1,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nDimW,PETSC_DETERMINE,&
             PETSC_NULL_SCALAR,w_like2,ierr)
        call EChk(ierr,__FILE__,__LINE__)
     else 
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nw,nDimW,PETSC_DETERMINE,&
             PETSC_NULL_SCALAR,w_like1,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call VecCreateMPIWithArray(SUMB_COMM_WORLD,nw,nDimW,PETSC_DETERMINE,&
             PETSC_NULL_SCALAR,w_like2,ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     ! Create Pre-Conditioning Matrix
     totalCells = nCellsLocal(1_intType)*nTimeIntervalsSpectral
     allocate(nnzDiagonal(totalCells), nnzOffDiag(totalCells))

     call initialize_stencils
     if (viscous .and. viscPC) then
        n_stencil = N_visc_pc
        stencil => visc_pc_stencil
     else
        n_stencil = N_euler_pc
        stencil => euler_pc_stencil
     end if

     ! Note: Since we are using blocked matrices, we computed the
     ! non-zero BLOCKS as opposed to the non-zero values. 
     level = 1_intType
     call statePreAllocation(nnzDiagonal, nnzOffDiag, &
          totalCells, stencil ,n_stencil, level)
     if (PETSC_VERSION_MINOR == 2) then
        call MatCreateMPIBAIJ(SUMB_COMM_WORLD, nw, &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,                   &
             0, nnzOffDiag,                    &
             dRdWPre, ierr)
     else
        call MatCreateBAIJ(SUMB_COMM_WORLD, nw, &
             nDimW, nDimW,                     &
             PETSC_DETERMINE, PETSC_DETERMINE, &
             0, nnzDiagonal,                   &
             0, nnzOffDiag,                    &
             dRdWPre, ierr)
     end if

     call EChk(ierr,__FILE__,__LINE__)
     
     deallocate(nnzDiagonal,nnzOffDiag)

     ! Setup Matrix-Free dRdw matrix and its function
     call MatCreateMFFD(sumb_comm_world,nDimW,nDimW,&
          PETSC_DETERMINE,PETSC_DETERMINE,dRdw,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call MatMFFDSetFunction(dRdw,FormFunction_mf,ctx,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the mat_row_oriented option to false so that dense
     ! subblocks can be passed in in fortran column-oriented format
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatSetOption(dRdW   , MAT_ROW_ORIENTED,PETSC_FALSE, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     !  Create the linear solver context
     call KSPCreate(SUMB_COMM_WORLD, newtonKrylovKSP, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set operators for the solver
     call KSPSetOperators(newtonKrylovKSP,dRdw,dRdWPre, &
          DIFFERENT_NONZERO_PATTERN,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! DEBUGGING ONLY!
!      call KSPMonitorSet(newtonKrylovKSP,MyKSPMonitor, PETSC_NULL_OBJECT, &
!          PETSC_NULL_FUNCTION, ierr)
!      call EChk(ierr,__FILE__,__LINE__)
    
     NKSolverSetup = .True.
     NKSolveCount = 0
     !if(equations == RANSEquations)  then
     !   turbCoupled = .True.
     !   turbSegregated = .False.
     !end if

  end if
#endif
end subroutine setupNKsolver
