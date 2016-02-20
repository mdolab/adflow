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
  use NKSolverVars, only: dRdw, dRdwPre, ctx, wVec, rVec, deltaW, &
       NKsolvecount, nksolversetup, work, g, NKViscPC, newtonKrylovKSP
  use ADJointPetsc,  only: drdwpret
  implicit none
#define PETSC_AVOID_MPIF_H

#include "include/petscversion.h"
#if PETSC_VERSION_MINOR > 5
#include "petsc/finclude/petsc.h"
#else
#include "include/finclude/petsc.h"
#endif

  ! Working Variables
  integer(kind=intType) :: ierr, nDimw
  integer(kind=intType) , dimension(:), allocatable :: nnzDiagonal, nnzOffDiag
  integer(kind=intType) :: n_stencil
  integer(kind=intType), dimension(:, :), pointer :: stencil
  integer(kind=intType) :: level
  external FormFunction_mf

  if (.not. NKSolverSetup) then
     nDimW = nw * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

     call VecCreate(SUMB_COMM_WORLD, wVec, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
     call EChk(ierr, __FILE__, __LINE__)
 
     call VecSetBlockSize(wVec, nw, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecSetType(wVec, VECMPI, ierr) 
     call EChk(ierr, __FILE__, __LINE__)

     !  Create duplicates for residual and delta
     call VecDuplicate(wVec, rVec, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecDuplicate(wVec, deltaW, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Create the two additional work vectors for the line search:
     call VecDuplicate(wVec, g, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDuplicate(wVec, work, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Create Pre-Conditioning Matrix
     allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
          nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

     if (viscous .and. NKviscPC) then
        stencil => visc_pc_stencil
        n_stencil = N_visc_pc
     else
        stencil => euler_pc_stencil
        n_stencil = N_euler_pc
     end if

     level = 1
     call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nw, stencil, n_stencil, &
          level, .False.)
     call myMatCreate(dRdwPre, nw, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
          __FILE__, __LINE__)

     call matSetOption(dRdwPre, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     deallocate(nnzDiagonal, nnzOffDiag)

     ! Setup Matrix-Free dRdw matrix and its function
     call MatCreateMFFD(sumb_comm_world, nDimW, nDimW, &
          PETSC_DETERMINE, PETSC_DETERMINE, dRdw, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call MatMFFDSetFunction(dRdw, FormFunction_mf, ctx, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Set the mat_row_oriented option to false so that dense
     ! subblocks can be passed in in fortran column-oriented format
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call MatSetOption(dRdW, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     !  Create the linear solver context
     call KSPCreate(SUMB_COMM_WORLD, newtonKrylovKSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Set operators for the solver

     call KSPSetOperators(newtonKrylovKSP, dRdw, dRdwPre, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     NKSolverSetup = .True.
     NKSolveCount = 0
  end if
#endif

end subroutine setupNKsolver
