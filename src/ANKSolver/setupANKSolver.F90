subroutine setupANKsolver
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
  use ADjointVars , only: nCellsLocal
  use ANKSolverVars, only: dRdwPre, wVec, rVec, deltaW, ANK_solverSetup, ANK_KSP, ANK_iter
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

  ! Make sure we don't have memory for the approximate and exact
  ! Newton solvers kicking around at the same time.
  call destroyNKSolver()

  if (.not. ANK_solverSetup) then
     nDimW = nwf * nCellsLocal(1_intTYpe) * nTimeIntervalsSpectral

     call VecCreate(SUMB_COMM_WORLD, wVec, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecSetSizes(wVec, nDimW, PETSC_DECIDE, ierr)
     call EChk(ierr, __FILE__, __LINE__)
 
     call VecSetBlockSize(wVec, nwf, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecSetType(wVec, VECMPI, ierr) 
     call EChk(ierr, __FILE__, __LINE__)

     !  Create duplicates for residual and delta
     call VecDuplicate(wVec, rVec, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecDuplicate(wVec, deltaW, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Create Pre-Conditioning Matrix
     allocate(nnzDiagonal(nCellsLocal(1_intType)*nTimeIntervalsSpectral), &
          nnzOffDiag(nCellsLocal(1_intType)*nTimeIntervalsSpectral) )

     stencil => euler_pc_stencil
     n_stencil = N_euler_pc

     level = 1
     call statePreAllocation(nnzDiagonal, nnzOffDiag, nDimW/nwf, stencil, n_stencil, &
          level)
     call myMatCreate(dRdwPre, nwf, nDimW, nDimW, nnzDiagonal, nnzOffDiag, &
          __FILE__, __LINE__)

     call matSetOption(dRdwPre, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     deallocate(nnzDiagonal, nnzOffDiag)
     
     ! Set the mat_row_oriented option to false so that dense
     ! subblocks can be passed in in fortran column-oriented format
     call MatSetOption(dRdWPre, MAT_ROW_ORIENTED, PETSC_FALSE, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     !  Create the linear solver context
     call KSPCreate(SUMB_COMM_WORLD, ANK_KSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ! Set operators for the solver

     call KSPSetOperators(ANK_KSP, dRdwPre, dRdwPre, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ANK_solverSetup = .True.
     ANK_iter = 0
  end if
#endif

end subroutine setupANKsolver

