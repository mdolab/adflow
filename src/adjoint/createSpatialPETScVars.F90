subroutine createSpatialPETScVars
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matrices/vectors that are required for the total    *
  !     * derivative: Matrices: dRdx                                     *
  !     *              Vectors: dJdx                                     *
  !     *                                                                *
  !     ******************************************************************

  use ADjointPETSc, only : dRdx, PETScIerr, dJdx, Xvec
  use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use communication   ! myID, nProc
  use inputTimeSpectral !nTimeIntervalsSpectral
  use flowVarRefState ! 
  use inputADjoint    !ApproxPC

  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  !
  !     Local variables.
  !
  integer(kind=intType) :: nDimW, nDimX
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
#ifndef USE_NO_PETSC

  nDimW = nw * nCellsLocal*nTimeIntervalsSpectral
  nDimX = 3 * nNodesLocal*nTimeIntervalsSpectral

  !     ******************************************************************
  !     *                                                                *
  !     * Create matrix dRdx that is used to compute the total cost /    *
  !     * constraint function sensitivity with respect to the spatial    *
  !     * design variables 'x' as dIdx = dJdx - psi^T dRdx.              *
  !     *                                                                *
  !     * Matrix dRdx has size [nDimW,nDimX] and is generally            *
  !     * sparse for the coordinate design variables.                    *
  !     *                                                                *
  !     * The local dimensions are specified so that the spatial         *
  !     * coordinates x (a) are placed in the local processor. This has  *
  !     * to be consistent with the vectors dIdx and dJdx.               *
  !     *                                                                *
  !     ******************************************************************

  allocate( nnzDiagonal(nDimX), nnzOffDiag(nDimX) )
  ! Create the matrix dRdx.
  call drdxPreAllocation(nnzDiagonal,nnzOffDiag,nDimX)
  
  ! Note we are creating the TRANPOSE of dRdx. It is size dDimX by nDimW
  if (PETSC_VERSION_MINOR == 2) then
     call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,                 &
          nDimX, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          8, nnzDiagonal,     &
          8, nnzOffDiag,            &
          dRdx, PETScIerr)
  else
     call MatCreateAIJ(SUMB_PETSC_COMM_WORLD,                 &
          nDimX, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          8, nnzDiagonal,     &
          8, nnzOffDiag,            &
          dRdx, PETScIerr)
  end if
  call EChk(PETScIerr,__FILE__,__LINE__)

  deallocate( nnzDiagonal, nnzOffDiag )

  ! Set column major order for the matrix dRdx.
  call MatSetOption(dRdx, MAT_ROW_ORIENTED,PETSC_FALSE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! Vectors
  call VecCreate(SUMB_COMM_WORLD, dJdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  
  call VecSetSizes(dJdx, nDimX, PETSC_DECIDE, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)
  
  call VecSetBlockSize(dJdx, 3, PETScIerr)
  call EChk(PETScIerr, __FILE__, __LINE__)
  
  call VecSetType(dJdx, VECMPI, PETScIerr) 
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! xVec
  call VecDuplicate(dJdx,xVec,PETScIerr)
  call EChk(PETScierr,__FILE__,__LINE__)

#endif
end subroutine createSpatialPETScVars
