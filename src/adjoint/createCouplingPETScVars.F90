subroutine createCouplingPETScVars
#ifndef USE_NO_PETSC
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matries: dFdw,dFdx                                  *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc, only: dFdx, dFdw, PETScIerr
  use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use communication   ! myID, nProc
  use inputTimeSpectral !nTimeIntervalsSpectral
  use flowVarRefState
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  !
  !     Local variables.
  !
  integer(kind=intType)  :: nDimW, nDimS, nTS
  integer(kind=intType), dimension(:), allocatable :: nnzDiagonal, nnzOffDiag

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !


  nDimW = nw * nCellsLocal*nTimeIntervalsSpectral

  call getForceSize(nDimS,nTS)
  nDimS = nDimS * 3 *nTimeIntervalsSpectral! Multiply by 3 for each
                                           ! dof on each point

  ! Create dFdx and dFdw

  ! Each nodal force is contribed by (nominally) 4 quadrilateral cells
  ! surrounding it. The pressure on each of these cells depend on the
  ! average of the pressure of the cell above and the halo below. The
  ! 1st halo is computed by linear pressure extrapolation from the
  ! first two cells ABOVE the surface. The results in each node being
  ! affected by 8 cells. All of these cells on on-processors so there
  ! should be zero offdiag entries. For dFdx, the coordinates of each
  ! of the 4 quadrilateral affects the force, so this results in 9
  ! points spatial points affecting the force. Each coordiante has 3
  ! dimension which results in 3x3x3=27 nonzeros per row
  ! don't know where the non-zeros will end up

  ! Create the matrix dFdw

  allocate( nnzDiagonal(nDimS), nnzOffDiag(nDimS) )
  nnzDiagonal = 8*nw
  nnzOffDiag  = 8*nw! Make the off diagonal the same, since we
  if (PETSC_VERSION_MINOR == 2) then
     call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,&
          nDimS, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdw, PETScIerr)
  else
      call MatCreateAIJ(SUMB_PETSC_COMM_WORLD,&
          nDimS, nDimW,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdw, PETScIerr)
   end if
   call EChk(PETScIerr,__FILE__,__LINE__)

  ! Create the matrix dFdx
  nnzDiagonal = 27
  nnzOffDiag = 27
  if (PETSC_VERSION_MINOR == 2) then
     call MatCreateMPIAIJ(SUMB_PETSC_COMM_WORLD,&
          nDimS, nDimS,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdx, PETScIerr)
  else
     call MatCreateAIJ(SUMB_PETSC_COMM_WORLD,&
          nDimS, nDimS,                     &
          PETSC_DETERMINE, PETSC_DETERMINE, &
          0, nnzDiagonal,         &
          0, nnzOffDiag,            &
          dFdx, PETScIerr)
  end if

  call EChk(PETScIerr,__FILE__,__LINE__)
  deallocate( nnzDiagonal, nnzOffDiag )

  ! Set column major order for the matrix dFdw.
  call MatSetOption(dFdw, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  call MatSetOption(dFdx, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

#endif

end subroutine createCouplingPETScVars
