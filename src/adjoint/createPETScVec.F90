!
!     ******************************************************************
!     *                                                                *
!     * File:          createPETScVec.F90                              *
!     * Author:        C.A.(Sandy) Mader, Andre C. Marta               *
!     * Starting date: 01-14-2008                                      *
!     * Last modified: 02-02-2008                                      *
!     *                                                                *
!     ******************************************************************
!
subroutine createPETScVec
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the global RHS vector of the adjoint system of          *
  !     * equations, the solution vector and the residual vector         *
  !     * as PETSc vector objects, which are used to compute the adjoint *
  !     * solution. Also, create the vectors of partial and total        *
  !     * function gradients - dJda,dJdx - and - dIda,dIdx - respectively*
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use ADjointVars     ! nCellsLocal,
  use communication   ! myID, nProc 
  use inputTimeSpectral !nTimeIntervalsSpectral
  use flowVarRefState ! 
  use mdData          ! mdNSurfNodesCompact
  implicit none
 

#ifndef USE_NO_PETSC

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Vectors of length nDimw
  !     *                                                                *
  !     ******************************************************************

  call MatGetVecs(dRdwT,dJdW,PETSC_NULL_OBJECT,PETScIerr)
  call EChk(PETScIerr,__file__,__line__)

  call VecDuplicate(dJdW, psi, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)
  call VecSet(psi,PETScZero,PETScIerr)
  call EChk(PETScIerr,__file__,__line__)

  call VecDuplicate(dJdW, pvr, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Vectors of length nDimS
  !     *                                                                *
  !     ******************************************************************

  call MatGetVecs(dFdx,dJdx,PETSC_NULL_OBJECT,PETScIerr)
  call EChk(PETScIerr,__file__,__line__)


#endif




end subroutine createPETScVec
