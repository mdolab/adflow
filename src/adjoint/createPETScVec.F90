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
  integer(kind=intType) :: ndimS, nTS

#ifndef USE_NO_PETSC

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Vectors of length nDimw
  !     *                                                                *
  !     ******************************************************************

  ! Get dJdw and psi from one MatGetVecs Call
  call MatGetVecs(dRdwT,dJdW,psi,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! adjointRes is the same size as dJdw,psi
  call VecDuplicate(dJdW, adjointRes, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecDuplicate(dJdW, adjointRHS, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  !     ******************************************************************
  !     *                                                                *
  !     * Vectors of length nDimS
  !     *                                                                *
  !     ******************************************************************

  call getForceSize(nDimS,nTS)
  nDimS = nDimS * 3 *nTimeIntervalsSpectral! Multiply by 3 for each
                                           ! dof on each point

  call VecCreateMPI(SUMB_PETSC_COMM_WORLD,nDimS,PETSC_DETERMINE,dJdx,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecSetBlockSize(dJdx,3,PETScIerr)
  call EChk(PETScierr,__FILE__,__LINE__)
 
#endif

end subroutine createPETScVec
