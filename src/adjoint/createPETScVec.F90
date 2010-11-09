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
  use ADjointVars     ! nCellsLocal, nDesignExtra
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

  !     ******************************************************************
  !     *                                                                *
  !     * Create the vectors of partial and total function gradient      *
  !     * with respect to the extra design variables dJda and dIda,      *
  !     * respectively. Vector dJda has size [nDesignExtra].             *
  !     *                                                                *
  !     * The global dimension is specified so that the extra design     *
  !     * variables are scattered among the processors.                  *
  !     * This has to be consistent with the matrix dRda.                *
  !     *                                                                *
  !     ******************************************************************
  !
  call VecCreate(SUMB_PETSC_COMM_WORLD, dJda, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)
  call VecSetSizes(dJda, PETSC_DECIDE, &
       nDesignExtra, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)
  call VecSetType(dJda,"mpi",PETScIerr)
  call EChk(PETScIerr,__file__,__line__)

  call VecDuplicate(dJda, dIda, PETScIerr)
  call EChk(PETScIerr,__file__,__line__)

#endif

end subroutine createPETScVec
