subroutine createExtraPETScVars
  !
  !     ******************************************************************
  !     *                                                                *
  !     * Create the matries: dRda                                       *
  !     *                                                                *
  !     ******************************************************************
  !
  use ADjointPETSc
  use ADjointVars     ! nCellsLocal,nNodesLocal, nDesignExtra
  use communication   ! myID, nProc
  use inputTimeSpectral !nTimeIntervalsSpectral
  use flowvarrefstate
  implicit none
  !
  !     Local variables.
  !
  integer(kind=intType) :: nDimW

  !
  !     ******************************************************************
  !     *                                                                *
  !     * Begin execution.                                               *
  !     *                                                                *
  !     ******************************************************************
  !
#ifndef USE_NO_PETSC

  nDimW = nw * nCellsLocal*nTimeIntervalsSpectral
 
  ! dRda
  call MatCreateMPIDense(SUMB_PETSC_COMM_WORLD,nDimW,PETSC_DECIDE,&
       PETSC_DETERMINE,nDesignExtra,PETSC_NULL_SCALAR,dRda,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  ! Set column major order for the matrix dRda.
  call MatSetOption(dRda, MAT_ROW_ORIENTED,PETSC_TRUE, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

#endif
end subroutine createExtraPETScVars
 
