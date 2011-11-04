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

  ! Once again, PETSC is royally screwed up. You CANNOT use PETSC_NULL
  ! arguments. They simply do NOT work in Fortran. The PETSc
  ! documentation lies to you. We have to allocate our own data. 

  allocate(dRda_data(nDimw,nDesignExtra))


  call MatCreateMPIDense(SUMB_PETSC_COMM_WORLD,nDimW,PETSC_DECIDE,&
       PETSC_DETERMINE,nDesignExtra,dRda_data,dRda,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

#endif
end subroutine createExtraPETScVars
 
