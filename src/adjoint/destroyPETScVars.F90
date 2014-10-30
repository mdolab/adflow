! there are a number of differenet functions in destroyPETScVars. They
! coorespond to the different creation functions:

subroutine destroyPETScVars

  use ADjointPETSc
  use inputADjoint    
  use blockPointers
  use costFunctions
  use inputTimeSpectral
  implicit none

  integer(kind=intType) :: i, nLevels, sps
  
#ifndef USE_NO_PETSC

  ! Matrices
  call MatDestroy(dRdWT, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dRdwTShell, PETScIErr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  if (ApproxPC) then
     call MatDestroy(dRdWPreT, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
  end if

  nLevels = ubound(flowDoms, 2)
  if (preCondType == 'mg') then
     do i=2,nLevels
        call MatDestroy(coarsedRdwPreT(i), PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
        
        call MatDestroy(restrictionOperator(i), PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)

        call MatDestroy(prolongationOperator(i), PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)

     end do
     deallocate(coarsedRdwPreT, restrictionOperator, prolongationOperator)
  end if

  ! Vectors
  call VecDestroy(psi, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecDestroy(dJdW, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  
  call VecDestroy(adjointRHS, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  
  call VecDestroy(adjointRes, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dFcdw, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dFcdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(doAdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dFndFc, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dFdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dFdw, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  do sps=1,nTimeIntervalsSpectral
     do i=1,nFM
        call VecDestroy(FMw(i, sps), PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
     end do
  end do
  deallocate(FMw)

  call MatDestroy(dRdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecDestroy(dJdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecDestroy(xVec,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  do sps=1,ntimeIntervalsSpectral
     do i=1,nFM
        call VecDestroy(FMx(i, sps), PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
     end do
  end do
  deallocate(FMx)

  call vecDestroy(overArea, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call vecDestroy(fCell, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call vecDestroy(fNode, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dRda, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
  
  if (allocated(dRda_data)) then
     deallocate(dRda_data)
  end if
  
  call KSPDestroy(adjointKSP, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)
#endif
end subroutine destroyPETScVars

