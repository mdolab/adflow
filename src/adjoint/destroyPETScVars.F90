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
  if (adjointPETScVarsAllocated) then 

     ! Matrices
     call MatDestroy(dRdWT, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     if (ApproxPC) then
        call MatDestroy(dRdWPreT, PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
     end if
     
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
        
     if (.not. useMatrixFreedRdx) then
        call MatDestroy(dRdx, PETScIerr)
        call EChk(PETScIerr,__FILE__,__LINE__)
     end if
     
     call vecDestroy(overArea, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call vecDestroy(fCell, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call vecDestroy(fNode, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     
     call KSPDestroy(adjointKSP, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     adjointPETScVarsAllocated = .False.
  end if
#endif
end subroutine destroyPETScVars

