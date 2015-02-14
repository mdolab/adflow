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
     
     call KSPDestroy(adjointKSP, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
     adjointPETScVarsAllocated = .False.
  end if
#endif
end subroutine destroyPETScVars

