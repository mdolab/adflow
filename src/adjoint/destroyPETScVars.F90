! there are a number of differenet functions in destroyPETScVars. They
! coorespond to the different creation functions:

subroutine destroyStatePETScVars

  use ADjointPETSc
  use inputADjoint    !ApproxPC
  implicit none

#ifndef USE_NO_PETSC

  ! Matrices
  call MatDestroy(dRdWT, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  if (ApproxPC) then
     call MatDestroy(dRdWPreT, PETScIerr)
     call EChk(PETScIerr,__FILE__,__LINE__)
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
 
#endif

end subroutine destroyStatePETScVars

subroutine destroySpatialPETScVars

  use ADjointPETSc
  implicit none

#ifndef USE_NO_PETSC

  call MatDestroy(dRdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecDestroy(dJdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call VecDestroy(xVec,PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)


#endif

end subroutine destroySpatialPETScVars

subroutine destroyCouplingPETScVars

  use ADjointPETSc
  implicit none
  
#ifndef USE_NO_PETSC

  call MatDestroy(dFdw, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  call MatDestroy(dFdx, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

#endif

end subroutine destroyCouplingPETScVars

subroutine destroyExtraPETScVars
  use ADjointPETSc
  implicit none
#ifndef USE_NO_PETSC

  call MatDestroy(dRda, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

  deallocate(drda_data)

#endif
  
end subroutine destroyExtraPETScVars

subroutine destroyPETScKSP

  use ADjointPETSc
  implicit none
#ifndef USE_NO_PETSC
  ! KSP Context
  call KSPDestroy(adjointKSP, PETScIerr)
  call EChk(PETScIerr,__FILE__,__LINE__)

#endif

end subroutine destroyPETScKSP

