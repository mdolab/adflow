subroutine destroyANKsolver
#ifndef USE_NO_PETSC
  ! Destroy all the PETSc objects for the Newton-Krylov
  ! solver. 

  use ANKsolverVars
  implicit none
  integer(kind=intType) :: ierr
  
  if (ANK_SolverSetup) then
 
     call MatDestroy(dRdwPre, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDestroy(wVec, ierr)  
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDestroy(rVec, ierr) 
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDestroy(deltaW, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call KSPDestroy(ANK_KSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     ANK_SolverSetup = .False.
  end if
#endif
end subroutine destroyANKsolver
