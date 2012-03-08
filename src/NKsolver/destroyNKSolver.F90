subroutine destroyNKsolver
#ifndef USE_NO_PETSC
  ! Destroy all the PETSc objects for the Newton-Krylov
  ! solver. 

  use NKsolverVars
  implicit none
  integer(kind=intType) :: ierr
  
  if (NKSolverSetup) then
 
     call MatDestroy(dRdw,ierr) 
     call EChk(ierr,__FILE__,__LINE__)
     
     call MatDestroy(dRdwPre,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(wVec,ierr)  
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(rVec,ierr) 
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(deltaW,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(g,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(w_like1,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(w_like2,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(work,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call KSPDestroy(global_ksp,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     NKSolverSetup = .False.
  end if
#endif
end subroutine destroyNKsolver
