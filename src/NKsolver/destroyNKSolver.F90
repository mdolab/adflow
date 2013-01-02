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
     
     call VecDestroy(wBase,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(rBase,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(scaleVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(g,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(w_like1,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(w_like2,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(work,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call KSPDestroy(newtonKrylovKSP, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     NKSolverSetup = .False.
  end if
#endif
end subroutine destroyNKsolver

subroutine destroyNKsolver2
#ifndef USE_NO_PETSC
  ! Destroy all the PETSc objects for the Newton-Krylov
  ! solver. 

  use NKsolverVars
  implicit none
  integer(kind=intType) :: ierr
  
  if (NKSolverSetup) then
 
     call VecDestroy(wVec,ierr)  
     call EChk(ierr,__FILE__,__LINE__)
     
     call VecDestroy(rVec,ierr) 
     call EChk(ierr,__FILE__,__LINE__)

     call VecDestroy(rhs,ierr) 
     call EChk(ierr,__FILE__,__LINE__)
     
     call snesDestroy(snes, ierr)
     call EChk(ierr,__FILE__,__LINE__)

     call snesDestroy(psnes, ierr)
     call EChk(ierr,__FILE__,__LINE__)
     NKSolverSetup = .False.
  end if
#endif
end subroutine destroyNKsolver2
