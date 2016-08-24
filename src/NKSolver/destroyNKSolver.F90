subroutine destroyNKsolver
#ifndef USE_NO_PETSC
  ! Destroy all the PETSc objects for the Newton-Krylov
  ! solver. 

  use constants
  use NKsolverVars, only :dRdw, dRdwPre, dRdwPseudo, wVec, rVec, deltaW, g, &
       work, NK_KSP, NK_solverSetup
  use utils, only: EChk
  implicit none
  integer(kind=intType) :: ierr
  
  if (NK_solverSetup) then 
 
     call MatDestroy(dRdw, ierr) 
     call EChk(ierr, __FILE__, __LINE__)
     
     call MatDestroy(dRdwPre, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call MatDestroy(dRdwPseudo, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDestroy(wVec, ierr)  
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDestroy(rVec, ierr) 
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDestroy(deltaW, ierr)
     call EChk(ierr, __FILE__, __LINE__)
     
     call VecDestroy(g, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call VecDestroy(work, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     call KSPDestroy(NK_KSP, ierr)
     call EChk(ierr, __FILE__, __LINE__)

     NK_solverSetup = .False.
  end if
#endif
end subroutine destroyNKsolver
