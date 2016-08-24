subroutine destroyANKsolver
#ifndef USE_NO_PETSC
  ! Destroy all the PETSc objects for the Newton-Krylov
  ! solver. 

  use constants
  use ANKsolverVars, only :dRdwPre, wVec, rVec, deltaW, ANK_KSP, dRdwPreTurb, &
       wVecTurb, rVecTurb, deltaWTurb, ANK_KSPTurb, ANK_solverSetup, ANK_turbSetup
  use utils, only : EChk
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

     if (ANK_turbSetup) then 
        call MatDestroy(dRdwPreTurb, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDestroy(wVecTurb, ierr)  
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDestroy(rVecTurb, ierr) 
        call EChk(ierr, __FILE__, __LINE__)
        
        call VecDestroy(deltaWTurb, ierr)
        call EChk(ierr, __FILE__, __LINE__)
        
        call KSPDestroy(ANK_KSPTurb, ierr)
        call EChk(ierr, __FILE__, __LINE__)
     end if
     ANK_SolverSetup = .False.
     ANK_TurbSetup = .False.
  end if
#endif
end subroutine destroyANKsolver
