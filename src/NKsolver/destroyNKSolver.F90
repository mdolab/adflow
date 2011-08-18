subroutine destroyNKsolver_petsc

  ! Destroy all the PETSc objects for the Newton-Krylov
  ! solver. 

  use NKsolverVars
  implicit none
  integer(kind=intType) :: ierr
  
  if (NKSolverSetup) then

     ! We will destroy the PETSc variables created in setupNKsolver
     call SNESDestroy(snes,ierr) ! Also destroys the underlying ksp
     ! and pc contexts
     call EChk(ierr,__FILE__,__LINE__)

     call MatDestroy(dRdw,ierr);      call EChk(ierr,__FILE__,__LINE__)
     call MatDestroy(dRdwPre,ierr);   call EChk(ierr,__FILE__,__LINE__)
     call VecDestroy(wVec,ierr);      call EChk(ierr,__FILE__,__LINE__)
     call VecDestroy(rVec,ierr);      call EChk(ierr,__FILE__,__LINE__)

     NKSolverSetup = .False.
  end if
  
end subroutine destroyNKsolver_petsc


subroutine destroyNKsolver

  ! Destroy all the PETSc objects for the Newton-Krylov
  ! solver. 

  use NKsolverVars
  implicit none
  integer(kind=intType) :: ierr
  
  if (NKSolverSetup) then

     ! We will destroy the PETSc variables created in setupNKsolver
     call KSPDestroy(ksp,ierr) ! Also destroys the underlying PC

     ! and pc contexts
     call EChk(ierr,__FILE__,__LINE__)

     call MatDestroy(dRdw,ierr);      call EChk(ierr,__FILE__,__LINE__)
     call MatDestroy(dRdwPre,ierr);   call EChk(ierr,__FILE__,__LINE__)
     call VecDestroy(wVec,ierr);      call EChk(ierr,__FILE__,__LINE__)
     call VecDestroy(rVec,ierr);      call EChk(ierr,__FILE__,__LINE__)

     NKSolverSetup = .False.
  end if
  
end subroutine destroyNKsolver
