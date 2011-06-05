subroutine snes_monitor(snes,its,norm,ctx,ierr)
  use communication
  use precision 
  use iteration
  use inputIteration
  use NKsolverVars, only: ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,&
       snes_atol,itertot0,jacobian_lag,ksp_subspace
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  SNES snes
  KSP  ksp
  PetscInt its, ierr
  PetscReal norm
  PetscFortranAddr ctx(*) ! This is probably going to be empty

  integer(kind=intType) :: ksp_its,temp,nfuncs
  real(kind=realType) :: CFLNew,rhoRes,totalRRes

  call SNESGetKSP(snes,ksp,ierr);  call EChk(ierr,__FILE__,__LINE__)
  call KSPGetTolerances(ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ierr)
  ! Set the ksp_atol to snes_atol...this lets the ksp exit early if we
  ! have hit are desired non-linear tolerance
  ksp_atol = snes_atol

  ! Set the maximum iterations ot the subspace size...not really much
  ! point to keep going...might as well do another non-linear iteration
  ksp_max_it = ksp_subspace

  call KSPSetTolerances(ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ierr)

  ! We want to get the number of function evals 
  ! its == 1 AFTER the first iteration. Take this opportunity to set
  ! the preconditioner lag to the value we actually want
  if (its == 1) then
     ! Reset the value of the preconditionerLag to what we actually want. It
     ! had been set to -1 or -2 depending on if we wanted to recompute
     ! the preconditioner on the first entry or not. 

     call SNESSetLagJacobian(snes, jacobian_lag, ierr); call EChk(ierr,__FILE__,__LINE__)

  end if

  if (its > 0 .or. iterTot0 == 0) then
     call SNESGetNumberFunctionEvals(snes, nfuncs,ierr);  call EChk(ierr,__FILE__,__LINE__)  

     nFuncs = max(nfuncs,1)
     iterTot = iterTot0 + nFuncs

     call convergenceInfo

     ierr = 0
  end if
end subroutine snes_monitor


subroutine ts_monitor(pts,steps,t,wVec,ctx,ierr)

!subroutine snes_monitor(snes,its,norm,ctx,ierr)

  use communication
  use precision 
  use iteration
  use inputIteration
  use NKsolverVars, only: ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,&
       snes_atol,itertot0,jacobian_lag,ksp_subspace
  implicit none

#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
#include "include/finclude/petscts.h"
  TS   pts
  SNES snes
  KSP  ksp
  PetscInt steps,ierr
  PetscReal t
  Vec wVec
  PetscFortranAddr ctx(*) 
  integer(kind=intType) :: ksp_its,temp,nfuncs

  if (myid == 0) then
     print *,'Monitor! Step,t:',steps,t
  end if


  call TSGetSNES(pts,snes,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call SNESGetKSP(snes,ksp,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call KSPGetTolerances(ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ierr)
  ! Set the ksp_atol to snes_atol...this lets the ksp exit early if we
  ! have hit are desired non-linear tolerance
  ksp_atol = snes_atol

  ! Set the maximum iterations ot the subspace size...not really much
  ! point to keep going...might as well do another non-linear iteration
  ksp_max_it = ksp_subspace

  call KSPSetTolerances(ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ierr)

  ! We want to get the number of function evals 
  ! its == 1 AFTER the first iteration. Take this opportunity to set
  ! the preconditioner lag to the value we actually want
  if (steps == 1) then
     ! Reset the value of the preconditionerLag to what we actually want. It
     ! had been set to -1 or -2 depending on if we wanted to recompute
     ! the preconditioner on the first entry or not. 

     call SNESSetLagJacobian(snes, jacobian_lag, ierr)
     call EChk(ierr,__FILE__,__LINE__)

  end if

  if (steps > 0 .or. iterTot0 == 0) then
     call SNESGetNumberFunctionEvals(snes, nfuncs,ierr)
     call EChk(ierr,__FILE__,__LINE__)  

     nFuncs = max(nfuncs,1)
     iterTot = iterTot0 + nFuncs

     call convergenceInfo

     ierr = 0
  end if


end subroutine ts_monitor
