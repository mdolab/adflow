!
!      ******************************************************************
!      *                                                                *
!      * File:          NKsolver.F90                                    *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 11-27-2010                                      *
!      * Last modified: 11-27-2010                                      *
!      *                                                                *
!      ******************************************************************

subroutine NKsolver
  use communication
  use constants
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: snes,dRdw,dRdwPre,ctx,jacobian_lag,&
       snes_stol,snes_max_funcs,nksolversetup,rhoRes0,rhoresstart, &
       snes_rtol,snes_atol,totalR0,itertot0,wVec,rVec,reason,NKSolvedOnce
  use InputIO ! L2conv,l2convrel
  use inputIteration
  use monitor
  use killSignals
  use iteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  integer(kind=intTYpe) :: sns_max_its,ierr,snes_max_its,temp
  real(kind=realType) :: rhoRes,totalRRes,rhoRes1

  ! We are going to have to compute what the tolerances should be
  ! since we are going to be using the same convergence criteria as
  ! SUmb originally uses, that is L2Conv and L2ConvRel. This however,
  ! gets a little trickier, since the NKsolver will always be called
  ! after the RK solver has been run at least once to get a good
  ! starting point. 

  snes_stol = 1e-14
  snes_max_funcs = ncycles-iterTot
  ! Since we're only interested in the maximum function evals, just
  ! set the max its to the same value...it will always punch out on
  ! funcs before iterations in this case

  snes_max_its = ncycles-iterTot

  ! Determine the current level of convergence of the solution

  call getCurrentResidual(rhoRes,totalRRes)

  ! We need to compute two convergences: One coorsponding to L2ConvRel
  ! and one for L2Conv
     
  snes_rtol = (rhoResStart * L2ConvRel)/ rhoRes  ! Target / Current
     
  ! Absolute Tol is the original totalR * L2conv

  snes_atol = totalR0 * L2Conv

  call SNESSetTolerances(snes,snes_atol,snes_rtol,snes_stol,snes_max_its,&
       snes_max_funcs,ierr); call EChk(ierr,__FILE__,__LINE__)

  ! Note: the krylov linear solver options are set in FormJacobian

  ! Form the initial guess from the current w-vector
  call setwVec(wVec)

  ! Solve IT!
  call SNESSolve(snes,PETSC_NULL_OBJECT,wVec,ierr) ! PETSC_NULL_OBJECT
                                                   ! MAY GIVE MEMORY
                                                   ! LEAK!!!!!!

  NKSolvedOnce = .True.

  call EChk(ierr,__FILE__,__LINE__)
  !iterTot = iterTot0

  call SNESGetConvergedReason(snes,reason,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (reason == SNES_CONVERGED_FNORM_ABS .or. &
      reason == SNES_CONVERGED_FNORM_RELATIVE .or. &
      reason == SNES_CONVERGED_PNORM_RELATIVE) then
     routineFailed = .False.
  else
     routineFailed = .True. 
  end if

end subroutine NKsolver


subroutine NKsolverPseudo
  use communication
  use constants
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: pts,snes,dRdw,dRdwPre,ctx,jacobian_lag,&
       snes_stol,snes_max_funcs,nksolversetup,rhoRes0,rhoresstart, &
       snes_rtol,snes_atol,totalR0,itertot0,wVec,rVec,reason,NKSolvedOnce
  use InputIO ! L2conv,l2convrel
  use inputIteration
  use monitor
  use killSignals
  use iteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  integer(kind=intTYpe) :: sns_max_its,ierr,snes_max_its,temp,steps
  real(kind=realType) :: rhoRes,totalRRes,rhoRes1,ftime

  ! We are going to have to compute what the tolerances should be
  ! since we are going to be using the same convergence criteria as
  ! SUmb originally uses, that is L2Conv and L2ConvRel. This however,
  ! gets a little trickier, since the NKsolver will always be called
  ! after the RK solver has been run at least once to get a good
  ! starting point. 

  snes_stol = 1e-14
  snes_max_funcs = ncycles-iterTot
  ! Since we're only interested in the maximum function evals, just
  ! set the max its to the same value...it will always punch out on
  ! funcs before iterations in this case

  snes_max_its = ncycles-iterTot

  ! Determine the current level of convergence of the solution

  call getCurrentResidual(rhoRes,totalRRes)

  ! We need to compute two convergences: One coorsponding to L2ConvRel
  ! and one for L2Conv
     
  snes_rtol = (rhoResStart * L2ConvRel)/ rhoRes  ! Target / Current
     
  ! Absolute Tol is the original totalR * L2conv

  snes_atol = totalR0 * L2Conv
  !call TSgetSNES(pts,snes,ierr)
  
  
  !call SNESSetTolerances(snes,snes_atol,snes_rtol,snes_stol,snes_max_its,&
  !     snes_max_funcs,ierr); 
  !call EChk(ierr,__FILE__,__LINE__)

  ! Note: the krylov linear solver options are set in FormJacobian

  ! Form the initial guess from the current w-vector
  call setwVec(wVec)

  call TSSetup(pts,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Solve IT!
  steps=1000
  call TSStep(pts,steps,ftime,ierr)
  
  NKSolvedOnce = .True.

  call EChk(ierr,__FILE__,__LINE__)
  routineFailed = .False.

end subroutine NKsolverPseudo
