!
!      ******************************************************************
!      *                                                                *
!      * File:          NKsolver_custom.F90                             *
!      * Author:        Gaetan Kenway                                   *
!      * Starting date: 11-27-2010                                      *
!      * Last modified: 11-27-2010                                      *
!      *                                                                *
!      ******************************************************************

subroutine NKsolver_custom
#ifndef USE_NO_PETSC
  use communication
  use constants
  use inputTimeSpectral
  use flowVarRefState
  use ADjointVars , only: nCellsLocal
  use NKSolverVars, only: dRdw,dRdwPre,dRdwPseudo,diagV,ctx,jacobian_lag,&
       nksolversetup,rhoRes0,rhoresstart, &
       totalR0,totalRStart,itertot0,wVec,rVec,deltaW,reason,NKSolvedOnce,&
       ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ksp_subspace,reason

  use InputIO ! L2conv,l2convrel
  use inputIteration
  use monitor
  use killSignals
  use iteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Working Variables
  integer(kind=intType) :: iter,ierr,ksp_iterations
  integer(kind=intType) :: maxNonLinearIts
  real(kind=realType) :: norm,old_norm,rtol_last
  real(kind=realType) :: dt_ref,dt_min,alpha,beta
  ! maxNonLinearIts is (far) larger that necessary. The "iteration"
  ! limit is really set from the maxmimum number of funcEvals
  maxNonLinearIts = ncycles-iterTot
  routineFailed = .True.
  ksp_iterations = 0
  norm = 0.0
  old_norm=0.0
  rtol_last =0.0

  ! Set the inital wVec
  call setwVec(wVec)
  dt_ref = 1000
  dt_min = 100
  alpha = 50.0
  beta = 1.5
  ! Master Non-Linear Loop:
  NonLinearLoop: do iter =1,maxNonLinearIts

     ! Set the diagV Vector
     !call setdiagV(dt_ref)

     ! Determine if if we need to form the Preconditioner:
     if (mod(iter-1,jacobian_lag) == 0) then
        call FormJacobian_custom()

        ! Update dt_ref only when jac is assembled and not on first
        ! iter:
        if (not(iter == 1)) then
           dt_ref = max(alpha*(norm/totalRStart)**(-beta),dt_min)
           !dt_ref = 10
        end if
!         if (myid == 0) then
!            print *,'dt_ref:',dt_ref
!         end if
     else
        ! Else just call assmebly begin/end on dRdW
        call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     ! Always call assembly begin/end on the pseudo shell matrix:
     call MatAssemblyBegin(dRdwPseudo,MAT_FINAL_ASSEMBLY,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     call MatAssemblyEnd(dRdwPseudo,MAT_FINAL_ASSEMBLY,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the BaseVector of the matrix-free matrix:
     call MatMFFDSetBase(dRdW,wVec,PETSC_NULL_OBJECT,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Set the new W, evalulate the residual and set it in PETSc Vec.
     call setW(wVec)
     call computeResidualNK()
     call setRVec(rVec)
 
     ! Compute the norm of rVec for use in EW Criteria
     old_norm = norm
     rtol_last = ksp_rtol
     call VecNorm(rVec,NORM_2,norm,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! The extra 1 is from the above compueResidualNK
     if (iter .ne. 1) then
        iterTot = iterTot + ksp_iterations + 1 
        call convergenceInfo
     end if
     
     ! Check to see if we're converged: We need to check if we've meet
     ! L2Conv or L2ConvRel
     if (norm / totalR0 < L2Conv) then
        routineFailed = .False.
        exit NonLinearLoop
     end if
     
     if (norm / totalRStart < L2ConvRel) then
        routineFailed = .False.
        exit NonLinearLoop
     end if

     ! Check to see if we've done too many function Evals:
     if (iterTot > ncycles) then
        exit NonLinearLoop
     end if

     ! Get the EW Forcing tolerance ksp_rtol
     call getEWTol(iter,norm,old_norm,rtol_last,ksp_rtol)

     ! Set all tolerances for linear solve:
     ksp_atol = 1e-12
     ksp_max_it = ksp_subspace
     call KSPSetTolerances(ksp,ksp_rtol,ksp_atol,ksp_div_tol,ksp_max_it,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Actually do the Linear Krylov Solve
     call KSPSolve(ksp,rVec,deltaW,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get convergence reason:
     call KSPGetConvergedReason(ksp,reason,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     if (myid == 0 .and. reason < 0) then
        !print *,'KSP Reason:',reason
     end if

     ! Blindly set the new wVec according to: w = w - deltaW
     ! We really should do a line search here by rights:
     ! VecAXPY:  y = alpha x + y. 
     call VecAXPY(wVec,-1.0,deltaW,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get the number of iterations to use with Convergence Info
     call KSPGetIterationNumber(ksp,ksp_iterations,ierr)
     call EChk(ierr,__FILE__,__LINE__)
  
  end do NonLinearLoop
     
  ! Not really anything else to do...

  NKSolvedOnce = .True.
#endif
end subroutine NKsolver_custom

subroutine getEWTol(iter,norm,old_norm,rtol_last,rtol)

  use precision
  implicit none

   ! There are the default EW Parameters from PETSc. They seem to work well
   !version:           2
   !rtol_0:  0.300000000000000     
   !rtol_max:  0.900000000000000     
   !gamma:   1.00000000000000     
   !alpha:   1.61803398874989     
   !alpha2:   1.61803398874989     
   !threshold:  0.100000000000000     

  integer(kind=intType) :: iter
  real(kind=realType), intent(in) :: norm,old_norm,rtol_last
  real(kind=realType), intent(out) :: rtol
  real(kind=realType) :: rtol_0,rtol_max,gamma,alpha,alpha2,threshold,stol

  rtol_0    = 0.30
  rtol_max  = 0.9
  gamma     = 1.0
  alpha     = (1+sqrt(5.0))/2.0
  alpha2    = (1+sqrt(5.0))/2.0
  threshold = 0.10

  if (iter == 1) then
     rtol = rtol_0
  else
     ! We use version 2:
     rtol = gamma*(norm/old_norm)**alpha
     stol = gamma*(rtol_last)**alpha

     if (stol > threshold) then
        rtol = max(rtol,stol)
     end if
     
     ! Safeguard: avoid rtol greater than one
     rtol = min(rtol,rtol_max)
  end if
 
end subroutine getEWTol

subroutine setdiagV(dt_pseudo)
#ifndef USE_NO_PETSC
  use flowVarRefState
  use inputTimeSpectral
  use blockPointers
  use NKSolverVars, only: diagV

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input
  real(kind=realType) :: dt_pseudo

  ! Working
  integer(kind=intType) :: nn,sps,i,j,k,ierr
  real(kind=realType) :: vals(nw),dt_loc,L_loc
  
 
  spectralLoop: do sps=1,nTimeIntervalsSpectral
     domainLoop: do nn=1,nDom
        ! Set the pointers to this block.
        call setPointersAdj(nn, 1, sps)

        do k=2,kl
           do j=2,jl
              do i=2,il
                 ! Set the I/dt term in diagV according to CFL_pseudo
!                  L_loc = vol(i,j,k)**(1/3)
!                  dt_loc = (L_loc + 1)/L_loc
!                  vals(:) = -dt_loc*dt_pseudo

                 vals(:) = -1/(dt_pseudo * dtl(i,j,k))
                 call VecSetValuesBlocked(diagV,1,globalCell(i,j,k),vals,&
                      INSERT_VALUES,ierr)
                 call EChk(ierr,__FILE__,__LINE__)
              end do
           end do
        end do
     end do domainLoop
  end do spectralLoop

 
  call VecAssemblyBegin(diagV,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecAssemblyEnd(diagV,ierr)
  call EChk(ierr,__FILE__,__LINE__)
#endif
end subroutine setdiagV


