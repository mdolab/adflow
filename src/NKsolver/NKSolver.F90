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
  use NKSolverVars, only: dRdw,dRdwPre,jacobian_lag,&
       totalR0,totalRStart,wVec,rVec,deltaW,reason,global_ksp,reason,&
       ksp_rtol,ksp_atol,ksp_max_it,ksp_subspace,ksp_div_tol,&
       nksolvedonce,times,petsccomm

  use InputIO ! L2conv,l2convrel
  use inputIteration
  use monitor
  use killSignals
  use iteration
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! PETSc Variables:
  Vec g,work

  ! Working Variables
  integer(kind=intType) :: iter,ierr,ksp_iterations
  integer(kind=intType) :: maxNonLinearIts,nfevals
  real(kind=realType) :: norm,old_norm,rtol_last
  real(kind=realType) :: fnorm,ynorm,gnorm

  ! maxNonLinearIts is (far) larger that necessary. The "iteration"
  ! limit is really set from the maxmimum number of funcEvals
  maxNonLinearIts = ncycles-iterTot
  routineFailed = .True.
  ksp_iterations = 0
  norm = 0.0
  old_norm=0.0
  rtol_last =0.0
  nfevals = 0

  times(10) = 0.0
  times(20) = 0.0

  ! Set the inital wVec
  call setwVec(wVec)

  ! Create the two additional work vectors for the line search:
  call VecDuplicate(wVec,g,ierr); call EChk(ierr,__FILE__,__LINE__)
  call VecDuplicate(wVec,work,ierr);  call EChk(ierr,__FILE__,__LINE__)
 
  ! Evaluate the residual before we start and copy the value into g

  if (petscComm) then
     call setW_ghost(wVec)
     call computeResidualNK2()
  else
     call setW(wVec)
     call computeResidualNK()
  end if
  call setRVec(rVec)

  call vecCopy(rVec,g,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  iterTot = iterTot + 1 ! Add this function evaluation
  
  ! Master Non-Linear Loop:
  NonLinearLoop: do iter= 1,maxNonLinearIts
     
     ! Increment the function evals from the Krylov Iterations and the
     ! line search iterations
     if (iter .ne. 1) then
        iterTot = iterTot + ksp_iterations + nfevals 
        call convergenceInfo
     end if

     ! Use the result from the last line search
     call vecCopy(g,rVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Determine if if we need to form the Preconditioner: 
     if (mod(iter-1,jacobian_lag) == 0) then
        call FormJacobian()
     else
        ! Else just call assmebly begin/end on dRdW
        call MatAssemblyBegin(dRdw,MAT_FINAL_ASSEMBLY,ierr)
        call EChk(ierr,__FILE__,__LINE__)
        call MatAssemblyEnd(dRdw,MAT_FINAL_ASSEMBLY,ierr)
        call EChk(ierr,__FILE__,__LINE__)
     end if

     ! Set the BaseVector of the matrix-free matrix:
     call MatMFFDSetBase(dRdW,wVec,PETSC_NULL_OBJECT,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Compute the norm of rVec for use in EW Criteria
     old_norm = norm
     rtol_last = ksp_rtol
     call VecNorm(rVec,NORM_2,norm,ierr)
     call EChk(ierr,__FILE__,__LINE__)

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
     ! Set absolve tolerance so we don't go past our target:
     ksp_atol = totalR0*L2Conv
     ksp_max_it = min(ksp_subspace,ncycles-iterTot)

     call KSPSetTolerances(global_ksp,ksp_rtol,ksp_atol,ksp_div_tol,&
          ksp_max_it,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Actually do the Linear Krylov Solve
     call KSPSolve(global_ksp,rVec,deltaW,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get convergence reason:
     call KSPGetConvergedReason(global_ksp,reason,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Get the number of iterations to use with Convergence Info
     call KSPGetIterationNumber(global_ksp,ksp_iterations,ierr)
     call EChk(ierr,__FILE__,__LINE__)

     ! Increment iterTot from krylov iterations
     iterTot = iterTot + ksp_iterations

     ! Check to see if we've done too many function Evals:
     if (iterTot >= ncycles) then
        call convergenceInfo
        iterTot = ncycles
        exit NonLinearLoop
     end if

     ! Linesearching:
     if (.True.) then! Check for type of line search:
        call LSCubic(wVec,rVec,g,deltaW,work,fnorm,ynorm,gnorm,nfevals)
     else ! No Linesearch, just accept the new step
        call LSNone(wVec,rVec,g,deltaW,work,nfevals)
     end if

     ! Increment the function evals from the line search iterations
     iterTot = iterTot + nfevals 

     ! Check to see if we've done too many function Evals:
     if (iterTot >= ncycles) then
        call convergenceInfo
        iterTot = ncycles
        exit NonLinearLoop
     end if

     ! Copy the work vector to wVec
     call VecCopy(work,wVec,ierr)
     call EChk(ierr,__FILE__,__LINE__)
     
     ! Print current convergence info
     call convergenceInfo

  end do NonLinearLoop
     
  ! Not really anything else to do...

  NKSolvedOnce = .True.

  ! Destroy the additional two vecs:
  call VecDestroy(g,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  call VecDestroy(work,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  
  print *,'myid, times(10)', myid,times(10)
  print *,'myid, times(20)', myid,times(20)

end subroutine NKsolver

subroutine LSCubic(x,f,g,y,w,fnorm,ynorm,gnorm,nfevals)
  use precision 
  use communication
  use NKSolverVars, only: dRdw,petsccomm

  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"

  ! Input/Output
  Vec x,f,g,y,w
  !x 	- current iterate
  !f 	- residual evaluated at x
  !y 	- search direction
  !w 	- work vector -> On output, new iterate
  !g    - residual evaluated at new iterate y

  real(kind=realType) :: fnorm,gnorm,ynorm
  real(kind=realType) :: alpha
  logical :: flag
  integer(kind=intType) :: nfevals
  !   Note that for line search purposes we work with with the related
  !   minimization problem:
  !      min  z(x):  R^n -> R,
  !   where z(x) = .5 * fnorm*fnorm, and fnorm = || f ||_2.
  !         

  real(kind=realType) :: initslope,lambdaprev,gnormprev,a,b,d,t1,t2,rellength
  real(kind=realType) :: minlambda,lambda,lambdatemp

  integer(kind=intType) :: ierr

  ! Set some defaults:
  alpha		= 1.e-4
  minlambda     = 1.e-12
  nfevals = 0
  flag = .True. 

  ! Compute the two norms we need:
  call VecNorm(y,NORM_2,ynorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecNorm(f,NORM_2,fnorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  call VecMaxPointwiseDivide(y,x,rellength,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  minlambda = minlambda/rellength ! Fix this
  call MatMult(dRdw,y,w,ierr)
  call EChk(ierr,__FILE__,__LINE__)
  nfevals = nfevals + 1

  call VecDot(f,w,initslope,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  if (initslope > 0.0)  then
     initslope = -initslope
  end if

  if (initslope == 0.0) then
     initslope = -1.0
  end if

  call VecWAXPY(w,-1.0,y,x,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute Function:

  if (petscComm) then
     call setW_ghost(w)
     call computeResidualNK2()
  else
     call setW(w)
     call computeResidualNK()
  end if
  call setRVec(g)

  nfevals = nfevals + 1

  call VecNorm(g,NORM_2,gnorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Sufficient reduction 
  if (.5*gnorm*gnorm <= .5*fnorm*fnorm + alpha*initslope) then
     goto 100
  end if

  ! Fit points with quadratic 
  lambda     = 1.0
  lambdatemp = -initslope/(gnorm*gnorm - fnorm*fnorm - 2.0*initslope)
  lambdaprev = lambda
  gnormprev  = gnorm

  if (lambdatemp > .5*lambda) then
     lambdatemp = .5*lambda
  end if

  if (lambdatemp <= .1*lambda) then
     lambda = .1*lambda
  else                 
     lambda = lambdatemp
  end if

  call VecWAXPY(w,-lambda,y,x,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute new function again:
  if (petscComm) then
     call setW_ghost(w)
     call computeResidualNK2()
  else
     call setW(w)
     call computeResidualNK()
  end if
  call setRVec(g)

  nfevals = nfevals + 1

  call VecNorm(g,NORM_2,gnorm,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Sufficient reduction 
  if (.5*gnorm*gnorm <= .5*fnorm*fnorm + lambda*alpha*initslope) then
     goto 100
  end if

  ! Fit points with cubic 
  cubic_loop: do while (.True.) 
    if (lambda <= minlambda) then 
       flag = .False.
       exit cubic_loop
    end if
    t1 = .5*(gnorm*gnorm - fnorm*fnorm) - lambda*initslope
    t2 = .5*(gnormprev*gnormprev  - fnorm*fnorm) - lambdaprev*initslope

    a  = (t1/(lambda*lambda) - t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev)
    b  = (-lambdaprev*t1/(lambda*lambda) + lambda*t2/(lambdaprev*lambdaprev))/(lambda-lambdaprev)
    d  = b*b - 3*a*initslope
    if (d < 0.0) then
       d = 0.0
    end if

    if (a == 0.0) then
       lambdatemp = -initslope/(2.0*b)
    else
       lambdatemp = (-b + sqrt(d))/(3.0*a)
    end if

    lambdaprev = lambda
    gnormprev  = gnorm

    if (lambdatemp > .5*lambda)  then
       lambdatemp = .5*lambda
    end if
    if (lambdatemp <= .1*lambda) then
       lambda = .1*lambda
    else           
       lambda = lambdatemp
    end if

    call  VecWAXPY(w,-lambda,y,x,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Compute new function again:
    if (petscComm) then
       call setW_ghost(w)
       call computeResidualNK2()
    else
       call setW(w)
       call computeResidualNK()
    end if
    call setRVec(g)
    nfevals = nfevals + 1

    call VecNorm(g,NORM_2,gnorm,ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Is reduction enough?
    if (.5*gnorm*gnorm <= .5*fnorm*fnorm + lambda*alpha*initslope) then
       exit cubic_loop
  end if
 end do cubic_loop

100 continue

 ! Optional user-defined check for line search step validity */

end subroutine LSCubic

subroutine LSNone(x,f,g,y,w,nfevals)
  use precision 
  use communication
  use NKSolverVars, only: dRdw
  implicit none
#define PETSC_AVOID_MPIF_H
#include "include/finclude/petsc.h"
  
  ! Input/Output
  Vec x,f,g,y,w
  !x 	- current iterate
  !f 	- residual evaluated at x
  !y 	- search direction
  !w 	- work vector -> On output, new iterate
  !g    - residual evaluated at new iterate y

  integer(kind=intType) :: nfevals
  integer(kind=intType) :: ierr

  ! We just accept the step and compute the new residual at the new iterate
  nfevals = 0
  call VecWAXPY(w,-1.0,y,x,ierr)
  call EChk(ierr,__FILE__,__LINE__)

  ! Compute new function:
  call setW_ghost(w)
  call computeResidualNK2()
  call setRVec(g)
  nfevals = nfevals + 1

end subroutine LSNone


! subroutine setdiagV(dt_pseudo)

!   use flowVarRefState
!   use inputTimeSpectral
!   use blockPointers
!   use NKSolverVars, only: diagV

!   implicit none
! #define PETSC_AVOID_MPIF_H
! #include "include/finclude/petsc.h"

!   ! Input
!   real(kind=realType) :: dt_pseudo

!   ! Working
!   integer(kind=intType) :: nn,sps,i,j,k,ierr
!   real(kind=realType) :: vals(nw),dt_loc,L_loc
  
 
!   spectralLoop: do sps=1,nTimeIntervalsSpectral
!      domainLoop: do nn=1,nDom
!         ! Set the pointers to this block.
!         call setPointersAdj(nn, 1, sps)

!         do k=2,kl
!            do j=2,jl
!               do i=2,il
!                  ! Set the I/dt term in diagV according to CFL_pseudo
! !                  L_loc = vol(i,j,k)**(1/3)
! !                  dt_loc = (L_loc + 1)/L_loc
! !                  vals(:) = -dt_loc*dt_pseudo

!                  vals(:) = -1/(dt_pseudo * dtl(i,j,k))
!                  call VecSetValuesBlocked(diagV,1,globalCell(i,j,k),vals,&
!                       INSERT_VALUES,ierr)
!                  call EChk(ierr,__FILE__,__LINE__)
!               end do
!            end do
!         end do
!      end do domainLoop
!   end do spectralLoop

 
!   call VecAssemblyBegin(diagV,ierr)
!   call EChk(ierr,__FILE__,__LINE__)
!   call VecAssemblyEnd(diagV,ierr)
!   call EChk(ierr,__FILE__,__LINE__)
! end subroutine setdiagV

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
